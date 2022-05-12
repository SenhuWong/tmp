#include "TopologyHolderStrategy.h"
#include "UnstructIntegrator.h"
#include "mpi.h"
#include <iostream>
#include <cmath>
TopologyHolderStrategy::TopologyHolderStrategy(UnstructTopologyHolder *hder,int nequ)
    : d_hder(hder)
{
    d_dim = d_hder->d_dim;
    d_nmesh = d_hder->d_nmesh;
    d_NEQU = nequ;
    commSendBuffer = new double *[d_hder->relatedProcs.size()];
    commRecvBuffer = new double *[d_hder->relatedProcs.size()];
    for (int i = 0; i < d_hder->relatedProcs.size(); i++)
    {
        commSendBuffer[i] = new double[d_hder->sendCellPtr(i, d_nmesh) * d_NEQU];
        commRecvBuffer[i] = new double[d_hder->recvCellPtr(i, d_nmesh) * d_NEQU];
    }
    std::cout<<"Constructor of TopologyHoldersTRATEGy finished\n";
}

TopologyHolderStrategy::~TopologyHolderStrategy()
{
}

void TopologyHolderStrategy::set_fs_variable(double Ma, double AOA, double density, double pressure, double eigenLen)
{
    // Ma is 0.8
    fs_mach = Ma;
    // AOA is 1.25
    fs_AOA = PI * AOA / 180;

    // gamma is 1.4
    // density is 1.225
    // pressure is 101325.0

    fs_density = density;
    fs_pressure = pressure;
    fs_eigen_Length = eigenLen;

    fs_soundSpeed = sqrtf64(Gamma * fs_pressure / fs_density);
    fs_velocity_magnitude = fs_mach * fs_soundSpeed;
    fs_velocity_components[0] = fs_velocity_magnitude * cosf64(fs_AOA);
    fs_velocity_components[1] = fs_velocity_magnitude * sinf64(fs_AOA);

    fs_E = 0.5 * fs_density * std::pow(fs_velocity_magnitude, 2) + fs_pressure / (Gamma - 1);
}

void TopologyHolderStrategy::update_nondim_variable()
{
    // u v
    fsnd_velocity_components[0] = fs_velocity_components[0] / (fs_soundSpeed / sqrtf64(Gamma));
    fsnd_velocity_components[1] = fs_velocity_components[1] / (fs_soundSpeed / sqrtf64(Gamma));
    // p d e
    fsnd_pressure = fs_pressure / fs_pressure;
    fsnd_density = fs_density / fs_density;
    fsnd_soundSpeed = sqrtf64(Gamma);
    fsnd_E = fs_E * Gamma / pow(fs_soundSpeed, 2);

    fs_primVar[0] = fsnd_density;
    fs_primVar[1] = fsnd_pressure;
    fs_primVar[2] = fsnd_velocity_components[0];
    fs_primVar[3] = fsnd_velocity_components[1];
}

void TopologyHolderStrategy::AllCellCommunication(double** value)
{
    // if (cur_proc == 0)
    //     std::cout << "Within block communication";
    for (int i = 0; i < d_hder->d_nmesh; i++)
    {

        auto iter = d_hder->sendBegin(i); // cur_ub.d_to_send.begin();
        for (int j = 0; j < d_hder->relatedProcs.size(); j++)
        {
            for (int k = d_hder->sendCellPtr(j, i); k < d_hder->sendCellPtr(j, i + 1); k++)
            {

                for (int l = 0; l < d_NEQU; l++)
                {

                    commSendBuffer[j][d_NEQU * k + l] = value[i][l+d_NEQU*iter->d_localId];
                }
                iter++;
            }
        }
    }

    std::vector<MPI_Request> recvRequests;
    std::vector<MPI_Request> sendRequests;
    std::vector<MPI_Status> recvStatus;
    recvRequests.resize(d_hder->relatedProcs.size());
    sendRequests.resize(d_hder->relatedProcs.size());
    recvStatus.resize(d_hder->relatedProcs.size());
    // recorder << "-----------------------------------------\n";

    for (int i = 0; i < d_hder->relatedProcs.size(); i++)
    {

        MPI_Irecv(commRecvBuffer[i], d_NEQU * d_hder->recvCellPtr(i, d_nmesh), MPI_DOUBLE, d_hder->relatedProcs[i], 0, MPI_COMM_WORLD, &(recvRequests[i]));
        MPI_Isend(commSendBuffer[i], d_NEQU * d_hder->sendCellPtr(i, d_nmesh), MPI_DOUBLE, d_hder->relatedProcs[i], 0, MPI_COMM_WORLD, &(sendRequests[i]));
    }
    MPI_Waitall(d_hder->relatedProcs.size(), recvRequests.data(), recvStatus.data());

    for (int i = 0; i < d_hder->d_nmesh; i++)
    {
        auto iter = d_hder->recvBegin(i);
        for (int j = 0; j < d_hder->relatedProcs.size(); j++)
        {
            for (int k = d_hder->recvCellPtr(j, i); k < d_hder->recvCellPtr(j, i + 1); k++)
            {
                // recorder << iter->d_localId << '\n';
                for (int l = 0; l < d_NEQU; l++)
                {
                    value[i][l+d_NEQU*iter->d_localId] = commRecvBuffer[j][k * d_NEQU + l];
                    // recorder << commRecvBuffer[j][k * d_NEQU + l] << " ";
                }
                // recorder << '\n';
                iter++;
            }
        }
    }
}

void TopologyHolderStrategy::NearCellCommunication(double** value)
{
    int temp_int = 0;
    for (int i = 0;i < d_hder->d_nmesh; i++)
    {
        for (int j = 0;j< d_hder->relatedProcs.size();j++)
        {
            for(int k = d_hder->sendNearCellPtr(j,i);k<d_hder->sendNearCellPtr(j,i+1);k++)
            {
                auto iter = d_hder->nearSend(i,temp_int++);
                for(int l = 0;l<d_NEQU;l++)
                {
                    commSendBuffer[j][d_NEQU *k +l] = value[i][l+d_NEQU*iter->d_localId];
                }
            }
        }
    }
    std::vector<MPI_Request> recvRequests;
    std::vector<MPI_Request> sendRequests;
    std::vector<MPI_Status> recvStatus;
    recvRequests.resize(d_hder->relatedProcs.size());
    sendRequests.resize(d_hder->relatedProcs.size());
    recvStatus.resize(d_hder->relatedProcs.size());

    for (int i = 0; i < d_hder->relatedProcs.size(); i++)
    {
        MPI_Irecv(commRecvBuffer[i],d_NEQU* d_hder->recvNearCellPtr(i, d_nmesh), MPI_DOUBLE, d_hder->relatedProcs[i], 0, MPI_COMM_WORLD, &(recvRequests[i]));
        MPI_Isend(commSendBuffer[i],d_NEQU* d_hder->sendNearCellPtr(i, d_nmesh), MPI_DOUBLE, d_hder->relatedProcs[i], 0, MPI_COMM_WORLD, &(sendRequests[i]));
    }
    MPI_Waitall(d_hder->relatedProcs.size(), recvRequests.data(), recvStatus.data());
    temp_int = 0;
    for (int i = 0; i < d_hder->d_nmesh; i++)
    {
        for (int j = 0;j < d_hder->relatedProcs.size(); j++)
        {
            for (int k = d_hder->recvNearCellPtr(j,i);k<d_hder->recvNearCellPtr(j,i+1);k++)
            {
                auto iter = d_hder->nearRecv(i,temp_int++);
                for(int l = 0;l<d_NEQU;l++)
                {
                    value[i][l+d_NEQU*iter->d_localId] = commRecvBuffer[j][d_NEQU*k+l];
                }
                
            }
        }
    }

}

void TopologyHolderStrategy::RemoteCellCommunication(double** value)
{
    int temp_int = 0;
    for (int i = 0;i < d_hder->d_nmesh; i++)
    {
        for (int j = 0;j< d_hder->relatedProcs.size();j++)
        {
            for(int k = d_hder->sendRemoteCellPtr(j,i);k<d_hder->sendRemoteCellPtr(j,i+1);k++)
            {
                auto iter = d_hder->remoteSend(i,temp_int++);
                for(int l = 0;l<d_NEQU;l++)
                {
                    commSendBuffer[j][d_NEQU *k +l] = value[i][l+d_NEQU*iter->d_localId];
                }
            }
        }
    }
    std::vector<MPI_Request> recvRequests;
    std::vector<MPI_Request> sendRequests;
    std::vector<MPI_Status> recvStatus;
    recvRequests.resize(d_hder->relatedProcs.size());
    sendRequests.resize(d_hder->relatedProcs.size());
    recvStatus.resize(d_hder->relatedProcs.size());

    for (int i = 0; i < d_hder->relatedProcs.size(); i++)
    {
        MPI_Irecv(commRecvBuffer[i],d_NEQU* d_hder->recvRemoteCellPtr(i, d_nmesh), MPI_DOUBLE, d_hder->relatedProcs[i], 0, MPI_COMM_WORLD, &(recvRequests[i]));
        MPI_Isend(commSendBuffer[i],d_NEQU* d_hder->sendRemoteCellPtr(i, d_nmesh), MPI_DOUBLE, d_hder->relatedProcs[i], 0, MPI_COMM_WORLD, &(sendRequests[i]));
    }
    MPI_Waitall(d_hder->relatedProcs.size(), recvRequests.data(), recvStatus.data());
    temp_int = 0;
    for (int i = 0; i < d_hder->d_nmesh; i++)
    {
        for (int j = 0;j < d_hder->relatedProcs.size(); j++)
        {
            for (int k = d_hder->recvRemoteCellPtr(j,i);k<d_hder->recvRemoteCellPtr(j,i+1);k++)
            {
                auto iter = d_hder->remoteRecv(i,temp_int++);
                for(int l = 0;l<d_NEQU;l++)
                {
                    value[i][l+d_NEQU*iter->d_localId] = commRecvBuffer[j][d_NEQU*k+l];
                }                
            }
        }
    }
}


