#include "TopologyHolderStrategy.h"

#include "UnstructIntegrator.h"
#include "mpi.h"
#include <iostream>
#include <hdf5.h>
#include <algorithm>
#include <vector>
class PressureHolder2D
{
    public:

    GeomElements::vector3d<2,double> center;
    double normal_y;
    double cp;
    PressureHolder2D(double x,double y,double normal,double cp_in)
    {
        center[0] = x;
        center[1] = y;
        normal_y = normal;
        cp = cp_in;
    }
    PressureHolder2D(PressureHolder2D& another)
    {
        center[0] = (another.center)[0];
        center[1] = (another.center)[1];
        normal_y = another.normal_y;
        cp = another.cp;
    }
    PressureHolder2D(const PressureHolder2D& another)
    {
        center = another.center;
        normal_y = another.normal_y;
        cp = another.cp;
    }

};
bool rankWithXY(PressureHolder2D &e1,PressureHolder2D &e2);


void TopologyHolderStrategy::test_unwantedSweep(int** UnwantedInd,int* nUnwantedInd)
{
    testing_flag = new int *[d_hder->d_nmesh];
    
    for (int i = 0; i < d_hder->d_nmesh; i++)
    {
        testing_flag[i] = new int[d_hder->nCells(i)];
        for (int j = 0; j < d_hder->nCells(i); j++)
        {
            testing_flag[i][j] = -1;
        }
    }
    for (int i = 0; i < d_hder->d_nmesh; i++)
    {
        for(int j = 0;j<nUnwantedInd[i];j++)
        {
            testing_flag[i][UnwantedInd[i][j]] = 0;
        }
    }
    for (int i = 0; i < d_hder->d_nmesh; i++)
    {
        writeTiogaFormat("Unwanted0", cur_proc, i, testing_flag[i]);
    }
    int nfringe = 2;

    for (int k = 0;k<nfringe;k++)
    {
        for(int i = 0;i<d_hder->d_nmesh;i++)
        {
            auto& curBlk = d_hder->blk2D[i];
            for(int j = 0;j<d_hder->nCells(i);j++)
            {
                if(testing_flag[i][j] == k)
                {
                    auto& curCell = curBlk.d_localCells[j];
                    for(int l = 0;l<curCell.edge_size();l++)
                    {
                        auto& curEdge = curBlk.d_localEdges[curCell.edgeInd(l)];
                        int lC = curEdge.lCInd();
                        int rC = curEdge.rCInd();
                        if(lC ==j)
                        {
                            if(testing_flag[i][j]==k)
                            {
                                if(rC>=0)
                                {
                                    if(testing_flag[i][rC]==-1)
                                    {
                                        testing_flag[i][rC] = k+1;
                                    }
                                }
                            }
                        }
                        else if(rC ==j)
                        {
                            if(testing_flag[i][j]==k)
                            {
                                if(testing_flag[i][lC]==-1)
                                {
                                    testing_flag[i][lC] = k+1;
                                }
                            }

                        }
                        else
                        {
                            std::cout<<lC<<","<<rC<<","<<j<<","<<curCell.edge_size()<<'\n';
                            std::cout<<curCell.pointInd(0)<<" "<<curCell.pointInd(1)<<'\n';
                            throw std::runtime_error("cell is not on either side of edge");
                        }
                    }
                }
            }
        }
    }
    for (int i = 0; i < d_hder->d_nmesh; i++)
    {
        writeTiogaFormat("Unwanted"+std::to_string(nfringe), cur_proc, i, testing_flag[i]);
    }
    // Allocate data storage
    for(int i = 0;i<d_hder->d_nmesh;i++)
    {
        delete[] testing_flag[i];
    }
    delete[] testing_flag;
}

void TopologyHolderStrategy::test_communication()
{
    std::ofstream recorder;
    std::string localFilename = "recorderTestComm" + std::to_string(cur_proc);
    recorder.open(localFilename);
    testing_flag = new int *[d_hder->d_nmesh];
    int **receiver_flag = new int *[d_hder->d_nmesh];
    for (int i = 0; i < d_hder->d_nmesh; i++)
    {
        receiver_flag[i] = new int[d_hder->nCells(i)];
        testing_flag[i] = new int[d_hder->nCells(i)];
        for (int j = 0; j < d_hder->nCells(i); j++)
        {
            receiver_flag[i][j] = -1;
            testing_flag[i][j] = -1;
        }
    }
    // std::cout << "Sending flag correction\n";
    for (int i = 0; i < d_hder->d_nmesh; i++)
    {
        int temp_count = 0;
        int temp_validCount = 0;
        for (auto iter = d_hder->sendBegin(i); iter != d_hder->sendEnd(i); iter++)
        {
            temp_count++;
            if (iter->d_localId >= 0)
                temp_validCount++;
            testing_flag[i][iter->d_localId] = iter->d_localId % 3 + 1;
            recorder << iter->d_localId << '\n';
            // std::cout<<iter->d_localId<<'\n';
        }
        for (auto iter = d_hder->recvBegin(i); iter != d_hder->recvEnd(i); iter++)
        {
            testing_flag[i][iter->d_localId] = iter->d_layer_level;
            receiver_flag[i][iter->d_localId] = iter->d_layer_level;//-1;
        }
        if (cur_proc == 0)
        {
            std::cout << "Temp count is " << temp_count << " Vliad is " << temp_validCount << '\n';
        }
    }
    for (int i = 0; i < d_hder->d_nmesh; i++)
    {
        writeTiogaFormat("BeforeExchangeTiogaStyle", cur_proc, i, testing_flag[i]);
        writeTiogaFormat("BeforeExchangeReceiverFlag", cur_proc, i, receiver_flag[i]);
    }
    // Allocate data storage
    int **send_buffer = new int *[d_hder->relatedProcs.size()];
    for (int i = 0; i < d_hder->relatedProcs.size(); i++)
    {
        send_buffer[i] = new int[d_hder->sendCellPtr(i, d_nmesh)];
    }
    // std::cout << '\n';
    for (int i = 0; i < d_hder->d_nmesh; i++)
    {
        auto iter = d_hder->sendBegin(i); // cur_ub.d_to_send.begin();
        for (int j = 0; j < d_hder->relatedProcs.size(); j++)
        {
            for (int k = d_hder->sendCellPtr(j, i); k < d_hder->sendCellPtr(j, i + 1); k++)
            {
                send_buffer[j][k] = testing_flag[i][iter->d_localId];
                // std::cout<<iter->d_localId<<'\n';
                if (send_buffer[j][k] == -1)
                {
                    std::cout << "Sender is -1\n";
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
    recorder << "-----------------------------------------\n";
    int **recv_buffer = new int *[d_hder->relatedProcs.size()];
    for (int i = 0; i < d_hder->relatedProcs.size(); i++)
    {
        recv_buffer[i] = new int[d_hder->recvCellPtr(i, d_nmesh)];
        MPI_Irecv(recv_buffer[i], d_hder->recvCellPtr(i, d_nmesh), MPI_INT, d_hder->relatedProcs[i], 0, MPI_COMM_WORLD, &(recvRequests[i]));
        MPI_Isend(send_buffer[i], d_hder->sendCellPtr(i, d_nmesh), MPI_INT, d_hder->relatedProcs[i], 0, MPI_COMM_WORLD, &(sendRequests[i]));
    }
    MPI_Waitall(d_hder->relatedProcs.size(), recvRequests.data(), recvStatus.data());
    for (int i = 0; i < d_hder->d_nmesh; i++)
    {
        auto iter = d_hder->recvBegin(i);
        for (int j = 0; j < d_hder->relatedProcs.size(); j++)
        {
            for (int k = d_hder->recvCellPtr(j, i); k < d_hder->recvCellPtr(j, i + 1); k++)
            {
                if (iter->d_localId < 0)

                    std::cout << iter->d_localId << "\n";
                recorder << iter->d_localId << '\n';
                testing_flag[i][iter->d_localId] = recv_buffer[j][k];
                iter++;
            }
        }
    }
    for (int i = 0; i < d_hder->d_nmesh; i++)
    {
        writeTiogaFormat("AfterExchangeTiogaStyle", cur_proc, i, testing_flag[i]);
    }
    recorder.close();
}

void TopologyHolderStrategy::test_partialcomm()
{
    testing_flag = new int *[d_hder->d_nmesh];
    for (int i = 0; i < d_hder->d_nmesh; i++)
    {
        testing_flag[i] = new int[d_hder->nCells(i)];
        for (int j = 0; j < d_hder->nCells(i); j++)
        {
            testing_flag[i][j] = -1;
        }
    }
    for (int i = 0; i < d_hder->d_nmesh; i++)
    {
        for(int j = 0;j<d_hder->nRemoteSend(i);j++)
        {
            auto iter = d_hder->remoteSend(i,j);
            testing_flag[i][iter->d_localId] = iter->d_localId % 3 + 1;
        }
        for(int j = 0;j<d_hder->nRemoteRecv(i);j++)
        {
            auto iter = d_hder->remoteRecv(i,j);
            testing_flag[i][iter->d_localId] = -2;
        }
    }
    for (int i = 0; i < d_hder->d_nmesh; i++)
    {
        writeTiogaFormat("BeforeExchangeRemoteComm", cur_proc, i, testing_flag[i]);
    }
    // Allocate data storage
    int **send_buffer = new int *[d_hder->relatedProcs.size()];
    for (int i = 0; i < d_hder->relatedProcs.size(); i++)
    {
        send_buffer[i] = new int[d_hder->sendRemoteCellPtr(i, d_nmesh)];
    }
    int temp_int = 0;
    for (int i = 0;i < d_hder->d_nmesh; i++)
    {
        for (int j = 0;j< d_hder->relatedProcs.size();j++)
        {
            for(int k = d_hder->sendRemoteCellPtr(j,i);k<d_hder->sendRemoteCellPtr(j,i+1);k++)
            {
                auto iter = d_hder->remoteSend(i,temp_int++);
                send_buffer[j][k] = testing_flag[i][iter->d_localId];
            }
        }
    }

    std::vector<MPI_Request> recvRequests;
    std::vector<MPI_Request> sendRequests;
    std::vector<MPI_Status> recvStatus;
    recvRequests.resize(d_hder->relatedProcs.size());
    sendRequests.resize(d_hder->relatedProcs.size());
    recvStatus.resize(d_hder->relatedProcs.size());
    int **recv_buffer = new int *[d_hder->relatedProcs.size()];
    for (int i = 0; i < d_hder->relatedProcs.size(); i++)
    {
        recv_buffer[i] = new int[d_hder->recvRemoteCellPtr(i, d_nmesh)];
        MPI_Irecv(recv_buffer[i], d_hder->recvRemoteCellPtr(i, d_nmesh), MPI_INT, d_hder->relatedProcs[i], 0, MPI_COMM_WORLD, &(recvRequests[i]));
        MPI_Isend(send_buffer[i], d_hder->sendRemoteCellPtr(i, d_nmesh), MPI_INT, d_hder->relatedProcs[i], 0, MPI_COMM_WORLD, &(sendRequests[i]));
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
                testing_flag[i][iter->d_localId] = recv_buffer[j][k];
            }
        }
    }
    for (int i = 0; i < d_hder->d_nmesh; i++)
    {
        writeTiogaFormat("AfterExchangeRemoteComm", cur_proc, i, testing_flag[i]);
    }
    for(int i = 0;i<d_hder->relatedProcs.size();i++)
    {
        delete[] send_buffer[i];
        delete[] recv_buffer[i];
    }
    delete[] send_buffer;
    delete[] recv_buffer;
    

    for (int i = 0; i < d_hder->d_nmesh; i++)
    {
        for (int j = 0; j < d_hder->nCells(i); j++)
        {
            testing_flag[i][j] = -1;
        }
    } 
    for (int i = 0; i < d_hder->d_nmesh; i++)
    {
        for(int j = 0;j<d_hder->nNearSend(i);j++)
        {
            auto iter = d_hder->nearSend(i,j);
            testing_flag[i][iter->d_localId] = iter->d_localId % 3 + 1;
        }
        for(int j = 0;j<d_hder->nNearRecv(i);j++)
        {
            auto iter = d_hder->nearRecv(i,j);
            testing_flag[i][iter->d_localId] = -2;
        }
    }
    for (int i = 0; i < d_hder->d_nmesh; i++)
    {
        writeTiogaFormat("BeforeExchangeNearComm", cur_proc, i, testing_flag[i]);
    }
    // Allocate data storage
    send_buffer = new int *[d_hder->relatedProcs.size()];
    for (int i = 0; i < d_hder->relatedProcs.size(); i++)
    {
        send_buffer[i] = new int[d_hder->sendNearCellPtr(i, d_nmesh)];
    }
    temp_int = 0;
    for (int i = 0;i < d_hder->d_nmesh; i++)
    {
        for (int j = 0;j< d_hder->relatedProcs.size();j++)
        {
            for(int k = d_hder->sendNearCellPtr(j,i);k<d_hder->sendNearCellPtr(j,i+1);k++)
            {
                auto iter = d_hder->nearSend(i,temp_int++);
                send_buffer[j][k] = testing_flag[i][iter->d_localId];
            }
        }
    }

    recvRequests.resize(d_hder->relatedProcs.size());
    sendRequests.resize(d_hder->relatedProcs.size());
    recvStatus.resize(d_hder->relatedProcs.size());
    recv_buffer = new int *[d_hder->relatedProcs.size()];
    for (int i = 0; i < d_hder->relatedProcs.size(); i++)
    {
        recv_buffer[i] = new int[d_hder->recvNearCellPtr(i, d_nmesh)];
        MPI_Irecv(recv_buffer[i], d_hder->recvNearCellPtr(i, d_nmesh), MPI_INT, d_hder->relatedProcs[i], 0, MPI_COMM_WORLD, &(recvRequests[i]));
        MPI_Isend(send_buffer[i], d_hder->sendNearCellPtr(i, d_nmesh), MPI_INT, d_hder->relatedProcs[i], 0, MPI_COMM_WORLD, &(sendRequests[i]));
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
                testing_flag[i][iter->d_localId] = recv_buffer[j][k];
            }
        }
    }
    for (int i = 0; i < d_hder->d_nmesh; i++)
    {
        writeTiogaFormat("AfterExchangeNearComm", cur_proc, i, testing_flag[i]);
    }

    for(int i = 0;i<d_hder->relatedProcs.size();i++)
    {
        delete[] send_buffer[i];
        delete[] recv_buffer[i];
    }
    delete[] send_buffer;
    delete[] recv_buffer;

    for(int i = 0;i<d_hder->d_nmesh;i++)
    {
        delete[] testing_flag[i];
    }
    delete[] testing_flag;

}


void TopologyHolderStrategy::writeTiogaFormat(const std::string &filename, int proc, int meshTag, int *icelldata)
{
    std::string local_filename = filename + std::to_string(meshTag) + std::to_string(proc);
    std::ofstream fout;
    fout.open(local_filename);
    if (!fout.is_open())
    {
        std::cout << "Open file failure\n";
        return;
    }
    fout << "TITLE =\"Tioga output\"\n";
    fout << "VARIABLES=\"X\",\"Y\",";
    if (d_dim == 3)
    {
        fout << "\"Z\",";
    }
    fout << "\"IBLANK\"\n";
    int nnodes = d_hder->nPoints(meshTag);
    int ncells = d_hder->nCells(meshTag);

    fout << "ZONE T=\"VOL_MIXED\",N=" << nnodes << " E=" << ncells;
    fout << " ET = QUADRILATERAL, F = FEBLOCK\n";
    fout << "VARLOCATION =  (";
    for (int i = 0; i < d_dim; i++)
    {
        fout << i + 1 << "=NODAL,";
    }
    fout << d_dim + 1 << "=CELLCENTERED)\n";
    for (int j = 0; j < d_dim; j++)
    {
        for (int i = 0; i < nnodes; i++)
        {
            fout << d_hder->blk2D[meshTag].d_localPoints[i][j] << '\n';
        }
    }
    for (int i = 0; i < ncells; i++)
    {
        fout << icelldata[i] << '\n';
    }
    for (int i = 0; i < ncells; i++)
    {
        for (int j = 0; j < d_hder->blk2D[meshTag].d_localCells[i].size(); j++)
        {
            fout << d_hder->blk2D[meshTag].d_localCells[i].pointInd(j) + 1 << '\t';
        }
        for (int j = 0; j < 4 - d_hder->blk2D[meshTag].d_localCells[i].size(); j++)
        {
            fout << d_hder->blk2D[meshTag].d_localCells[i].pointInd(d_hder->blk2D[meshTag].d_localCells[i].size() - 1) + 1 << '\t';
        }
        fout << '\n';
    }
    fout.close();
    return;
}


void TopologyHolderStrategy::writeCellData(const std::string &filename, int proc, int meshTag, double **dcelldata)
{
    
    std::string local_filename = filename + "_"+std::to_string(meshTag+1)+"_" + std::to_string(proc+1);
    std::ofstream fout;
    fout.open(local_filename,std::ios::out);
    
    
    if (!fout.is_open())
    {
        std::cout << "Open file failure\n";
        return;
    }
    
    fout << "TITLE =\"Tioga output\"\n";
    fout << "VARIABLES=\"X\",\"Y\",";
    if(d_dim==3)
    {
        fout << "\"Z\",";
    }
    fout << "\"Density\",\"Pressure\",\"U\",\"V\",\n";
    if(d_dim==3)
    {
        fout << "\"W\",";
    }
    fout << "\"Cp\",\"Ma\"";
    int nnodes = d_hder->nPoints(meshTag);
    int ncells = d_hder->nCells(meshTag);
    fout << "ZONE T=\"VOL_MIXED\",N=" << nnodes << " E=" << ncells;
    fout << " ET = QUADRILATERAL, F = FEBLOCK\n";
    fout << "VARLOCATION =  (";
    
    for (int i = 0; i < d_dim; i++)
    {
        fout << i + 1 << "=NODAL,";
    }
    for (int i = 0; i < d_NEQU + 1; i++)
    {
        fout << d_dim + 1 + i << "=CELLCENTERED,";
    }
    fout << d_dim + d_NEQU + 2 << "=CELLCENTERED)\n";
    for (int j = 0; j < d_dim; j++)
    {
        for (int i = 0; i < nnodes; i++)
        {
            fout << d_hder->blk2D[meshTag].d_localPoints[i][j] << '\n';
        }
    }
    for (int i = d_dim; i < d_dim + d_NEQU; i++)
    {
        for (int j = 0; j < d_hder->nCells(meshTag); j++)
        {
            fout << dcelldata[meshTag][i - d_dim+d_NEQU*j] << '\n';
        }
    }
    double fs_velocityMag = (fsnd_velocity_components[0]*fsnd_velocity_components[0]+fsnd_velocity_components[1]*fsnd_velocity_components[1]);

    for(int j = 0;j <d_hder->nCells(meshTag);j++)
    {
        double Cp = 2*(dcelldata[meshTag][1 + d_NEQU*j] - fsnd_pressure)/fs_velocityMag;
        fout << Cp <<'\n';
    }
    for(int j =0;j < d_hder->nCells(meshTag);j++)
    {
        double c = sqrtf64(dcelldata[meshTag][1+d_NEQU*j]/dcelldata[meshTag][0+d_NEQU*j]);
        double velocityMag = 0;
        for(int k = 0;k<d_dim;k++)
        {
            velocityMag += dcelldata[meshTag][2+k+d_NEQU*j]*dcelldata[meshTag][2+k+d_NEQU*j];
        }
        velocityMag = sqrtf64(velocityMag);
        fout <<velocityMag/c<<'\n';
    }

    for (int i = 0; i < ncells; i++)
    {
        for (int j = 0; j < d_hder->blk2D[meshTag].d_localCells[i].size(); j++)
        {
            fout << d_hder->blk2D[meshTag].d_localCells[i].pointInd(j) + 1 << '\t';
        }
        for (int j = 0; j < 4 - d_hder->blk2D[meshTag].d_localCells[i].size(); j++)
        {
            fout << d_hder->blk2D[meshTag].d_localCells[i].pointInd(d_hder->blk2D[meshTag].d_localCells[i].size() - 1) + 1 << '\t';
        }
        fout << '\n';
    }
    fout.close();
    return;
}

void TopologyHolderStrategy::writeNOdeData(const std::string &filename, int proc, int meshTag, double **dnodedata)
{
    std::string local_filename = filename + "_"+std::to_string(meshTag+1)+"_" + std::to_string(proc+1);
    std::ofstream fout;
    fout.open(local_filename,std::ios::out);
    
    
    if (!fout.is_open())
    {
        std::cout << "Open file failure\n";
        return;
    }
    
    fout << "TITLE =\"Tioga output\"\n";
    fout << "VARIABLES=\"X\",\"Y\",";
    if(d_dim==3)
    {
        fout << "\"Z\",";
    }
    for(int i = 0;i<d_NEQU-1;i++)
    {
        fout << "\"Variable"<<i+1<<"\",";
    }
    fout << "\"Variable\"";
    //fout << "\"Density\",\"Pressure\",\"U\",\"V\",\n";
    //if(d_dim==3)
    //{
        //fout << "\"W\",";
    //}
    // fout << "\"Cp\",\"Ma\"";
    int nnodes = d_hder->nPoints(meshTag);
    int ncells = d_hder->nCells(meshTag);
    fout << "ZONE T=\"VOL_MIXED\",N=" << nnodes << " E=" << ncells;
    fout << " ET = QUADRILATERAL, F = FEBLOCK\n";
    fout << "VARLOCATION =  (";
    
    for (int i = 0; i < d_dim; i++)
    {
        fout << i + 1 << "=NODAL,";
    }
    for (int i = 0; i < d_NEQU - 1; i++)
    {
        fout << d_dim + 1 + i << "=NODAL,";
    }
    fout << d_dim + d_NEQU << "=NODAL)\n";
    for (int j = 0; j < d_dim; j++)
    {
        for (int i = 0; i < nnodes; i++)
        {
            fout << d_hder->blk2D[meshTag].d_localPoints[i][j] << '\n';
        }
    }
    for (int i = 0; i < d_NEQU; i++)
    {
        for (int j = 0; j < d_hder->nPoints(meshTag); j++)
        {
            fout << dnodedata[meshTag][i+d_NEQU*j] << '\n';
        }
    }
    for (int i = 0; i < ncells; i++)
    {
        for (int j = 0; j < d_hder->blk2D[meshTag].d_localCells[i].size(); j++)
        {
            fout << d_hder->blk2D[meshTag].d_localCells[i].pointInd(j) + 1 << '\t';
        }
        for (int j = 0; j < 4 - d_hder->blk2D[meshTag].d_localCells[i].size(); j++)
        {
            fout << d_hder->blk2D[meshTag].d_localCells[i].pointInd(d_hder->blk2D[meshTag].d_localCells[i].size() - 1) + 1 << '\t';
        }
        fout << '\n';
    }
    fout.close();
    return;

}

void TopologyHolderStrategy::writeRawCellData(const std::string &filename, int proc, int meshTag, double **dcelldata)
{
    
    std::string local_filename = filename + "_"+std::to_string(meshTag+1)+"_" + std::to_string(proc+1);
    std::ofstream fout;
    fout.open(local_filename,std::ios::out);
    
    
    if (!fout.is_open())
    {
        std::cout << "Open file failure\n";
        return;
    }
    
    fout << "TITLE =\"Tioga output\"\n";
    fout << "VARIABLES=\"X\",\"Y\",";
    if(d_dim==3)
    {
        fout << "\"Z\",";
    }
    fout << "\"Value\"\n";
    int nnodes = d_hder->nPoints(meshTag);
    int ncells = d_hder->nCells(meshTag);
    fout << "ZONE T=\"VOL_MIXED\",N=" << nnodes << " E=" << ncells;
    fout << " ET = QUADRILATERAL, F = FEBLOCK\n";
    fout << "VARLOCATION =  (";
    
    for (int i = 0; i < d_dim; i++)
    {
        fout << i + 1 << "=NODAL,";
    }
    fout << d_dim + + 1 << "=CELLCENTERED)\n";
    for (int j = 0; j < d_dim; j++)
    {
        for (int i = 0; i < nnodes; i++)
        {
            fout << d_hder->blk2D[meshTag].d_localPoints[i][j] << '\n';
        }
    }

    for (int j = 0; j < d_hder->nCells(meshTag); j++)
    {
        fout << dcelldata[meshTag][j] << '\n';
    }

    for (int i = 0; i < ncells; i++)
    {
        for (int j = 0; j < d_hder->blk2D[meshTag].d_localCells[i].size(); j++)
        {
            fout << d_hder->blk2D[meshTag].d_localCells[i].pointInd(j) + 1 << '\t';
        }
        for (int j = 0; j < 4 - d_hder->blk2D[meshTag].d_localCells[i].size(); j++)
        {
            fout << d_hder->blk2D[meshTag].d_localCells[i].pointInd(d_hder->blk2D[meshTag].d_localCells[i].size() - 1) + 1 << '\t';
        }
        fout << '\n';
    }
    fout.close();
}

void TopologyHolderStrategy::outPutCp(std::string& filename, int mesh_ind, int step)
{
    double** U = getU();
    auto& curBlk = d_hder->blk2D[mesh_ind];
    // int
    double* p_average = new double[d_hder->nPoints(mesh_ind)];
    int* count_average = new int[d_hder->nPoints(mesh_ind)];
    for(int i = 0;i<d_hder->nPoints(mesh_ind);i++)
    {
        count_average[i] = 0;
        p_average[i] = 0;
    }
    for(int i = 0;i<d_hder->nCells(mesh_ind);i++)
    {
        auto& curCell = curBlk.d_localCells[i];

        for(int k = 0;k<curCell.size();k++)
        {
            int nodeInd = curCell.pointInd(k);
            
            ++count_average[nodeInd];
            p_average[nodeInd] += U[mesh_ind][1+d_NEQU*i];

        }
    }
    double fs_rhoV2 = fsnd_density*(fsnd_velocity_components[0]*fsnd_velocity_components[0]+fsnd_velocity_components[1]*fsnd_velocity_components[1]);
    
    for(int i = 0;i<d_hder->nPoints(mesh_ind);i++)
    {
        p_average[i] = p_average[i]/count_average[i];
        // std::cout<<count_average[i]<<'\n';
        p_average[i] = (p_average[i] - fsnd_pressure)*2/fs_rhoV2;
    }

    std::set<int,std::less<int>>* wallNodeInds = new std::set<int,std::less<int>>[1];
    

    for(int k = 0;k<d_hder->nEdges(mesh_ind);k++)
    {
        auto& curEdge = curBlk.d_localEdges[k];
        int rC = curEdge.rCInd();
        if(rC==GeomElements::edge3d<2>::BoundaryType::WALL)
        {
            for(int j = 0;j<curEdge.size();j++)
            {
                
                int ptId = curEdge.pointInd(j);
                if(curEdge.normal_vector()[1]>0)//Upper half
                {
                    //std::cout<<"Proc "<<cur_proc<<" found one lower half"<<'\n';
                    count_average[ptId] = -std::abs(count_average[ptId]);
                }
                wallNodeInds->insert(ptId);
            }
        }
    }
    int nlocalWall = wallNodeInds->size();
    //Xloc,Yloc,Count,
    double* realData;
    if(nlocalWall>0)
    {
        realData = new double[nlocalWall*4];
    }
    int m = 0;
    for(auto iter = wallNodeInds->begin();iter!=wallNodeInds->end();iter++)
    {
        int nodeId = *iter;
        auto curPoint = curBlk.d_localPoints[nodeId].d_pos;
        realData[m++] = curPoint[0];
        realData[m++] = curPoint[1];
        realData[m++] = count_average[nodeId];
        realData[m++] = p_average[nodeId];
    }

    // std::ofstream fout0;
    // fout0.open("Cp0"+std::to_string(cur_proc));
    // for(int i = 0;i<m/4;i++)
    // {
    //     fout0<<realData[4*i]<<"\t"<<realData[4*i+1]<<"\t"<<realData[4*i+2]<<"\t"<<realData[4*i+3]<<'\n';
    // }
    // fout0.close();

    int* nEachLocalWall = new int[num_proc+1];

    MPI_Allgather(&nlocalWall,1,MPI_INT,nEachLocalWall+1,1,MPI_INT,MPI_COMM_WORLD);

    nEachLocalWall[0] = 0;
    for(int i = 0;i<num_proc;i++)
    {
        nEachLocalWall[i+1] = nEachLocalWall[i] + nEachLocalWall[i+1];
    }
    //I'm using HDF5 Parallel IO feature for this shit.
    
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    



    hsize_t dims1[2];

    hid_t acc_tpl = H5Pcreate(H5P_FILE_ACCESS);

    herr_t status = H5Pset_fapl_mpio(acc_tpl,comm,info);
    
    hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC,H5P_DEFAULT,acc_tpl);

    status = H5Pclose(acc_tpl);
    
    dims1[0] = nEachLocalWall[num_proc];
    dims1[1] = 4;

    hid_t dataspace_id = H5Screate_simple(2,dims1,NULL);

    hid_t dataset_id = H5Dcreate(file_id,("Cp"+std::to_string(step)).c_str(),H5T_NATIVE_DOUBLE,dataspace_id,
                                H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

    //Now each Proc has its own slice to fill in.
    
    if(nEachLocalWall[cur_proc+1]-nEachLocalWall[cur_proc]!=0)
    {
        hsize_t count[2] = {nlocalWall,4};
        hsize_t offset[2] = {nEachLocalWall[cur_proc],0};
        hsize_t stride[2] = {1,1};
        hsize_t block[2] = {1,1};
        //Need to write in

        hid_t memspace_id = H5Screate_simple(2,count,NULL);

        status = H5Sselect_hyperslab(dataspace_id,H5S_SELECT_SET, offset, stride, count, block);

        status = H5Dwrite(dataset_id,H5T_NATIVE_DOUBLE,memspace_id,dataspace_id,H5P_DEFAULT,realData);

        status = H5Sclose(memspace_id);
    }
    status = H5Sclose(dataspace_id);
    status = H5Dclose(dataset_id);
    status = H5Fclose(file_id);
    MPI_Barrier(MPI_COMM_WORLD);
    if(cur_proc==0)
    {
        file_id = H5Fopen(filename.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
        dataset_id = H5Dopen(file_id,("Cp"+std::to_string(step)).c_str(),H5P_DEFAULT);

        dataspace_id =H5Dget_space(dataset_id);
        double* buffer = new double[nEachLocalWall[num_proc]*4];
        status = H5Dread(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,dataspace_id,H5P_DEFAULT,buffer);
        std::vector<PressureHolder2D> holderVector;
        for(int i = 0;i<nEachLocalWall[num_proc];i++)
        {
            holderVector.push_back(PressureHolder2D(buffer[4*i],buffer[4*i+1],buffer[4*i+2],buffer[4*i+3]));   
        }
        std::sort(holderVector.begin(),holderVector.end(),rankWithXY);
        //Remove the redundant ones.
        bool TheSame = true;
        int cur,next;
        int theSameCount = 0;
        m = 0;
        for(int i = 0;i<holderVector.size();i++)
        {
            cur = i;
            //Add cur to the buffer
            buffer[m++] = holderVector[cur].center[0];
            buffer[m++] = holderVector[cur].center[1];
            buffer[m++] = holderVector[cur].normal_y;
            buffer[m++] = holderVector[cur].cp;
            next = (i+1)%(holderVector.size());
            //Compare their x_loc,If the same 
            TheSame=true;
            while(TheSame)
            {
                if(fabsf64(holderVector[cur].center[0]-holderVector[next].center[0])>0.0001)
                {
                    TheSame =false;
                }
                else
                {
                    if(fabsf64(holderVector[next].normal_y)>fabsf64(holderVector[cur].normal_y))
                    {
                        buffer[m-1] = holderVector[next].cp;
                        buffer[m-2] = holderVector[next].normal_y;
                        cur = next;
                        next = next+1;
                    }
                    else
                    {
                        next = next+1;
                    }
                    i++;
                }
            }
        }
        int nCp = m/4;
        std::ofstream fout;
        fout.open("Cp_"+std::to_string(step));
        fout<<"Title = \"Dussin_exp\"\n";
        fout<<"VARIABLES = \"X\",\"Cp\"\n";
        fout<<"ZONE T=\"Dussin_exp\", I = "<<nCp+1<<", F=POINT\n";

        for(int i = 0;i<nCp;i++)
        {
            fout<<buffer[4*i]<<"\t"<<buffer[4*i+3]<<'\n';//<<"\t"<<buffer[4*i+2]<<"\t"<<buffer[4*i+3]<<'\n';
        }
        fout<<buffer[4*0]<<"\t"<<buffer[4*0+3]<<'\n';
        fout.close();

        

        for(int i = 0;i<nEachLocalWall[num_proc];i++)
        {
            buffer[4*i] = holderVector[i].center[0];
            buffer[4*i+1] = holderVector[i].center[1];
            buffer[4*i+2] = holderVector[i].normal_y;
            buffer[4*i+3] = holderVector[i].cp;
        }
        
        status = H5Dwrite(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,dataspace_id,H5P_DEFAULT,buffer);

        status = H5Sclose(dataspace_id);
        status = H5Dclose(dataset_id);
        status = H5Fclose(file_id);
        delete[] buffer;
    }
    delete[] wallNodeInds;
    delete[] p_average;
    delete[] count_average;
    if(nlocalWall>0)
    {
        delete[] realData;
        
    }
    delete[] nEachLocalWall;
}