#include "Euler2D.h"
#include "HLLCFluxStrategy.h"
#include "UnstructIntegrator.h"
#include "Vankatakrishnan_Limiter.h"
double computeEdgeNormal(int dim, double xv[4][3], int nvert, double xnorm[3], bool counter_clock = false);
#include <math.h>
#include "mpi.h"
#include <hdf5.h>
#include <algorithm>
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
//Increasing ranking: e1<e2
bool rankWithXY(PressureHolder2D &e1,PressureHolder2D &e2)
{
    //If e1's normaly is negative(downward, upper face) and e2's normaly is positive (upward lower face)
    //return false means that upper face is e1 > e2
    if(e1.normal_y<0 and e2.normal_y>=0) return false;
    else if(e1.normal_y>=0 and e2.normal_y < 0) return true;
    else
    {
        if(e1.normal_y>0)//Both upper
        {
            if(e1.center[0] < e2.center[0])
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        else
        {
            if(e1.center[0] < e2.center[0])
            {
                return false;
            }
            else
            {
                return true;
            }
        }
    }
}

Euler2D::Euler2D(UnstructTopologyHolder *hder)
    : TopologyHolderStrategy(hder)
{
    // Communication set
    cur_proc = d_hder->curproc;
    num_proc = d_hder->numproc;
    // Topology holder set
    d_limiter = new Vankatakrishnan_Limiter(d_hder, this);
#ifdef DEBUG
    d_limiter->setComm(cur_proc, num_proc);

#endif // DEBUG
    d_fluxComputer = new HLLCFluxStrategy(d_hder, this);
    d_fluxComputer->setComm(num_proc, cur_proc);
    d_nmesh = d_hder->d_nmesh;
    // Allocate Data Storage for this problem
    registerTopology();
    // This should be done with an input_db
    set_fs_variable(0.75, 3, 1.225, 101325, 1);
    //
    update_nondim_variable();
}

void Euler2D::registerTopology()
{
    // std::cout<<"Topology registered\n";
    // Initialize Conservative Variables
    U = new double *[d_nmesh];
    // Update Edge Variable
    U_edge = new double **[d_nmesh];
    // With Edge Variables compute gradient
    gradU = new GeomElements::vector3d<2, double> **[d_nmesh];
    Spectrum_cell = new double *[d_nmesh];
    // Compute left and right value of Conservative Variables
    UL = new double **[d_nmesh];
    UR = new double **[d_nmesh];
    Spectrum_edge = new double *[d_nmesh];
    // Compute flux and fluxSum
    Flux_edge = new double **[d_nmesh];
    Residual = new double *[d_nmesh];

    // Compute Local dt
    // dt = new double **[d_nmesh];

    // Allocate space for sndBuffer
    // One for each communication requiring proc
    commSendBuffer = new double *[d_hder->relatedProcs.size()];
    commRecvBuffer = new double *[d_hder->relatedProcs.size()];
    for (int i = 0; i < d_hder->relatedProcs.size(); i++)
    {
        commSendBuffer[i] = new double[d_hder->sendCellPtr(i, d_nmesh) * d_NEQU];
        commRecvBuffer[i] = new double[d_hder->recvCellPtr(i, d_nmesh) * d_NEQU];
    }
    // normal_vector = new double *[d_nmesh];
    for (int i = 0; i < d_nmesh; i++)
    {
        U[i] = new double [d_NEQU*d_hder->nCells(i)];
        Residual[i] = new double [d_NEQU*d_hder->nCells(i)];
        gradU[i] = new GeomElements::vector3d<2, double> *[d_NEQU];
        Spectrum_cell[i] = new double[d_hder->nCells(i)];

        U_edge[i] = new double *[d_NEQU];
        UL[i] = new double *[d_NEQU];
        UR[i] = new double *[d_NEQU];
        Spectrum_edge[i] = new double[d_hder->nEdges(i)];

        Flux_edge[i] = new double *[d_NEQU];
        for (int j = 0; j < d_NEQU; j++)
        {
            // U[i][j] = new double[d_hder->nCells(i)];
            // Residual[i][j] = new double[d_hder->nCells(i)];
            gradU[i][j] = new GeomElements::vector3d<2, double>[d_hder->nCells(i)];
            U_edge[i][j] = new double[d_hder->nEdges(i)];
            UL[i][j] = new double[d_hder->nEdges(i)];
            UR[i][j] = new double[d_hder->nEdges(i)];

            Flux_edge[i][j] = new double[d_hder->nEdges(i)];
        }
    }
}

void Euler2D::set_fs_variable(double Ma, double AOA, double density, double pressure, double eigenLen)
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

void Euler2D::update_nondim_variable()
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

void Euler2D::test_unwantedSweep(int** UnwantedInd,int* nUnwantedInd)
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

void Euler2D::test_communication()
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
        // auto &cur_ub = d_hder->blk2D[i];
        // std::cout<<"size of send cell is "<<d_hder->blk2D[i].d_to_recv.size()<<"on mesh "<<i<<" proc "<<cur_proc<<'\n';
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
        // auto& cur_ub = d_hder->blk2D[i];
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
        // auto& cur_ub = d_hder->blk2D[i];
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

void Euler2D::test_partialcomm()
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
    // std::cout << "Sending flag correction\n";
    for (int i = 0; i < d_hder->d_nmesh; i++)
    {
        // auto &cur_ub = d_hder->blk2D[i];
        // std::cout<<"size of send cell is "<<d_hder->blk2D[i].d_to_recv.size()<<"on mesh "<<i<<" proc "<<cur_proc<<'\n';
        // int temp_count = 0;
        // int temp_validCount = 0;
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


void Euler2D::writeTiogaFormat(const std::string &filename, int proc, int meshTag, int *icelldata)
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

// Initialize Data to t0
void Euler2D::initializeData()
{
    if (cur_proc == 0)
        std::cout << "Initializing Data\n";
    for (int i = 0; i < d_nmesh; i++)
    {
        for (int k = 0; k < d_hder->nCells(i); k++)
        {
            for (int j = 0; j < d_NEQU; j++)
            {
                U[i][j+d_NEQU*k] = fs_primVar[j];
            }
        }
    }
    for (int k = 0; k < d_NEQU; k++)
    {
        std::cout << fs_primVar[k] << '\n';
    }
    // writeCellData("RightAfterInit", cur_proc, 0, W);
}

void Euler2D::ReconstructionInitial()
{
    // Grad should be set to 0
    for (int i = 0; i < d_nmesh; i++)
    {

        for (int j = 0; j < d_NEQU; j++)
        {
            for (int k = 0; k < d_hder->nCells(i); k++)
            {
                gradU[i][j][k].reset();
            }
        }
    }
    // Spectrum of cell should be set to 0
    for (int i = 0; i < d_nmesh; i++)
    {
        for (int k = 0; k < d_hder->nCells(i); k++)
        {
            Spectrum_cell[i][k] = 0;
            for (int j = 0; j < d_NEQU; j++)
            {
                Residual[i][j+d_NEQU*k] = 0;
            }
        }
    }
    // Limiter related should be set
    d_limiter->init();
}

// Boundary and normal edge's conservative update
void Euler2D::updateEdgeValues()
{
    // if (cur_proc == 0)
    //     std::cout << "updateEdgeVariables\n";
    GeomElements::vector3d<2, double> Leftvelocity;
    GeomElements::vector3d<2, double> Rightvelocity;

    GeomElements::vector3d<2, double> normalVec;

    double normalVelocityComponentLeft;
    GeomElements::vector3d<2, double> tangenVelocityLeft;

    double normalVelocityComponentRight;
    GeomElements::vector3d<2, double> tangenVelocityRight;

    for (int i = 0; i < d_nmesh; i++)
    {
        auto &curBlk = d_hder->blk2D[i];
        for (int k = 0; k < curBlk.d_nEs; k++)
        {
            auto &curEdge = curBlk.d_localEdges[k];
            int lC = curEdge.lCInd();
            int rC = curEdge.rCInd();

            if (rC >= 0) // Normal
            {
                for (int j = 0; j < d_NEQU; j++)
                {
                    U_edge[i][j][k] = (U[i][j+d_NEQU*lC] + U[i][j+d_NEQU*rC]) / 2;
                }
            }
            else if (rC == GeomElements::edge3d<2>::BoundaryType::WALL) // Wall
            {
                Leftvelocity[0] = U[i][2+d_NEQU*lC];
                Leftvelocity[1] = U[i][3+d_NEQU*lC];
                normalVec = curEdge.normal_vector();
                normalVelocityComponentLeft = Leftvelocity.dot_product(normalVec);
                tangenVelocityLeft = Leftvelocity - normalVec * normalVelocityComponentLeft;
                if (tangenVelocityLeft.dot_product(normalVec) > 1.0e-11)
                {
                    throw std::runtime_error("Something wrong with Vector computation\n");
                }
                U_edge[i][0][k] = U[i][0+d_NEQU*lC];
                U_edge[i][1][k] = U[i][1+d_NEQU*lC];
                U_edge[i][2][k] = tangenVelocityLeft[0];
                U_edge[i][3][k] = tangenVelocityLeft[1];
            }
            else if (rC == GeomElements::edge3d<2>::BoundaryType::FARFIELD) // Far field
            {
                // This is the derivation
                Leftvelocity[0] = U[i][2+d_NEQU*lC];
                Leftvelocity[1] = U[i][3+d_NEQU*lC];

                normalVec = curEdge.normal_vector();
                normalVelocityComponentLeft = Leftvelocity.dot_product(normalVec);
                tangenVelocityLeft = Leftvelocity - normalVec * normalVelocityComponentLeft;

                Rightvelocity[0] = fsnd_velocity_components[0];
                Rightvelocity[1] = fsnd_velocity_components[1];

                normalVec = curEdge.normal_vector();
                normalVelocityComponentRight = Rightvelocity.dot_product(normalVec);
                tangenVelocityRight = (Rightvelocity - normalVec * normalVelocityComponentRight);

                GeomElements::vector3d<2, double> velo = Leftvelocity;
                double rho_left = U[i][0+d_NEQU*lC];
                double p_left = U[i][1+d_NEQU*lC];

                if (abs(normalVelocityComponentRight) < fsnd_soundSpeed) // Sub sonic
                {
                    double Vnb = (normalVelocityComponentLeft + normalVelocityComponentRight + 2 * (sqrt(Gamma * p_left / rho_left) - fsnd_soundSpeed) / (Gamma - 1)) / 2;
                    double Cb = (normalVelocityComponentLeft - normalVelocityComponentRight + 2 * (fsnd_soundSpeed + sqrt(Gamma * p_left / rho_left)) / (Gamma - 1)) * (Gamma - 1) / 4;
                    if (normalVelocityComponentRight >= 0) // Out
                    {
                        // double Vnb = (normalVelocityComponentLeft + normalVelocityComponentRight + 2 * (sqrt(Gamma * p_left / rho_left) - fsnd_soundSpeed) / (Gamma - 1)) / 2;
                        // double Cb = (normalVelocityComponentLeft - normalVelocityComponentRight + 2 * (sqrt(Gamma * p_left / rho_left) + fsnd_soundSpeed) / (Gamma - 1)) * (Gamma - 1) / 4;
                        //   s should equal to Soo
                        double Sb = (p_left / pow(rho_left, Gamma));
                        GeomElements::vector3d<2, double> EdgeVelocity = tangenVelocityLeft + normalVec * Vnb;
                        U_edge[i][0][k] = pow((Cb * Cb / (Gamma * Sb)), (1.0 / (Gamma - 1)));
                        U_edge[i][1][k] = U_edge[i][0][k] * Cb * Cb / Gamma;
                        U_edge[i][2][k] = EdgeVelocity[0];
                        U_edge[i][3][k] = EdgeVelocity[1];
                    }
                    else
                    {
                        // double Vnb = (normalVelocityComponentLeft + normalVelocityComponentRight + 2 * (fsnd_soundSpeed - sqrt(Gamma * p_left / rho_left)) / (Gamma - 1)) / 2;
                        // double Cb = (normalVelocityComponentRight - normalVelocityComponentLeft + 2 * (sqrt(Gamma * p_left / rho_left) + fsnd_soundSpeed) / (Gamma - 1)) * (Gamma - 1) / 4;
                        double Sb = (fsnd_pressure / pow(fsnd_density, Gamma));
                        GeomElements::vector3d<2, double> Edgevelocity = tangenVelocityRight + normalVec * Vnb;
                        U_edge[i][0][k] = pow((Cb * Cb / (Gamma * Sb)), (1.0 / (Gamma - 1)));
                        U_edge[i][1][k] = U_edge[i][0][k] * Cb * Cb / Gamma;
                        U_edge[i][2][k] = Edgevelocity[0];
                        U_edge[i][3][k] = Edgevelocity[1];
                    }
                }
                else
                {
                    if (normalVelocityComponentRight >= 0)
                    {
                        U_edge[i][0][k] = U[i][0+d_NEQU*lC];
                        U_edge[i][1][k] = U[i][1+d_NEQU*lC];
                        U_edge[i][2][k] = U[i][2+d_NEQU*lC];
                        U_edge[i][3][k] = U[i][3+d_NEQU*lC];
                    }
                    else
                    {
                        U_edge[i][0][k] = fsnd_density;
                        U_edge[i][1][k] = fsnd_pressure;
                        U_edge[i][2][k] = fsnd_velocity_components[0];
                        U_edge[i][3][k] = fsnd_velocity_components[1];
                    }
                }
            }
            else
            {
                throw std::runtime_error("Undefined Boundary Type\n");
            }
        }
    }
}
// Gradient computation for each cell and LeftRight Value Computation for each edge
// Since there ain't no cell to edge ptr, we have to make two loops

void Euler2D::solveGradient()
{
    for (int i = 0; i < d_nmesh; i++)
    {
        auto &curBlk = d_hder->blk2D[i];
        for (int k = 0; k < d_hder->nEdges(i); k++)
        {
            auto &curEdge = curBlk.d_localEdges[k];
            int lC = curEdge.lCInd();
            int rC = curEdge.rCInd();
            for (int j = 0; j < d_NEQU; j++)
            {
                gradU[i][j][lC] = gradU[i][j][lC] + curEdge.normal_vector() * (curEdge.area() * U_edge[i][j][k]);
            }
            if (rC >= 0)
            {
                for (int j = 0; j < d_NEQU; j++)
                {
                    gradU[i][j][rC] = gradU[i][j][rC] - curEdge.normal_vector() * (curEdge.area() * U_edge[i][j][k]);
                }
            }
        }
        for (int k = 0; k < d_hder->nCells(i); k++)
        {
            auto &curCell = curBlk.d_localCells[k];
            for (int j = 0; j < d_NEQU; j++)
            {
                gradU[i][j][k] = gradU[i][j][k] / curCell.volume();
            }
        }
    }
}

void Euler2D::solveSpectralRadius()
{
    GeomElements::vector3d<2, double> velocity_edge;
    for (int i = 0; i < d_nmesh; i++)
    {
        auto &curBlk = d_hder->blk2D[i];
        for (int k = 0; k < d_hder->nEdges(i); k++)
        {
            auto &curEdge = curBlk.d_localEdges[k];
            velocity_edge[0] = U_edge[i][2][k];
            velocity_edge[1] = U_edge[i][3][k];
            double c = sqrt(Gamma * U_edge[i][1][k] / U_edge[i][0][k]);
            double Vn = velocity_edge.dot_product(curEdge.normal_vector());
            int lC = curEdge.lCInd();
            int rC = curEdge.rCInd();
            Spectrum_edge[i][k] = (std::abs<double>(Vn) + c) * curEdge.area();
            Spectrum_cell[i][lC] += Spectrum_edge[i][k];
            if (rC >= 0)
            {
                Spectrum_cell[i][rC] += Spectrum_edge[i][k];
            }
        }
    }
}

void Euler2D::updateEdgeLeftRightValues()
{
    // if (cur_proc == 0)
    //     std::cout << "Compute Gradient with left and right\n";
    //#define GRAD_RECORDER
    const double EPSILON = 1.0e-5;


    // Compute Limiter for each cell
    // This should be done by a Limiter class
    d_limiter->computeLimiter();
    // d_limiter->write_Lmitera();
    bool hasTrouble = false;
    // Compute left and right value for each Edge
    for (int i = 0; i < d_nmesh; i++)
    {
        auto &curBlk = d_hder->blk2D[i];
        for (int k = 0; k < curBlk.d_nEs; k++)
        {

            auto &curEdge = curBlk.d_localEdges[k];
            int lc = curEdge.lCInd();
            int rc = curEdge.rCInd();

            if (lc >= 0 and rc >= 0)
            {
                for (int j = 0; j < d_NEQU; j++)
                {
                    UL[i][j][k] = UL[i][j][k] * d_limiter->getLimiter(i, j, lc) + U[i][j+d_NEQU*lc];
                }
                for (int j = 0; j < d_NEQU; j++)
                {
                    UR[i][j][k] = UR[i][j][k] * d_limiter->getLimiter(i, j, rc) + U[i][j+d_NEQU*rc];
                }
            }
            else
            {
                if (lc < 0)
                {
                    throw std::runtime_error("Left Cell can't be NULL\n");
                }
                if (rc == GeomElements::edge3d<2>::BoundaryType::WALL) // Wall, set to be the cell variables
                {
                    for (int j = 0; j < d_NEQU; j++)
                    {
                        UL[i][j][k] = UR[i][j][k] = U[i][j+d_NEQU*lc];
                    }
                }
                else if (rc == GeomElements::edge3d<2>::BoundaryType::FARFIELD) // Far_field, set to be the face variables
                {
                    for (int j = 0; j < d_NEQU; j++)
                    {
                        UL[i][j][k] = UR[i][j][k] = U_edge[i][j][k];
                    }
                }
                else
                {
                    std::cout << lc << "," << rc << '\n';
                    std::cout << "SOmething wrong with left and right cell\n";
                    std::cin.get();
                }
            }
        }
    }
    if (hasTrouble)
    {
        throw std::runtime_error("Limiter Error");
    }
}
// Update flux with the corresponding left and right value
void Euler2D::updateFlux()
{
    // if (cur_proc == 0)
    //     std::cout << "Update flux\n";
    d_fluxComputer->computeFlux();
    // First clear FLuxsum from last iteration
    for (int i = 0; i < d_nmesh; i++)
    {
        for (int j = 0; j < d_NEQU; j++)
        {
            for (int k = 0; k < d_hder->nCells(i); k++)
            {
                Residual[i][j+d_NEQU*k] = 0;
            }
        }
    }
    // std::ofstream cellRecorder;
    // cellRecorder.open("RecordFluxCell" + std::to_string(cur_proc));
    // std::ofstream recorder;
    // recorder.open("RecordFluxEdge" + std::to_string(cur_proc));
    // Then update the Residual
    for (int i = 0; i < d_nmesh; i++)
    {
        auto &curBlk = d_hder->blk2D[i];

        for (int k = 0; k < d_hder->nEdges(i); k++)
        {

            auto &curEdge = curBlk.d_localEdges[k];
            int lc = curEdge.lCInd();
            int rc = curEdge.rCInd();
            // if (i == 0)
            // {
            //     recorder << k << " th edge with left and right being " << lc << " " << rc << '\n';
            //     for (int j = 0; j < d_NEQU; j++)
            //     {
            //         recorder << Fi[i][j][k] << " ";
            //     }
            //     recorder << '\n';
            // }

            if (lc >= 0) // Left plus right minus
            {
                for (int j = 0; j < d_NEQU; j++)
                {
                    Residual[i][j+d_NEQU*lc] += Flux_edge[i][j][k];
                }
            }
            else
            {
                std::cout << "left cell ind can't be negative\n";
                std::cin.get();
            }
            if (rc >= 0)
            {
                for (int j = 0; j < d_NEQU; j++)
                {
                    Residual[i][j+d_NEQU*rc] -= Flux_edge[i][j][k];
                }
            }
        }
    }
}

void Euler2D::withinBlockCommunication()
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

                    commSendBuffer[j][d_NEQU * k + l] = U[i][l+d_NEQU*iter->d_localId];
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
        // auto& cur_ub = d_hder->blk2D[i];
        auto iter = d_hder->recvBegin(i);
        for (int j = 0; j < d_hder->relatedProcs.size(); j++)
        {
            for (int k = d_hder->recvCellPtr(j, i); k < d_hder->recvCellPtr(j, i + 1); k++)
            {
                // recorder << iter->d_localId << '\n';
                for (int l = 0; l < d_NEQU; l++)
                {
                    U[i][l+d_NEQU*iter->d_localId] = commRecvBuffer[j][k * d_NEQU + l];
                    // recorder << commRecvBuffer[j][k * d_NEQU + l] << " ";
                }
                // recorder << '\n';
                iter++;
            }
        }
    }
}

void Euler2D::writeCellData(const std::string &filename, int proc, int meshTag, double **dcelldata)
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
    fout << "\"Density\",\"MomentX\",\"MomentY\",\"Energy\",\"Pressure\"\n";
    int nnodes = d_hder->nPoints(meshTag);
    int ncells = d_hder->nCells(meshTag);
    fout << "ZONE T=\"VOL_MIXED\",N=" << nnodes << " E=" << ncells;
    fout << " ET = QUADRILATERAL, F = FEBLOCK\n";
    fout << "VARLOCATION =  (";
    
    for (int i = 0; i < d_dim; i++)
    {
        fout << i + 1 << "=NODAL,";
    }
    for (int i = 0; i < d_NEQU; i++)
    {
        fout << d_dim + 1 + i << "=CELLCENTERED,";
    }
    fout << d_dim + d_NEQU + 1 << "=CELLCENTERED)\n";
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
    for (int j = 0; j < d_hder->nCells(meshTag); j++)
    {
        double VM = 0;
        for (int l = 0; l < d_dim; l++)
        {
            VM += powf64(dcelldata[meshTag][l + 1+d_NEQU*j] / dcelldata[meshTag][0+d_NEQU*j], 2);
        }
        double pressure = (Gamma - 1) * (dcelldata[meshTag][d_dim + 1+d_NEQU*j] - 0.5 * dcelldata[meshTag][0+d_NEQU*j] * VM);
        fout << pressure << '\n';
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

void Euler2D::writeUnacceptableCell(const std::string &filename, int proc, int meshTag, int **dcelldata)
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
    fout << "\"Flag\"\n";
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
    return;
}


void Euler2D::RemoteCellCommunication(double** value)
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

void Euler2D::NearCellCommunication(double** value)
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

void Euler2D::outPutCp(std::string& filename, int mesh_ind)
{
    auto& curBlk = d_hder->blk2D[mesh_ind];
    std::set<int,std::less<int>>* recvCellSet = curBlk.getRecvIndSet();
    std::vector<double> raw_data;
    for(int k = 0;k<d_hder->nEdges(mesh_ind);k++)
    {
        auto& curEdge = curBlk.d_localEdges[k];
        int lC = curEdge.lCInd();
        int rC = curEdge.rCInd();
        if(rC==GeomElements::edge3d<2>::BoundaryType::WALL)
        {
            //Check if lC is inside the recvCell
            if(recvCellSet->find(lC)==recvCellSet->end())
            {
                double cp = (U[mesh_ind][1+d_NEQU*lC] - fsnd_pressure)*2/(fsnd_density*(fsnd_velocity_components[0]*fsnd_velocity_components[0]+fsnd_velocity_components[1]*fsnd_velocity_components[1]));
                raw_data.push_back(curEdge.center()[0]);
                raw_data.push_back(curEdge.center()[1]);
                raw_data.push_back(curEdge.normal_vector()[1]);
                raw_data.push_back(cp);
            }
        }
    }
    int nlocalWall = raw_data.size()/4;
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

    hid_t dataset_id = H5Dcreate(file_id,"Cp",H5T_NATIVE_DOUBLE,dataspace_id,
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

        status = H5Dwrite(dataset_id,H5T_NATIVE_DOUBLE,memspace_id,dataspace_id,H5P_DEFAULT,raw_data.data());

        status = H5Sclose(memspace_id);
    }
    status = H5Sclose(dataspace_id);
    status = H5Dclose(dataset_id);
    status = H5Fclose(file_id);

    MPI_Barrier(MPI_COMM_WORLD);
    if(cur_proc==0)
    {
        file_id = H5Fopen(filename.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
        dataset_id = H5Dopen(file_id,"Cp",H5P_DEFAULT);

        dataspace_id =H5Dget_space(dataset_id);
        double* buffer = new double[nEachLocalWall[num_proc]*4];
        status = H5Dread(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,dataspace_id,H5P_DEFAULT,buffer);
        std::vector<PressureHolder2D> holderVector;
        for(int i = 0;i<nEachLocalWall[num_proc];i++)
        {
            holderVector.push_back(PressureHolder2D(buffer[4*i],buffer[4*i+1],buffer[4*i+2],buffer[4*i+3]));
        }
        std::sort(holderVector.begin(),holderVector.end(),rankWithXY);
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

    delete[] nEachLocalWall;

}