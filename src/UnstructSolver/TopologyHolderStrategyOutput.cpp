#include "TopologyHolderStrategy.h"

#include "UnstructIntegrator.h"
#include "mpi.h"
#include <iostream>

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

