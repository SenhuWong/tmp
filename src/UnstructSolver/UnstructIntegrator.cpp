#include "UnstructIntegrator.h"
#include "TopologyHolderStrategy.h"
//#include "TiogaFeeder.h"
#include "UnstructFeeder.h"
#include <set>
UnstructTopologyHolder::UnstructTopologyHolder(UnstructFeeder *fder)
{
    
    d_nmesh = fder->d_nmesh;
    d_dim = fder->d_dim;

    numproc = fder->num_proc;
    curproc = fder->cur_proc;

    blk2D = fder->UBs2;
    blk3D = fder->UBs3;
    
    if (d_dim == 2)
    {
        std::cout << "MetaData is computed\n";
        for (int i = 0; i < d_nmesh; i++)
        {
            blk2D[i].computeEdgeCellMetaData();
        }
    }
    else if (d_dim == 3)
    {
    }
    else
    {
        throw std::runtime_error("Undefined Dimension\n");
    }
}

void UnstructTopologyHolder::arrange_communication()
{
    // First modify the proc index relationship in UnstructuredBlock
    relatedProcs.clear();
    std::vector<int> nblocks_perProc(numproc, 0);
    for (int i = 0; i < numproc; i++)
    {
        if (d_dim == 2)
        {
            for (int j = 0; j < d_nmesh; j++)
            {
                auto &cur_ub = blk2D[j];

                if (cur_ub.d_proc2Index->find(i) != cur_ub.d_proc2Index->end())
                {
                    relatedProcs.push_back(i);
                    break;
                }
            }
        }
        else if (d_dim == 3)
        {
        }
    }
    if (d_dim == 2)
    {
        for (int i = 0; i < d_nmesh; i++)
        {
            auto &cur_ub = blk2D[i];
            cur_ub.resetProcId(relatedProcs,d_splitter_level);
        }
    }
    // Then calculate how many cells there is needed for each proc
    // like this [Proc0([block0,block1,block2])][Porc1]
    // Therefore we need [nrelatdProc]*[nblock]+1 ptrs to indicate
    if (recvCell_ProcBlockPtr){delete[] recvCell_ProcBlockPtr;}
    if (sendCell_ProcBlockPtr){delete[] sendCell_ProcBlockPtr;}
    recvCell_ProcBlockPtr = new int *[relatedProcs.size()];
    sendCell_ProcBlockPtr = new int *[relatedProcs.size()];
    for (int i = 0; i < relatedProcs.size(); i++)
    {
        recvCell_ProcBlockPtr[i] = new int[d_nmesh + 1];
        sendCell_ProcBlockPtr[i] = new int[d_nmesh + 1];
    }
    // like this [Proc0([block0,block1,block2])][Porc1]
    // Therefore we need [nrelatdProc]*[nblock]+1 ptrs to indicate
    std::vector<int> siOffset(relatedProcs.size(), 0);
    std::vector<int> riOffset(relatedProcs.size(), 0);
    for (int i = 0; i < relatedProcs.size(); i++)
    {
        for (int j = 0; j < d_nmesh; j++)
        {
            recvCell_ProcBlockPtr[i][j] = riOffset[i];
            sendCell_ProcBlockPtr[i][j] = siOffset[i];
            if (d_dim == 2)
            {
                auto &cur_ub = blk2D[j];
                riOffset[i] += cur_ub.getNRecvCellOnProcInd(i);
                siOffset[i] += cur_ub.getNSendCellOnProcInd(i);
            }
        }
        recvCell_ProcBlockPtr[i][d_nmesh] = riOffset[i];
        sendCell_ProcBlockPtr[i][d_nmesh] = siOffset[i];
    }

    if (remoteRecvCell_ProcBlockPtr){delete[] remoteRecvCell_ProcBlockPtr;}
    if (remoteSendCell_ProcBlockPtr){delete[] remoteSendCell_ProcBlockPtr;}
    if (nearRecvCell_ProcBlockPtr){delete[] nearRecvCell_ProcBlockPtr;}
    if (nearSendCell_ProcBlockPtr){delete[] nearSendCell_ProcBlockPtr;}
    remoteRecvCell_ProcBlockPtr = new int*[relatedProcs.size()];
    remoteSendCell_ProcBlockPtr = new int*[relatedProcs.size()];
    nearRecvCell_ProcBlockPtr = new int*[relatedProcs.size()];
    nearSendCell_ProcBlockPtr = new int*[relatedProcs.size()];
    for(int i = 0;i<relatedProcs.size();i++)
    {
        remoteRecvCell_ProcBlockPtr[i] = new int[d_nmesh+1];
        remoteSendCell_ProcBlockPtr[i] = new int[d_nmesh+1];
        nearRecvCell_ProcBlockPtr[i] = new int[d_nmesh+1];
        nearSendCell_ProcBlockPtr[i] = new int[d_nmesh+1];
    }

    std::fill(siOffset.begin(),siOffset.end(),0);
    std::fill(riOffset.begin(),riOffset.end(),0);

    for(int i = 0;i<relatedProcs.size();i++)
    {
        for(int j = 0;j<d_nmesh;j++)
        {
            remoteRecvCell_ProcBlockPtr[i][j] = riOffset[i];
            remoteSendCell_ProcBlockPtr[i][j] = siOffset[i];
            if(d_dim==2)
            {
                auto& cur_ub = blk2D[j];
                riOffset[i] += cur_ub.getNRemoteRecvCellOnProcInd(i);
                siOffset[i] += cur_ub.getNRemoteSendCellOnProcInd(i);
            }
        }
        remoteRecvCell_ProcBlockPtr[i][d_nmesh] = riOffset[i];
        remoteSendCell_ProcBlockPtr[i][d_nmesh] = siOffset[i];
    }

    std::fill(siOffset.begin(),siOffset.end(),0);
    std::fill(riOffset.begin(),riOffset.end(),0);

    for(int i = 0;i<relatedProcs.size();i++)
    {
        for(int j = 0;j<d_nmesh;j++)
        {
            nearRecvCell_ProcBlockPtr[i][j] = riOffset[i];
            nearSendCell_ProcBlockPtr[i][j] = siOffset[i];
            if(d_dim==2)
            {
                auto& cur_ub = blk2D[j];
                riOffset[i] += cur_ub.getNNearRecvCellOnProcInd(i);
                siOffset[i] += cur_ub.getNNearSendCellOnProcInd(i);
            }
        }
        nearRecvCell_ProcBlockPtr[i][d_nmesh] = riOffset[i];
        nearSendCell_ProcBlockPtr[i][d_nmesh] = siOffset[i];
    }
    
    // if (curproc != 0)
    //     return;
    // for (int i = 0; i < relatedProcs.size(); i++)
    // {
    //     std::cout << "Index of related Proc is " << i << " ";
    //     for (int j = 0; j < d_nmesh + 1; j++)
    //     {
    //         std::cout << recvCell_ProcBlockPtr[i][j] << '\t';
    //     }
    //     std::cout << '\n';
    // }
    // for (int i = 0; i < relatedProcs.size(); i++)
    // {
    //     std::cout << "Index of related Proc is " << i << " ";
    //     for (int j = 0; j < d_nmesh + 1; j++)
    //     {
    //         std::cout << sendCell_ProcBlockPtr[i][j] << '\t';
    //     }
    //     std::cout << '\n';
    // }
}

void UnstructTopologyHolder::writeGridMetaData(const std::string &filename, int meshTag)
{
    std::ofstream fout;

    std::string localEdgeFilename = filename + "_edge_" + std::to_string(meshTag) + std::to_string(curproc);
    fout.open(localEdgeFilename);
    GeomElements::edge3d<2>::bindPoints(blk2D[meshTag].d_localPoints);
    for (int i = 0; i < blk2D[meshTag].d_nEs; i++)
    {
        fout << "Edge " << i << '\n';
        for (int j = 0; j < blk2D[meshTag].d_localEdges[i].size(); j++)
        {
            auto &curNode = blk2D[meshTag].d_localPoints[blk2D[meshTag].d_localEdges[i].pointInd(j)];
            //std::cout<<blk2D[meshTag].d_localEdges->pointInd(j)<<'\n';

            fout << curNode[0] << " " << curNode[1] << '\n';
        }
        auto curNorm = blk2D[meshTag].d_localEdges[i].normal_vector();
        auto curCenter = blk2D[meshTag].d_localEdges[i].center();
        fout << "area: " << blk2D[meshTag].d_localEdges[i].area() << '\n';
        fout << "center: " << curCenter[0] << " , " << curCenter[1] << '\n';
        fout << "normal: " << curNorm[0] << " , " << curNorm[1] << '\n';
    }
    fout.close();

    std::string localCellFilename = filename + "_cell_" + std::to_string(meshTag) + std::to_string(curproc);

    fout.open(localCellFilename);
    GeomElements::cell3d<2>::bindPoints(blk2D[meshTag].d_localPoints);
    for (int i = 0; i < blk2D[meshTag].d_nCs; i++)
    {
        fout << "Cell " << i << '\n';
        for (int j = 0; j < blk2D[meshTag].d_localCells[i].size(); j++)
        {
            auto &curNode = blk2D[meshTag].d_localPoints[blk2D[meshTag].d_localCells[i].pointInd(j)];
            fout << curNode[0] << " " << curNode[1] << '\n';
        }
        auto curCenter = blk2D[meshTag].d_localCells[i].center();
        fout << "volume: " << blk2D[meshTag].d_localCells[i].volume() << '\n';
        fout << "center: " << curCenter[0] << " , " << curCenter[1] << '\n';
    }
    fout.close();
}