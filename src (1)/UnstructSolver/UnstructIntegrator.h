//This class holds all the local MeshBlock and maintains the edge2cell connection provided by BlockInfo.
//It manages uses a flux computation strategy to alocate, solve and communicate the flux on the edge.
//The cell index is different from that in TIOGA.Whenever a cell want to locate a node it needs to be transfered first into TIOGA style so
//the vconn could be accessed correctly.An edge 's left right cell is obtained directly in this, but still needs a map to TIOGA style if it is used in tioga.
//
#pragma once
#include <vector>
#include "UnstructBlock.h"

//Todo:: To enable more efficient communication at LU-SGS, we might want to add a filter for each level of communication cells.
//It is impossible to list all combination, so the more practicle way is to take in a splitting level_num, then 3 groups of communication cell is formed
//One that is below(or all), one that is above and all.
class UnstructFeeder;
class UnstructTopologyHolder
{
    friend class TopologyHolderStrategy;
public:
    int d_nmesh = 0;
    int d_dim = -1;
    //std::vector<int> mtags;
    //std::map<int,int> tag_iblk_map;//Tag to index of the block
    int numproc = -1;
    int curproc = -1;
    
    //New added feature
    //
    int d_splitter_level = 1;
    int* recvCellInd_level[2] = {NULL};
    int* sendCellInd_level[2] = {NULL};

    std::vector<int> relatedProcs;
public:
    //Number of mesh with number of cells,edges and node on each mesh.
    
    UnstructBlock2D<2>* blk2D = NULL;
    UnstructBlock2D<3>* blk3D = NULL;
    //Find the aggregated related CommProcs;
    int** recvCell_ProcBlockPtr = NULL;
    int** sendCell_ProcBlockPtr = NULL;

    int** nearRecvCell_ProcBlockPtr = NULL;
    int** nearSendCell_ProcBlockPtr = NULL;
    
    int** remoteRecvCell_ProcBlockPtr = NULL;
    int** remoteSendCell_ProcBlockPtr = NULL;

public:
    UnstructTopologyHolder(UnstructFeeder *fder); //Maybe later we might need a reading constructor
    void arrange_communication();
    void initialize_integrator();
    int nPoints(int blkid)
    {
        if(d_dim==2)
        {
            return blk2D[blkid].nPoints();
        }
        else if(d_dim ==3)
        {
            return -1;
        }
        else
        {
            return -1;
        }
    }
    int nCells(int blkid)
    {
        if(d_dim==2)
        {
            return blk2D[blkid].nCells();
        }
        else if(d_dim==3)
        {
            return -1;
        }
        else
        {
            return -1;
        }
        
    }
    int nEdges(int blkid)
    {
        if(d_dim==2)
        {
            return blk2D[blkid].nEdges();
        }
        else if(d_dim==3)
        {
            return -1;
        }
        else
        {
            return -1;
        }
    }
    int sendRemoteCellPtr(int procInd,int blockInd)
    {
        return remoteSendCell_ProcBlockPtr[procInd][blockInd]; 
    }
    int recvRemoteCellPtr(int procInd,int blockInd)
    {
        return remoteRecvCell_ProcBlockPtr[procInd][blockInd];
    }

    int sendNearCellPtr(int procInd,int blockInd)
    {
        return nearSendCell_ProcBlockPtr[procInd][blockInd];
    }

    int recvNearCellPtr(int procInd,int blockInd)
    {
        return nearRecvCell_ProcBlockPtr[procInd][blockInd];
    }

    int sendCellPtr(int procInd,int blockInd)
    {
        return sendCell_ProcBlockPtr[procInd][blockInd];
    }
    int recvCellPtr(int procInd,int blockInd)
    {
        return recvCell_ProcBlockPtr[procInd][blockInd];
    }
    double EdgeNormalComponent(int blockInd,int edgeInd,int dirInd)
    {
        if(d_dim==2)
        {
            return blk2D[blockInd].d_localEdges[edgeInd].normal_vector()[dirInd];
        }
        else if(d_dim==3)
        {
            return -1;
        }
    }
    double EdgeArea(int blockInd,int edgeInd)
    {
        if(d_dim==2)
        {
            return blk2D[blockInd].d_localEdges[edgeInd].area();
        }
        else if(d_dim==3)
        {
            return -1;
        }
    }
    double CellVolume(int blockInd,int cellInd)
    {
        if(d_dim==2)
        {
            return blk2D[blockInd].d_localCells[cellInd].volume();
        }
        else if(d_dim==3)
        {
            return -1;
        }
    }
    std::vector<CommunicationCell>::iterator sendBegin(int blkid)
    {
        if(d_dim==2)
        {
            return blk2D[blkid].sendCommCells().begin();
        }
        else
        {
            return blk2D[blkid].sendCommCells().begin();
        }
    }
    std::vector<CommunicationCell>::iterator sendEnd(int blkid)
    {
        if(d_dim==2)
        {
            return blk2D[blkid].sendCommCells().end();
        }
        else
        {
            return blk2D[blkid].sendCommCells().end();
        }
    }
    std::vector<CommunicationCell>::iterator recvBegin(int blkid)
    {
        if(d_dim==2)
        {
            return blk2D[blkid].recvCommCells().begin();
        }
        else
        {
            return blk2D[blkid].recvCommCells().begin();
        }
    }
    std::vector<CommunicationCell>::iterator recvEnd(int blkid)
    {
        if(d_dim==2)
        {
            return blk2D[blkid].recvCommCells().end();
        }
        else
        {
            return blk2D[blkid].recvCommCells().end();
        }
    }
    
    int nRemoteSend(int blkid)
    {
        if(d_dim==2)
        {
            return blk2D[blkid].getNRemoteSendCell();
        }
    }

    int nNearSend(int blkid)
    {
        if(d_dim==2)
        {
            return blk2D[blkid].getNNearSendCell();
        }
    }

    int nRemoteRecv(int blkid)
    {
        if(d_dim==2)
        {
            return blk2D[blkid].getNRemoteRecvCell();
        }
    }

    int nNearRecv(int blkid)
    {
        if(d_dim==2)
        {
            return blk2D[blkid].getNNearRecvCell();
        }
    }

    int remoteSendInd(int blkid,int cellId_filter)
    {
        if(d_dim==2)
        {
            return blk2D[blkid].d_sendRemote->getCellInd(cellId_filter);
        }
    }

    int remoteRecvInd(int blkid,int cellId_filter)
    {
        if(d_dim==2)
        {
            return blk2D[blkid].d_recvRemote->getCellInd(cellId_filter);
        }
    }

    int nearSendInd(int blkid,int cellId_filter)
    {
        if(d_dim==2)
        {
            return blk2D[blkid].d_sendNear->getCellInd(cellId_filter);
        }
    }

    int nearRecvInd(int blkid,int cellId_filter)
    {
        if(d_dim==2)
        {
            return blk2D[blkid].d_recvNear->getCellInd(cellId_filter);
        }
    }


    std::vector<CommunicationCell>::iterator remoteSend(int blkid,int cellId_filter)
    {
        if(d_dim==2)
        {
            return blk2D[blkid].sendCommCells().begin() + remoteSendInd(blkid,cellId_filter);
        }
    }

    std::vector<CommunicationCell>::iterator remoteRecv(int blkid,int cellId_filter)
    {
        if(d_dim==2)
        {
            return blk2D[blkid].recvCommCells().begin() + remoteRecvInd(blkid,cellId_filter);
        }
    }

    std::vector<CommunicationCell>::iterator nearSend(int blkid,int cellId_filter)
    {
        if(d_dim==2)
        {
            return blk2D[blkid].sendCommCells().begin() + nearSendInd(blkid,cellId_filter);
        }
    }

    std::vector<CommunicationCell>::iterator nearRecv(int blkid,int cellId_filter)
    {
        if(d_dim==2)
        {
            return blk2D[blkid].recvCommCells().begin() + nearRecvInd(blkid,cellId_filter);
        }

    }

    void writeGridMetaData(const std::string& filename,int mesh);
    void write(const std::string& filename)
    {
        for(int i =0;i<d_nmesh;i++)
        {
            writeGridMetaData(filename,i);
        }
    }
};