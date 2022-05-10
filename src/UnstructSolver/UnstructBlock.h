#pragma once
#include "../toolBox/cell3d_int.h"
#include "../toolBox/edge3d_int.h"
#include "commCell.h"
#include <map>
#include <set>
#include <fstream>

class UnstructBlock2D
{
public:
    // Geometric Elements;
    int d_nPs = 0;
    int d_nCs = 0;
    int d_nEs = 0;

    int d_nFluxedEs = 0;


    GeomElements::point3d<2> *d_localPoints = NULL;
    GeomElements::cell3d<2> *d_localCells = NULL;
    GeomElements::edge3d<2> *d_localEdges = NULL;
    // Communication Cells;
    // CommunicationProcs on this mesh only is documented at read time.
    // proc id will be modified when communication proc of all meshes is considered.
    int global_nRelatedProcs = 0;
    int *nsendCellEachProc = NULL;
    int *nrecvCellEachProc = NULL;

    // The new added near and remote cell info.
    // int *nnearSendCellEachProc = NULL;
    // int *nnearRecvCellEachProc = NULL;
    // int nnearSend = 0;
    // int *nearSendCellInds = NULL;
    // int nnearRecv = 0;
    // int *nearRecvCellInds = NULL;
    
    CommunicationLayerFilter * d_sendRemote;
    CommunicationLayerFilter * d_sendNear;
    CommunicationLayerFilter * d_recvRemote;
    CommunicationLayerFilter * d_recvNear;

    // int *nremoteSendCellEachProc = NULL;
    // int *nremoteRecvCellEachProc = NULL;
    // int nremoteSend = 0;
    // int *remoteSendCellInds = NULL;
    // int nremoteRecv = 0;
    // int *remoteRecvCellInds = NULL;

    std::vector<CommunicationCell> d_to_send;
    std::vector<CommunicationCell> d_to_recv;
    std::map<int, int> *d_proc2Index = NULL; // From proc to local index(after calling reset, From proc to global index)
    std::map<int, int> *d_index2Proc = NULL; // From local index to proc
    //std::vector<int> localInd2AggregatedInd;
    void write_grd(const std::string &filename, int meshtag, int cur_proc)
    {
        std::string localFilename = filename + std::to_string(meshtag) + std::to_string(cur_proc) + ".dat";
        std::ofstream fout;

        fout.open(localFilename);
        if (!fout.is_open())
        {
            std::cout << "Open file failure\n";
            return;
        }
        // Assume that the vconn is organized in a Tecplot friendly style
        fout << "TITLE =\"Tioga output\"\n";
        fout << "VARIABLES=\"X\",\"Y\",";
        fout << "ZONE T=\"VOL_MIXED\",N=" << d_nPs << " E=" << d_nCs;

        fout << " ET = QUADRILATERAL, F = FEPOINT\n";

        for (int i = 0; i < d_nPs; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                fout << d_localPoints[i][j] << " ";
            }
            fout << '\n';
        }
        for (int i = 0; i < d_nCs; i++)
        {
            for (int j = 0; j < d_localCells[i].size(); j++)
            {
                fout << d_localCells[i].pointInd(j) + 1 << " ";
            }
            for (int j = 0; j < 4 - d_localCells[i].size(); j++)
            {
                fout << d_localCells[i].pointInd(d_localCells[i].size() - 1) + 1 << " ";
            }
            fout << '\n';
        }
        fout.close();
    }

    void write_grdSendRecv(const std::string &filename, int meshtag, int cur_proc)
    {
        std::string localFilename = filename + std::to_string(meshtag) + std::to_string(cur_proc) + ".dat";
        std::ofstream fout;

        fout.open(localFilename);
        if (!fout.is_open())
        {
            std::cout << "Open file failure\n";
            return;
        }
        // Assume that the vconn is organized in a Tecplot friendly style
        fout << "TITLE =\"Tioga output\"\n";
        fout << "VARIABLES=\"X\",\"Y\",\"SendRecv\",";
        fout << "ZONE T=\"VOL_MIXED\",N=" << d_nPs << " E=" << d_nCs;

        fout << " ET = QUADRILATERAL, F = FEBLOCK\n";
        fout << "VARLOCATION =  (1=NODAL,2=NODAL,3=CELLCENTERED)\n";
        
        for (int j = 0; j < 2; j++)
        {
            for (int i = 0; i < d_nPs; i++)
            {
                fout << d_localPoints[i][j] << "\n";
            }
        }
        int* sendRecvFlag = new int[d_nCs];
        for(int i = 0;i<d_nCs;i++)
        {
            sendRecvFlag[i] = -1;
        }
        for(auto iter = d_to_send.begin();iter!=d_to_send.end();iter++)
        {
            sendRecvFlag[iter->d_localId] = 1;
        }
        for(auto iter = d_to_recv.begin();iter!=d_to_recv.end();iter++)
        {
            sendRecvFlag[iter->d_localId] = 2;
        }
        for(int i = 0;i<d_nCs;i++)
        {
            fout << sendRecvFlag[i]<<'\n';
        }
        for (int i = 0; i < d_nCs; i++)
        {
            for (int j = 0; j < d_localCells[i].size(); j++)
            {
                fout << d_localCells[i].pointInd(j) + 1 << " ";
            }
            for (int j = 0; j < 4 - d_localCells[i].size(); j++)
            {
                fout << d_localCells[i].pointInd(d_localCells[i].size() - 1) + 1 << " ";
            }
            fout << '\n';
        }
        fout.close();
    }

    void writeEdgeLeftRight(const std::string &filename, int meshtag, int cur_proc)
    {
        std::cout << "WriteEdgeLeftRight called\n";
        // std::cin.get();
        std::ofstream fout;
        std::string local_filename = filename + std::to_string(meshtag) + std::to_string(cur_proc);
        fout.open(local_filename);
        for (int i = 0; i < d_nEs; i++)
        {
            fout << "i" << d_localEdges[i].lCInd() << " " << d_localEdges[i].rCInd() << '\n';
            // fout<<W[meshtag][0][d_localEdges[i].lCInd()]<<" "<<W[meshtag][0][d_localEdges[i].rCInd()];
        }
        fout.close();
    }

    void writeEdge(const std::string &filename, int meshtag, int cur_proc)
    {
        std::ofstream fout;
        std::string local_filename = filename + std::to_string(meshtag) + std::to_string(cur_proc);
        fout.open(local_filename);
        int *flag = new int[d_nCs];
        for (int i = 0; i < d_nCs; i++)
        {
            flag[i] = -1;
        }
        for (int i = 0; i < d_nEs; i++)
        {
            auto &curEdge = d_localEdges[i];
            int lc = d_localEdges[i].lCInd();
            int rc = d_localEdges[i].rCInd();
            if (rc == GeomElements::edge3d<2>::BoundaryType::WALL) // Wall
            {
                flag[lc] += -2;
            }
            else if (rc == GeomElements::edge3d<2>::BoundaryType::FARFIELD)
            {
                flag[lc] += -3;
            }
            else
            {
                flag[lc] += 1;
                flag[rc] += 1;
            }
        }
        fout << "TITLE =\"Tioga output\"\n";
        fout << "VARIABLES=\"X\",\"Y\",\"FLAG\"";
        fout << "ZONE T=\"VOL_MIXED\",N=" << d_nPs << " E=" << d_nCs;

        fout << " ET = QUADRILATERAL, F = FEBLOCK\n";
        fout << "VARLOCATION =  (1=NODAL,2=NODAL,3=CELLCENTERED)\n";
        for (int j = 0; j < 2; j++)
        {
            for (int i = 0; i < d_nPs; i++)
            {

                fout << d_localPoints[i][j] << "\n";
            }
        }
        for (int i = 0; i < d_nCs; i++)
        {
            fout << flag[i] << '\n';
        }
        for (int i = 0; i < d_nCs; i++)
        {
            for (int j = 0; j < d_localCells[i].size(); j++)
            {
                fout << d_localCells[i].pointInd(j) + 1 << " ";
            }
            for (int j = 0; j < 4 - d_localCells[i].size(); j++)
            {
                fout << d_localCells[i].pointInd(d_localCells[i].size() - 1) + 1 << " ";
            }
            fout << '\n';
        }
        fout.close();
        delete[] flag;
    }
    
    void clear_external()
    {
        if (d_localPoints)
        {
            delete[] d_localPoints;
            d_localPoints = NULL;
        }
        if (d_localCells)
        {
            delete[] d_localCells;
            d_localCells = NULL;
        }
        if (d_localEdges)
        {
            delete[] d_localEdges;
            d_localEdges = NULL;
        }
        if (nsendCellEachProc)
        {
            delete[] nsendCellEachProc;
        }
        if (nrecvCellEachProc)
        {
            delete[] nrecvCellEachProc;
        }
        if (!d_to_send.empty())
        {
            d_to_send.clear();
        }
        if (!d_to_recv.empty())
        {
            d_to_recv.clear();
        }
        
        
    }
    UnstructBlock2D()
    {
        d_sendRemote = new CommunicationLayerFilter[1];
        d_recvRemote = new CommunicationLayerFilter[1];
        d_sendNear = new CommunicationLayerFilter[1];
        d_recvNear = new CommunicationLayerFilter[1];
    }
    ~UnstructBlock2D()
    {
        //Clear pointers set from outside sources.
        clear_external();
        //Clear pointers allocated in constructor.
        if (d_sendRemote){delete[] d_sendRemote;}
        if (d_recvRemote){delete[] d_recvRemote;}
        if (d_sendNear){delete[] d_sendNear;}
        if (d_recvNear){delete[] d_recvNear;}
    }
    
    void computeEdgeCellMetaData()
    {
        bindCurPoints();
        
        for (int i = 0; i < d_nCs; i++)
        {
            d_localCells[i].computeMetaData();
        }
        for (int i = 0; i < d_nEs; i++)
        {
            d_localEdges[i].computeMetaData();
        }
    }

    void bindCurPoints()
    {
        GeomElements::cell3d<2>::bindPoints(d_localPoints);
        GeomElements::edge3d<2>::bindPoints(d_localPoints);
    }

    void setLocalStructure(GeomElements::point3d<2> *localP,
                           GeomElements::cell3d<2> *localC,
                           GeomElements::edge3d<2> *localE,
                           int nlocalP,
                           int nlocalC,
                           int nlocalE,
                           int nFluxedE)
    {
        if (d_localPoints)
        {
            delete[] d_localPoints;
            d_localPoints = NULL;
        }
        if (d_localCells)
        {
            delete[] d_localCells;
            d_localCells = NULL;
        }
        if (d_localEdges)
        {
            delete[] d_localEdges;
            d_localEdges = NULL;
        }
        d_nPs = nlocalP;
        d_nCs = nlocalC;
        d_nEs = nlocalE;
        d_localPoints = localP;
        d_localCells = localC;
        d_localEdges = localE;
        d_nFluxedEs = nFluxedE;
    }

    void setLocalCommunication(std::map<int, int> *proc2index, std::map<int, int> *index2Proc)
    {
        d_proc2Index = proc2index;
        d_index2Proc = index2Proc;
    }
    void pushSendCell(const CommunicationCell &sendC)
    {
        // std::cout<<"Send cell pushed\n";
        d_to_send.push_back(sendC);
    }

    void pushRecvCell(const CommunicationCell &recvC)
    {
        std::cout << "Recv cell pushed\n";
        d_to_recv.push_back(recvC);
        std::cout << d_to_recv.size() << '\n';
    }
    // According to the given global view of realtedProcs, change the proc2Ind and count how many cells there is on each ind.
    void resetProcId(const std::vector<int> &aggregated_relatedProcs, int level_splitter)
    {
        // The proc to ind is changed to be the globalized ind.
        global_nRelatedProcs = aggregated_relatedProcs.size();
        // Make the map of global proc id to index.
        int ind = 0;
        for (auto iter = aggregated_relatedProcs.begin(); iter < aggregated_relatedProcs.end(); iter++)
        {
            auto to_search = d_proc2Index->find(*iter);
            if (to_search != d_proc2Index->end())
            {
                (*d_proc2Index)[to_search->first] = to_search->second + ind;
            }
            else
            {
                ind++;
            }
        }
        // Apply that map to each CommCell.
        for (auto iter = d_to_send.begin(); iter != d_to_send.end(); iter++)
        {
            int local_ind = iter->d_remoteProc;
            int global_ind = (*d_proc2Index)[(*d_index2Proc)[local_ind]];
            iter->d_remoteProc = global_ind;
        }
        for (auto iter = d_to_recv.begin(); iter != d_to_recv.end(); iter++)
        {
            int local_ind = iter->d_remoteProc;
            int global_ind = (*d_proc2Index)[(*d_index2Proc)[local_ind]];
            iter->d_remoteProc = global_ind;
        }

        if (nsendCellEachProc){delete[] nsendCellEachProc;}
        if (nrecvCellEachProc){delete[] nrecvCellEachProc;}
        
        nsendCellEachProc = new int[global_nRelatedProcs];
        nrecvCellEachProc = new int[global_nRelatedProcs];
        for (int i = 0; i < global_nRelatedProcs; i++)
        {
            nsendCellEachProc[i] = 0;
            nrecvCellEachProc[i] = 0;
        }
        
        d_sendRemote->setnProc(global_nRelatedProcs);
        d_recvRemote->setnProc(global_nRelatedProcs);
        d_sendNear->setnProc(global_nRelatedProcs);
        d_recvNear->setnProc(global_nRelatedProcs);
        int splitting_level = level_splitter;
        // nnearRecv = 0;
        // nremoteRecv = 0;
        for (auto iter = d_to_recv.begin(); iter != d_to_recv.end(); iter++)
        {
            if(iter->d_layer_level>splitting_level)
            {
                d_recvRemote->pushCellInd(iter-d_to_recv.begin());
                d_recvRemote->nCellCountPlus(iter->d_remoteProc);
            }
            else
            {
                d_recvNear->pushCellInd(iter-d_to_recv.begin());
                d_recvNear->nCellCountPlus(iter->d_remoteProc);
            }
            ++nrecvCellEachProc[iter->d_remoteProc];
        }

        for (auto iter = d_to_send.begin(); iter != d_to_send.end(); iter++)
        {
            if(iter->d_layer_level>splitting_level)
            {
                d_sendRemote->pushCellInd(iter-d_to_send.begin());
                d_sendRemote->nCellCountPlus(iter->d_remoteProc);
                
            }
            else
            {
                d_sendNear->pushCellInd(iter - d_to_send.begin());
                d_sendNear->nCellCountPlus(iter->d_remoteProc);
            }
            ++nsendCellEachProc[iter->d_remoteProc];
        }
    }
    
    int getNRemoteSendCell()
    {
        return d_sendRemote->layerSize();
    }

    int getNNearSendCell()
    {
        return d_sendNear->layerSize();
    }

    int getNRemoteRecvCell()
    {
        return d_recvRemote->layerSize();
    }

    int getNNearRecvCell()
    {
        return d_recvNear->layerSize();
    }

    int getNSendCellOnProcInd(int global_id)
    {
        return nsendCellEachProc[global_id];
    }
    int getNRecvCellOnProcInd(int global_id)
    {
        return nrecvCellEachProc[global_id];
    }
    int getNRemoteSendCellOnProcInd(int global_id)
    {
        return d_sendRemote->getCellCountOnProc(global_id);
    }
    int getNRemoteRecvCellOnProcInd(int global_id)
    {
        return d_recvRemote->getCellCountOnProc(global_id);
    }
    int getNNearSendCellOnProcInd(int global_id)
    {
        return d_sendNear->getCellCountOnProc(global_id);
    }
    int getNNearRecvCellOnProcInd(int global_id)
    {
        return d_recvNear->getCellCountOnProc(global_id);
    }

    //This is used only for output and therefore don't care about performance

    std::set<int,std::less<int>>* getRecvIndSet()
    {
        std::set<int,std::less<int>>* recvInds =  new std::set<int,std::less<int>>[1];
        for(auto iter = d_to_recv.begin();iter!=d_to_recv.end();iter++)
        {
            recvInds->insert(iter->d_localId);
        }
        return recvInds;
    } 

    

    // std::vector<CommunicationCell>::iterator remoteRecvBegin(int cellid)
    // {
    //     return d_to_recv.begin() + cellid;
    // }
    // std::vector<CommunicationCell>::iterator remoteSendBegin(int cellid)
    // {


    // }
    // std::vector<CommunicationCell>::iterator nearRecvBegin(int cellid)
    // {

    // }
    // std::vector<CommunicationCell>::iterator nearSendBegin(int cellid)
    // {

    // }


    std::vector<CommunicationCell> &getRecvCells()
    {
        return d_to_recv;
    }
    std::vector<CommunicationCell> &getSendCells()
    {
        return d_to_send;
    }
};

class UnstructBlock3D
{
public:
    std::vector<GeomElements::point3d<3>> d_localPoints;
    std::vector<GeomElements::cell3d<3>> d_localCells;
    std::vector<GeomElements::edge3d<3>> d_edges;
    CommunicationCell *d_to_send = NULL;
    CommunicationCell *d_to_recv = NULL;
};