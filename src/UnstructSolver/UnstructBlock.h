#pragma once
#include "../toolBox/cell3d_int.h"
#include "../toolBox/edge3d_int.h"
#include "commCell.h"
#include <map>
#include <set>
#include <fstream>

class UnstructBlock
{
protected:
    int d_nPs = 0;
    int d_nCs = 0;
    int d_nEs = 0;
    int d_nFluxedEs = 0;//This is for debugging at computation of Flux
    // Pointers to local point,edge and cell should be implemented at subclasses.


    // Communication related Procs for this mesh only is documented at read time.
    // proc id will be reset when communication proc of all meshes is formed.
    int global_nRelatedProcs = 0;
    int *nsendCellEachProc = NULL;
    int *nrecvCellEachProc = NULL;

     //All communication cells.
    std::vector<CommunicationCell> d_to_send;
    std::vector<CommunicationCell> d_to_recv;
public:
    std::map<int, int> *d_proc2Index = NULL; // From proc to local index(after calling reset, From proc to global index)
    std::map<int, int> *d_index2Proc = NULL; // From local index to proc

    //Filters to allow nearCell's sendRecv and remoteCell's sendRecv.
    CommunicationLayerFilter * d_sendRemote = NULL;
    CommunicationLayerFilter * d_sendNear = NULL;
    CommunicationLayerFilter * d_recvRemote = NULL;
    CommunicationLayerFilter * d_recvNear = NULL;

   
    
    
public:
    int nPoints(){return d_nPs;}
    int nCells(){return d_nCs;}
    int nEdges(){return d_nEs;}
    std::vector<CommunicationCell>& sendCommCells()
    {
        return d_to_send;
    }
    std::vector<CommunicationCell>& recvCommCells()
    {
        return d_to_recv;
    }

    UnstructBlock()
    {
        d_sendRemote = new CommunicationLayerFilter[1];
        d_recvRemote = new CommunicationLayerFilter[1];
        d_sendNear = new CommunicationLayerFilter[1];
        d_recvNear = new CommunicationLayerFilter[1];
    }
    ~UnstructBlock()
    {
        //Clear pointers set from outside sources(usually the templated ones).
        clear_external();
        //Clear pointers allocated in constructor.
        if (d_sendRemote){delete[] d_sendRemote;}
        if (d_recvRemote){delete[] d_recvRemote;}
        if (d_sendNear){delete[] d_sendNear;}
        if (d_recvNear){delete[] d_recvNear;}
        if (nsendCellEachProc){delete[] nsendCellEachProc;}
        if (nrecvCellEachProc){delete[] nrecvCellEachProc;}
    }
    virtual void clear_external()
    {

    }
    virtual void computeEdgeCellMetaData() = 0;


    virtual void setLocalStructure(void* localP,
                           void* localC,
                           void* localE,
                           int nlocalP,
                           int nlocalC,
                           int nlocalE,
                           int nFluxedE) = 0;

    void setLocalCommunication(std::map<int, int> *proc2index, std::map<int, int> *index2Proc)
    {
        d_proc2Index = proc2index;
        d_index2Proc = index2Proc;
    }
    void pushSendCell(const CommunicationCell &sendC)
    {
        d_to_send.push_back(sendC);
    }

    void pushRecvCell(const CommunicationCell &recvC)
    {
        d_to_recv.push_back(recvC);
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

    std::vector<CommunicationCell> &getRecvCells()
    {
        return d_to_recv;
    }
    std::vector<CommunicationCell> &getSendCells()
    {
        return d_to_send;
    }
};

template<int ndim>
class UnstructBlock2D: public UnstructBlock
{
public:
    GeomElements::point3d<ndim> *d_localPoints = NULL;
    GeomElements::cell3d<ndim> *d_localCells = NULL;
    GeomElements::edge3d<ndim> *d_localEdges = NULL;

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
            if (rc == GeomElements::edge3d<ndim>::BoundaryType::WALL) // Wall
            {
                flag[lc] += -2;
            }
            else if (rc == GeomElements::edge3d<ndim>::BoundaryType::FARFIELD)
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
    
    void clear_external()//Clear pointers from UnstructFeeder.
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
        if (!d_to_send.empty())
        {
            d_to_send.clear();
        }
        if (!d_to_recv.empty())
        {
            d_to_recv.clear();
        }
        if (d_proc2Index)
        {
            delete[] d_proc2Index;
            d_proc2Index = NULL;
        }
        if (d_index2Proc)
        {
            delete[] d_index2Proc;
            d_index2Proc = NULL;
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
        if (nsendCellEachProc){delete[] nsendCellEachProc;}
        if (nrecvCellEachProc){delete[] nrecvCellEachProc;}
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
        GeomElements::cell3d<ndim>::bindPoints(d_localPoints);
        GeomElements::edge3d<ndim>::bindPoints(d_localPoints);
    }

    void setLocalStructure(void* localP,
                           void* localC,
                           void* localE,
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
        d_localPoints = (GeomElements::point3d<ndim>*) localP;
        d_localCells = (GeomElements::cell3d<ndim>*) localC;
        d_localEdges = (GeomElements::edge3d<ndim>*)localE;
        d_nFluxedEs = nFluxedE;
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