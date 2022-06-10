#pragma once
//There are two types of communication cell,Ones that send and ones that receive.
#include <vector>
class CommunicationCell
{
    public:
    int d_layer_level = -1;
    int d_localId = -1;
    int d_remoteId = -1;
    int d_remoteProc = -1;
    CommunicationCell()
    {

    }
    CommunicationCell(int indHere,int whichProc,int indThere, int level)
    {
        d_localId = indHere;
        d_remoteId = indThere;
        d_remoteProc = whichProc;
        d_layer_level = level;
    }
    void setCommCell(int indHere,int whichProc,int indThere, int level)
    {
        d_localId = indHere;
        d_remoteId = indThere;
        d_remoteProc = whichProc;
        d_layer_level = level;
    }
};

class CommunicationLayerFilter
{
public:
    

    void setnProc(int nProc)
    {
        nCellEachProc.resize(nProc);
        std::fill(nCellEachProc.begin(),nCellEachProc.end(),0);
        filtered_inds.clear();
    }

    void nCellCountPlus(int ind)
    {
        ++nCellEachProc[ind];
    }

    void pushCellInd(int cellInd)
    {
        filtered_inds.push_back(cellInd);
    }

    int layerSize()
    {
        return filtered_inds.size();
    }


    int getCellInd(int f_ind)
    {
        return filtered_inds[f_ind];
    }

    int getCellCountOnProc(int ind)
    {
        return nCellEachProc[ind];
    }
private:
    std::vector<int> filtered_inds;
    std::vector<int> nCellEachProc;
};


bool rankWithSource(const CommunicationCell &e1,const CommunicationCell &e2);
bool rankWithDestin(const CommunicationCell &e1,const CommunicationCell &e2);