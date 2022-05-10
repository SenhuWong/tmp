
#pragma once
// This class should provide an interface for FLux computation for Euler
#include<stdio.h>
class Euler2D;
class UnstructTopologyHolder;
class FluxStrategy
{
protected:
    int num_proc;
    int cur_proc;
    int d_nmesh = -1;
    int d_NEQU = -1;
    int d_dim = -1;
    UnstructTopologyHolder *d_hder = NULL;
    Euler2D *d_hder_strategy = NULL;
    FluxStrategy()
    {
    }
    FluxStrategy(UnstructTopologyHolder *hder, Euler2D *hder_strategy);
    //     : d_nmesh(hder->d_nmesh), d_NEQU(hder_strategy->d_NEQU), d_dim(hder->d_dim), d_hder(hder),d_hder_strategy(hder_strategy)
    // {
        
    // }
    

public:
    virtual void computeFlux() = 0;
    virtual void computeFlux1() = 0;
    void setComm(int num,int cur)
    {
        num_proc = num;
        cur_proc = cur;

    }

};
