
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

public:
    virtual void computeFlux() = 0;
    void setComm(int num,int cur)
    {
        num_proc = num;
        cur_proc = cur;

    }

};
