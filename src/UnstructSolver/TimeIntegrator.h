#pragma once
// This class should provide an interface
#include "TopologyHolderStrategy.h"
class UnstructTopologyHolder;
class TimeIntegrator
{
protected:
    int d_nmesh = -1;
    int d_NEQU = -1;
    int d_dim = -1;
    UnstructTopologyHolder *d_hder = nullptr;
    TopologyHolderStrategy *d_hder_strategy = nullptr;
    TimeIntegrator();
    virtual void initializeSubStrategy();
    TimeIntegrator(UnstructTopologyHolder *hder, TopologyHolderStrategy *hder_strategy);
public:
    void initializeStrategy(UnstructTopologyHolder *hder,TopologyHolderStrategy* hder_strategy);
    
    virtual void initializeData() = 0;
    //Things need to be done before making a time step forward
    //virtual void preprocessUpdate() = 0;
    virtual void singleStep(int curStage) = 0;
    //virtual void advance() = 0;
    //virtual void postprocessUpdate() = 0;
    virtual bool converged() = 0;

};