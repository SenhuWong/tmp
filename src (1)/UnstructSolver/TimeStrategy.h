#pragma once
// This class should provide an interface
#include "TopologyHolderStrategy.h"
class UnstructTopologyHolder;
class TimeStrategy
{
private:
    TimeStrategy()
    {
    }

protected:
    int d_nmesh = -1;
    int d_NEQU = -1;
    int d_dim = -1;
    UnstructTopologyHolder *d_hder = nullptr;
    TopologyHolderStrategy *d_hder_strategy = nullptr;
    TimeStrategy(UnstructTopologyHolder *hder, TopologyHolderStrategy *hder_strategy);
public:
    virtual void initialize() = 0;
    //Things need to be done before making a time step forward
    //virtual void preprocessUpdate() = 0;
    virtual void singleStep(int curStage) = 0;
    //virtual void advance() = 0;
    //virtual void postprocessUpdate() = 0;
    virtual bool converged() = 0;

};