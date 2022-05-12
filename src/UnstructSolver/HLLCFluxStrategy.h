#pragma once
#include "FluxStrategy.h"

// Given 4 conservatives, calculate the pressure and sound speed
struct FlowParamHolder
{
    double density;
    double pressure;
    double soundSpeed;
    double velocity[3];
};
void computeFlowParams(double ***W, int mesh_ind, int cell_ind, int dim, FlowParamHolder *hder);
class HLLCFluxStrategy : public FluxStrategy
{
public:
    HLLCFluxStrategy(UnstructTopologyHolder *hder, TopologyHolderStrategy *hder_strategy)
        : FluxStrategy(hder, hder_strategy)
    {
    }
    void computeFlux();
};