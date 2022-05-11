#include "TimeStrategy.h"
#include "TopologyHolderStrategy.h"
#include "UnstructIntegrator.h"
TimeStrategy::TimeStrategy(UnstructTopologyHolder *hder, TopologyHolderStrategy *hder_strategy)
    : d_nmesh(hder_strategy->getNMesh()),
    d_NEQU(hder_strategy->getNEquation()), 
    d_dim(hder_strategy->getNDim()),
    d_hder(hder), d_hder_strategy(hder_strategy)
{
     
    std::cout<<"Entering time strategy\n";

}