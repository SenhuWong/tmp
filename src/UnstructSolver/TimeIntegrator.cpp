#include "TimeIntegrator.h"
#include "TopologyHolderStrategy.h"
#include "UnstructIntegrator.h"
TimeIntegrator::TimeIntegrator()
{

}

void TimeIntegrator::initializeSubStrategy()
{

}

void TimeIntegrator::initializeStrategy(UnstructTopologyHolder *hder,TopologyHolderStrategy* hder_strategy)
{
    d_hder = hder;
    d_hder_strategy = hder_strategy;

    d_nmesh = d_hder_strategy->getNMesh();
    d_NEQU = d_hder_strategy->getNEquation();
    d_dim = d_hder_strategy->getNDim();
    initializeSubStrategy();
}

   


TimeIntegrator::TimeIntegrator(UnstructTopologyHolder *hder, TopologyHolderStrategy *hder_strategy)
    : d_nmesh(hder_strategy->getNMesh()),
    d_NEQU(hder_strategy->getNEquation()), 
    d_dim(hder_strategy->getNDim()),
    d_hder(hder), d_hder_strategy(hder_strategy)
{
     
    std::cout<<"Entering time strategy\n";

}