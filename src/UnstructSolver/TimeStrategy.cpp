#include "TimeStrategy.h"
#include "Euler2D.h"
#include "UnstructIntegrator.h"
TimeStrategy::TimeStrategy(UnstructTopologyHolder *hder, Euler2D *hder_strategy)
    : d_nmesh(hder->d_nmesh), d_NEQU(hder_strategy->d_NEQU), d_dim(hder->d_dim), d_hder(hder), d_hder_strategy(hder_strategy)
{
    std::cout<<"Entering time strategy\n";

}