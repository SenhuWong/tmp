#include "FluxStrategy.h"
#include "Euler2D.h"
FluxStrategy::FluxStrategy(UnstructTopologyHolder *hder, Euler2D *hder_strategy)
    : d_nmesh(hder->d_nmesh), d_NEQU(hder_strategy->d_NEQU), d_dim(hder->d_dim), d_hder(hder), d_hder_strategy(hder_strategy)
{
    
}