#include "LimiterStrategy.h"
#include "Euler2D.h"

LimiterStrategy::LimiterStrategy(UnstructTopologyHolder *hder, Euler2D *hder_strategy)
    : d_nmesh(hder->d_nmesh), d_hder(hder), d_hder_strategy(hder_strategy)
{
    d_dim = d_hder->d_dim;
    d_NEQU = d_hder_strategy->d_NEQU;
}
LimiterStrategy::~LimiterStrategy()
{

}