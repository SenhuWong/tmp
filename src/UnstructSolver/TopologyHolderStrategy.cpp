#include "TopologyHolderStrategy.h"
#include <iostream>
TopologyHolderStrategy::TopologyHolderStrategy(UnstructTopologyHolder *hder)
    : d_hder(hder)
{
    //std::cout<<d_hder<<" TopologyHolderStrategy Constructor called\n";
    
}

TopologyHolderStrategy::~TopologyHolderStrategy()
{
}