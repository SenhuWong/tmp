#pragma once
#include<stdio.h>
//This should be an abstract base class for Euler Equation or N-S Equation advancing.
//We need to specify the following routines to perform data management on TopologyHolder:
//1.Allocate data storage for Conservative Variables' Data
//2.Allocate data storage for Flux Variables' Data
//3.Estimate and allocate buffer storage for cross-process communication.
//4.Register subStrategies for Flux computation(by passing left right values or something else) and so on.
//That's all I can plan for now
class UnstructTopologyHolder;
class TopologyHolderStrategy
{
protected:
    UnstructTopologyHolder* d_hder = NULL;
public:
    TopologyHolderStrategy(UnstructTopologyHolder* hder);
    ~TopologyHolderStrategy();

    virtual int getNEQU() = 0;

    virtual void registerTopology() = 0;
};