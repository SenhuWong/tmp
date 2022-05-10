#pragma once
// This class should provide an interface
#include <stdio.h>
class Euler2D;
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
    UnstructTopologyHolder *d_hder = NULL;
    Euler2D *d_hder_strategy = NULL;
    TimeStrategy(UnstructTopologyHolder *hder, Euler2D *hder_strategy);
    //     : d_nmesh(hder->d_nmesh), d_NEQU(hder_strategy->d_NEQU), d_dim(hder->d_dim), d_hder(hder), d_hder_strategy(hder_strategy)
    // {

    // }
public:
    virtual void initialize() = 0;
    //Things need to be done before making a time step forward
    //virtual void preprocessUpdate() = 0;
    virtual void singleStep(int curStage) = 0;
    //virtual void advance() = 0;
    //virtual void postprocessUpdate() = 0;
    virtual bool converged() = 0;

};