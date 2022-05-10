
#pragma once
#define LIMITER_K 5
// This should be a subclass for limiter .

// Looks like I have to write limiter one for unstruct and another for struct.
// Or if all limiter want to know max and min, we could implement the filling
// of max and min in Strategy.
// And after that is filled ,Wmax and Wmin could be directly used to store
// delta1,max and delta1,min.
// For every edge that we want to compute limiter, we need to pass in the delta2
#include "LimiterStrategy.h"

class Vankatakrishnan_Limiter : public LimiterStrategy
{
private:
    double*** Umax=NULL;
    double*** Umin=NULL;

public:
    Vankatakrishnan_Limiter(UnstructTopologyHolder *hder, Euler2D *hder_strategy);

    void computeLimiter(); // cell_ind starts from 0
    
    double getLimiter(int meshInd, int equInd, int cellInd);

    void write_Lmitera();

    void computeLeftRight();

    void init();

};