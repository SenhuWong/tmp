#pragma once
class Euler2D;

#include "UnstructIntegrator.h"
#define LIMITER_BIGVALUE 10.0;
//#include "Euler2D.h"
// A Valid limiter, given the TopologyHolder and TopologyHolderStrategy, it should be able to
// (1).compute for each cell its limiter

class LimiterStrategy
{
public:
    int d_nmesh = -1;
    int d_dim = 0;
    int d_nequ = 0;
    double ***mesh_var_cell_limit = NULL;
    UnstructTopologyHolder *d_hder = NULL;
    Euler2D *d_hder_strategy = NULL;
#define DEBUG
#ifdef DEBUG
    int num_proc = -1;
    int cur_proc = -1;
#endif // DEBUG
public:
#ifdef DEBUG
    void setComm(int cur, int num)
    {
        cur_proc = cur;
        num_proc = num;
    }
#endif // DEBUG

    LimiterStrategy(UnstructTopologyHolder *hder, Euler2D *hder_strategy);
    ~LimiterStrategy();
public:
    // Provided Left and Right Consrevative Storage, compute the Left and Right increment as limiter is computed.
    virtual void computeLimiter() = 0;
    // Get Limiter for cell given meshInd,equInd,cellInd.
    virtual double getLimiter(int meshInd, int equInd, int cellInd) = 0;

    virtual void init() = 0;
};