#pragma once
#include "FluxStrategy.h"
#include "LimiterStrategy.h"
#include <iostream>
#include "UnstructIntegrator.h"
#include "toolBox/edge3d_int.h"
// Given 4 conservatives, calculate the pressure and sound speed
struct FlowParamHolder
{
    double density;
    double pressure;
    double soundSpeed;
    double velocity[3];
};

class HLLCFluxStrategy : public FluxStrategy
{
public:
    double ** UL = nullptr;
    double ** UR = nullptr;
    LimiterStrategy* d_limiter = nullptr;
    HLLCFluxStrategy(UnstructTopologyHolder *hder, TopologyHolderStrategy *hder_strategy, LimiterStrategy* limiter)
        : FluxStrategy(hder, hder_strategy)
    {
        d_limiter = limiter;
        UL = new double *[d_nmesh];
        UR = new double *[d_nmesh];
        for(int i = 0;i<d_nmesh;i++)
        {
            UL[i] = new double [d_NEQU*d_hder->nEdges(i)];
            UR[i] = new double [d_NEQU*d_hder->nEdges(i)];
        }
    }

    void updateLeftRight()
    {
        d_limiter->init();
        const double EPSILON = 1.0e-5;
        // Compute Limiter for each cell
        // This should be done by a Limiter class
        d_limiter->computeLimiter(UL,UR);
        double** U = d_hder_strategy->getU();
        double** U_edge = d_hder_strategy->getUEdge();
        bool hasTrouble = false;
        // Compute left and right value for each Edge
        for (int i = 0; i < d_nmesh; i++)
        {
            auto &curBlk = d_hder->blk2D[i];
            for (int k = 0; k < curBlk.nEdges(); k++)
            {

                auto &curEdge = curBlk.d_localEdges[k];
                int lc = curEdge.lCInd();
                int rc = curEdge.rCInd();

                if (lc >= 0 and rc >= 0)
                {
                    for (int j = 0; j < d_NEQU; j++)
                    {
                        UL[i][j+d_NEQU*k] = UL[i][j+d_NEQU*k] * d_limiter->getLimiter(i, j, lc) + U[i][j+d_NEQU*lc];
                    }
                    for (int j = 0; j < d_NEQU; j++)
                    {
                        UR[i][j+d_NEQU*k] = UR[i][j+d_NEQU*k] * d_limiter->getLimiter(i, j, rc) + U[i][j+d_NEQU*rc];
                    }
                }
                else
                {
                    if (lc < 0)
                    {
                        throw std::runtime_error("Left Cell can't be NULL\n");
                    }
                    if (rc == GeomElements::edge3d<2>::BoundaryType::WALL) // Wall, set to be the cell variables
                    {
                        for (int j = 0; j < d_NEQU; j++)
                        {
                            UL[i][j+d_NEQU*k] = UR[i][j+d_NEQU*k] = U_edge[i][j+d_NEQU*k];
                        }
                    }
                    else if (rc == GeomElements::edge3d<2>::BoundaryType::FARFIELD) // Far_field, set to be the face variables
                    {
                        for (int j = 0; j < d_NEQU; j++)
                        {
                            UL[i][j+d_NEQU*k] = UR[i][j+d_NEQU*k] = U_edge[i][j+d_NEQU*k];
                        }
                    }
                    else
                    {
                        std::cout << lc << "," << rc << '\n';
                        std::cout << "SOmething wrong with left and right cell\n";
                        std::cin.get();
                    }
                }
            }
        }
        if (hasTrouble)
        {
            throw std::runtime_error("Limiter Error");
        }
        }
        void computeFlux();
}; 