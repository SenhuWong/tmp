#pragma once
#include "../toolBox/vector3d.h"
#include "TopologyHolderStrategy.h"
#include "Vankatakrishnan_Limiter.hpp"
#include "FluxStrategy.h"
// I am not ready to decide what's for name.
#include "SAMRAI/tbox/InputDatabase.h"
#include "TimeStrategy.h"
#include "UnstructIntegrator.h"
#include "BoundaryStrategy.h"
class Euler2D : public TopologyHolderStrategy
{
private:
    // Variables in CELLS  [blk][][cell]
    double **U = NULL; // Conservatives

    //GradU is to be used by both LimiterStrategy and LaminarFluxStrategy, so I leave it in Euler2D.
    GeomElements::vector3d<2, double> **gradU = NULL;
    
    double ** Residual = NULL; // Sum of edge's fluxes
    double **Spectrum_cell_c = NULL;

    // Variables in EDGES  [equ][blk][edge]

    //U_edge is used by both LImiterStrategy and LaminarFluxStrategy, so I leave it here in Euler2D.
    double ** U_edge = NULL;


    // LimiterStrategy:
    FluxStrategy *d_fluxComputer = NULL;
    BoundaryStrategy * d_boundaryHandler = NULL;
    

private:
    // These should be called in the constructor given inputDatabase
    //  Allocate Data Storage
    void registerTopology();

public:

    Euler2D(UnstructTopologyHolder *hder); // I will add an input_db to read the freestream variable value.

    // Initialize Data to t0
    void initializeData();
    // Perform Communication to fill in the blank of buffer

    void ReconstructionInitial();
    
    /*
    **  These following three methods are from solveReconstruction()'s parts
    **  
    */
    // Boundary and normal edge's conservative update
    void updateEdgeValues();

    void solveGradient();

    void solveSpectralRadius();

    void SolveTime(double** dt, double CFL)
    {
        for(int i = 0;i<d_nmesh;i++)
        {
            for(int j =0;j<d_hder->nCells(i);j++)
            {
                auto& curCell = d_hder->blk2D[i].d_localCells[j];
                dt[i][j] = CFL * curCell.volume()/Spectrum_cell_c[i][j];
            }
        }
    }

    // Update flux with the corresponding left and right value
    void updateFlux();
    inline double **getU() override
    {
        return U;
    }

    inline double **getUEdge() override
    {
        return U_edge;
    }

    inline void **getGradientPrimitive() override
    {
        return (void**)gradU;
    }

    

    inline double **getResidual() override
    {
        return Residual;
    }

    inline double **getSpectrumConvective() override
    {
        return Spectrum_cell_c;
    }

    inline bool isInvicid() override
    {
        return true;
    }


    // This is intended to be called by TimeStrategy,should be an interface.
    void preprocessAdvance(int istage)
    {
        ReconstructionInitial();

        // Communicate between buffered Cells.
        AllCellCommunication(getU());

        // Get conservative variables at edge
        updateEdgeValues();

        solveGradient();

        solveSpectralRadius();

        // Compute fluxsum
        updateFlux();
    }

    void preprocessAdvanceSerial(int istage)
    {
        ReconstructionInitial();
        updateEdgeValues();
        solveGradient();
        solveSpectralRadius();
        updateFlux();
    }
    
    void test_MultipleCommunication();

    void cleaning(double **W0)
    {
        for (int i = 0; i < d_hder->d_nmesh; i++)
        {
            for (auto iter = d_hder->recvBegin(i); iter != d_hder->recvEnd(i); iter++)
            {
                for (int j = 0; j < d_NEQU; j++)
                {
                    U[i][j+d_NEQU*iter->d_localId] = 0;
                    W0[i][j+d_NEQU*iter->d_localId] = 0;
                }
            }
        }
    }

    void postprocessAdvance()
    {
    }


    
//These are all about writing

    void outPutCp(std::string& filename, int mesh_ind);
};
