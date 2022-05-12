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
    // FreeStream,this should be read from an input file
    
    

    // int buff_size;
    //  Data storage for Cell Variables

    // Variables in CELLS  [blk][][cell]
    double **U = NULL; // Conservatives
    GeomElements::vector3d<2, double> **gradU = NULL;
    double ** UL = NULL;
    double ** UR = NULL;
    double ** Residual = NULL; // Sum of edge's fluxes
    double **Spectrum_cell = NULL;

    // Variables in EDGES  [equ][blk][edge]
    double ** U_edge = NULL;
    double ** Flux_edge = NULL;
    double **Spectrum_edge = NULL;

    // Used for storing Maximum and minimum neighbours for limiter
    // LimiterStrategy:
    LimiterStrategy *d_limiter = NULL;
    FluxStrategy *d_fluxComputer = NULL;
    BoundaryStrategy * d_boundaryHandler = NULL;
    

private:
    // These should be called in the constructor given inputDatabase
    //  Allocate Data Storage
    void registerTopology();

public:
    UnstructTopologyHolder *getIntegrator()
    {
        return d_hder;
    }

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
                dt[i][j] = CFL * curCell.volume()/Spectrum_cell[i][j];
            }
        }
    }

    
    // Gradient computation for each cell and LeftRight Value Computation for each edge
    void updateEdgeLeftRightValues();
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

    inline double **getUL() override
    {
        return UL;
    }

    inline double **getUR() override
    {
        return UR;
    }

    inline double **getFluxEdge() override
    {
        return Flux_edge;
    }

    inline double **getResidual() override
    {
        return Residual;
    }

    inline double **getSpectrum() override
    {
        return Spectrum_cell;
    }

    inline bool isInvicid() override
    {
        return true;
    }

    inline void getFreeStreamVars(double* out)
    {
        for(int i = 0;i <d_NEQU;i++)
        {
            out[i] = fs_primVar[i];
        }
        out[d_NEQU] = fsnd_soundSpeed;
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

        // Get Left right value at edge
        updateEdgeLeftRightValues();

        // Compute fluxsum
        updateFlux();
    }

    void preprocessAdvanceSerial(int istage)
    {
        ReconstructionInitial();
        updateEdgeValues();
        solveGradient();
        solveSpectralRadius();
        updateEdgeLeftRightValues();
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
