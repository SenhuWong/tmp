#pragma once

#include "UnstructIntegrator.h"
#include "TopologyHolderStrategy.h"
#include "../toolBox/vector3d.h"
#include "BoundaryStrategy.h"
#include "Vankatakrishnan_Limiter.hpp"
#include "FluxStrategy.h"
#include "TimeStrategy.h"
#include "LaminarFluxStrategy.h"

class ViscousFlow2D : public TopologyHolderStrategy
{
private:
    // Variables in CELLS  [blk][][cell]
    double **U = NULL; // Conservatives

    //GradU is to be used by both LimiterStrategy and LaminarFluxStrategy, so I leave it in Euler2D.
    GeomElements::vector3d<2, double> **gradU = NULL;
    GeomElements::vector3d<2, double> **gradT = NULL;
    
    double **Residual = NULL; // Sum of edge's fluxes
    double **Spectrum_cell_c = NULL;

    // Variables in EDGES  [equ][blk][edge]

    //U_edge is used by both LImiterStrategy and LaminarFluxStrategy, so I leave it here in Euler2D.
    double ** U_edge = NULL;


    // LimiterStrategy:
    FluxStrategy *d_cfluxComputer = NULL;
    LaminarFluxStrategy *d_vfluxComputer = NULL;
    
    BoundaryStrategy * d_boundaryHandler = NULL;

    //Not used yet
    
    
    
    
    double **Spectrum_cell_c = NULL;
    double **Spectrum_cell_v = NULL;

    double ** Residual = NULL; // Sum of edge's fluxes
    // Variables in EDGES  [equ][blk][edge]
    

private:
    // These should be called in the constructor given inputDatabase
    //  Allocate Data Storage
    void registerTopology();

public:
    UnstructTopologyHolder *getIntegrator()
    {
        return d_hder;
    }

    ViscousFlow2D(UnstructTopologyHolder *hder); // I will add an input_db to read the freestream variable value.

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
            auto& curBlk = d_hder->blk2D[i];
            for(int k =0;k<d_hder->nCells(i);k++)
            {
                auto& curCell = curBlk.d_localCells[k];
                dt[i][k] = CFL * curCell.volume()/(Spectrum_cell_c[i][k]+Spectrum_cell_v[i][k]);
            }
        }
    }

    double getMuOverPr() override
    {
        return d_vfluxComputer->getMuOverPr();
    }


    // Update flux with the corresponding left and right value
    void updateFlux();
    inline double **getU() override
    {
        return U;
    }

    inline double**getSpectrumViscour() override
    {
        return Spectrum_cell_v;
    }

    inline void **getGradientPrimitive() override
    {
        return (void**)gradU;
    }

    inline void **getGradientT() override
    {
        return (void**)gradT;
    }

    inline double **getResidual() override
    {
        return Residual;
    }

    inline double **getSpectrumConvective() override
    {
        return Spectrum_cell_c;
    }

    inline double **getMu() override
    {
        return d_vfluxComputer->getMu();
    }

    inline double **getMuEdge() override
    {
        return d_vfluxComputer->getMuEdge();
    }

    // This is intended to be called by TimeStrategy,should be an interface.
    void preprocessAdvance(int istage)
    {
        ReconstructionInitial();
        // Communicate between buffered Cells.
        AllCellCommunication(getU());
        // For each cell compute mu.

        // Get conservative variables at edge
        updateEdgeValues();

        solveGradient();

        //SpectralRadius need to add viscos terms, but I am not sure how.
        solveSpectralRadius();
        

        // writeEdgeValue("LeftValue"+std::to_string(istage)+"stage",WL);
        // writeEdgeValue("RightValue"+std::to_string(istage)+"stage",WR);

        // Get flux from left right value(or boundary)
        // And compute fluxsum
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
        //This should do stuff like:
        //compute mu with updated 
    }

   
    void outPutCp(std::string& filename, int mesh_ind);
};
