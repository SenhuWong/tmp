#pragma once

#include "UnstructIntegrator.h"
#include "TopologyHolderStrategy.h"
#include "../toolBox/vector3d.h"
#include "BoundaryStrategy.h"
#include "Vankatakrishnan_Limiter.hpp"
#include "FluxStrategy.h"
#include "TimeIntegrator.h"
#include "LaminarFluxStrategy.h"
#include "LUSGS_Integrator.h"
#include "SSTKWStrategy.h"


class ViscousFlow2D : public TopologyHolderStrategy
{
private:
    double **U = NULL;
    GeomElements::vector3d<2, double> **gradU = NULL;
    GeomElements::vector3d<2, double> **gradT = NULL;
    double ** U_edge = NULL;

    FluxStrategy *d_cfluxComputer = NULL;
    LaminarFluxStrategy *d_vfluxComputer = NULL;
    BoundaryStrategy * d_boundaryHandler = NULL;

    

    double **Spectrum_cell_c = NULL;
    double **Spectrum_cell_v = NULL;

    double ** Residual = NULL;
    
private:
    void registerTopology();
    int getNEquation() override
    {
        return 4;
    }

public:
    LUSGSIntegrator* turbModelTimeStrategy = NULL;
    SSTkomegaModel<2>* d_turbulentModel = NULL;

    //Interfaces from TopologyHolderStrategy:

    ViscousFlow2D()
    {

    }
    void initializeSubStrategy() override;
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

    void solveGradient(int istage);

    void solveSpectralRadius();


    

    void updateFlux(int istage);

    inline bool isInvicid() override
    {
        return false;
    }

    //Mandatory Interfaces(Getters): 
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
    //Override Interfaces(Getters):
    inline void **getGradientT() override
    {
        return (void**)gradT;
    }

    inline double**getSpectrumViscous() override
    {
        return Spectrum_cell_v;
    }

    inline double getMuCell(int curMesh,int curCell) override
    {
        if(isLaminar())
        {
            return d_vfluxComputer->getMuCell(curMesh,curCell); 
        }
        return d_vfluxComputer->getMuCell(curMesh,curCell) + d_turbulentModel->getMuCell(curMesh,curCell);
        
    }

    inline double getMuEdge(int curMesh,int curEdge) override
    {
        if(isLaminar())
        {
            return d_vfluxComputer->getMuEdge(curMesh,curEdge);
        }
        return d_vfluxComputer->getMuEdge(curMesh,curEdge) + d_turbulentModel->getMuEdge(curMesh,curEdge);

    }

    double getMuOverPrEdge(int curMesh,int curEdge) override
    {
        if(isLaminar())
        {
            return d_vfluxComputer->getMuOverPrEdge(curMesh,curEdge);
        }
        return d_vfluxComputer->getMuOverPrEdge(curMesh,curEdge) + d_turbulentModel->getMuOverPrEdge(curMesh,curEdge);
    }

    double getMuOverPrCell(int curMesh,int curCell) override
    {
        if(isLaminar())
        {
            return d_vfluxComputer->getMuOverPrCell(curMesh,curCell);
        }
        return d_vfluxComputer->getMuOverPrCell(curMesh,curCell) + d_turbulentModel->getMuOverPrCell(curMesh,curCell);
    }

    bool isLaminar()
    {
        return true;
    }


    //Mandatory Interfaces(Requested by timeStrategy)
    void preprocessAdvance(int istage)
    {
        // Communicate between buffered Cells.
        AllCellCommunication(getU());

        // Reconstruction Init
        ReconstructionInitial();

        updateEdgeValues();

        solveGradient(istage);

        d_vfluxComputer->computeMu();
        if(!isLaminar())
        {
            d_turbulentModel->SolveTurbulenceEquation();
        }
        //SpectralRadius need to add viscos terms, but I am not sure how.
        solveSpectralRadius();
        

        // writeEdgeValue("LeftValue"+std::to_string(istage)+"stage",WL);
        // writeEdgeValue("RightValue"+std::to_string(istage)+"stage",WR);

        // Get flux from left right value(or boundary)
        // And compute fluxsum
        updateFlux(istage);
    }

    void preprocessAdvanceSerial(int istage)
    {
        ReconstructionInitial();

        updateEdgeValues();

        solveGradient(istage);

        d_vfluxComputer->computeMu();

        //SpectralRadius need to add viscos terms, but I am not sure how.
        solveSpectralRadius();
        

        // writeEdgeValue("LeftValue"+std::to_string(istage)+"stage",WL);
        // writeEdgeValue("RightValue"+std::to_string(istage)+"stage",WR);

        // Get flux from left right value(or boundary)
        // And compute fluxsum
        updateFlux(istage);
    }

    void postprocessAdvanceSerial(int iStage)
    {
        if(!isLaminar())
        {
             std::cout<<"chao guoran!!!1\n";
            turbModelTimeStrategy->singleStepSerial(0);
        }
    }
    
    void SolveTime(double CFL)
    {
        for(int i = 0;i<d_nmesh;i++)
        {
            auto& curBlk = d_hder->blk2D[i];
            for(int k =0;k<d_hder->nCells(i);k++)
            {
                auto& curCell = curBlk.d_localCells[k];
                if(!isLaminar())
                {
                    d_stableDt[i][k] = CFL * curCell.volume()/(Spectrum_cell_c[i][k]+4.0*Spectrum_cell_v[i][k]);
                }
                else
                {
                    d_stableDt[i][k] = CFL * curCell.volume()/(Spectrum_cell_c[i][k]);
                }
                
            }
        }
    }

    void postprocessAdvance(int iStage)
    {
        if(!isLaminar())
        turbModelTimeStrategy->singleStep(0);
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

    
    void outPutCp(std::string& filename, int mesh_ind, int step);
};
