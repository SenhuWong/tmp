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

public:

    //Interfaces from TopologyHolderStrategy:


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

    void solveSpectralRadius()
    {
        double Spectrum_edge;
        GeomElements::vector3d<2, double> velocity_edge;
        for (int i = 0; i < d_nmesh; i++)
        {
            auto &curBlk = d_hder->blk2D[i];
            for (int k = 0; k < d_hder->nEdges(i); k++)
            {
                auto &curEdge = curBlk.d_localEdges[k];
                velocity_edge[0] = U_edge[i][2+d_NEQU*k];
                velocity_edge[1] = U_edge[i][3+d_NEQU*k];
                double c = sqrt(Gamma * U_edge[i][1+d_NEQU*k] / U_edge[i][0+d_NEQU*k]);
                double Vn = velocity_edge.dot_product(curEdge.normal_vector());
                int lC = curEdge.lCInd();
                int rC = curEdge.rCInd();
                Spectrum_edge = (std::abs(Vn) + c) * curEdge.area();
                Spectrum_cell_c[i][lC] += Spectrum_edge;
                if (rC >= 0)
                {
                    Spectrum_cell_c[i][rC] += Spectrum_edge;
                }
                auto& leftCell = curBlk.d_localCells[lC];
                Spectrum_edge = std::max<double>(4/(3.0*U_edge[i][0+d_NEQU*k]),Gamma/U_edge[i][0+d_NEQU*k])*d_vfluxComputer->getMuOverPrEdge(i,k)*curEdge.area()*curEdge.area();
                Spectrum_cell_v[i][lC] += Spectrum_edge / leftCell.volume();
                if(rC >= 0)
                {
                    auto& rightCell = curBlk.d_localCells[rC];
                    Spectrum_cell_v[i][rC] += Spectrum_edge /  rightCell.volume();
                }
            }
        }
    }


    double getMuOverPrEdge(int curMesh,int curEdge) override
    {
        return d_vfluxComputer->getMuOverPrEdge(curMesh,curEdge);
    }

    double getMuOverPrCell(int curMesh,int curCell) override
    {
        return d_vfluxComputer->getMuOverPrCell(curMesh,curCell);
    }

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

    inline double **getMu() override
    {
        return d_vfluxComputer->getMu();
    }

    inline double **getMuEdge() override
    {
        return d_vfluxComputer->getMuEdge();
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
        return;
        d_vfluxComputer->computeMu();
        writeCellData("PrimVar"+std::to_string(istage),cur_proc,0,U);
        writeRawCellData("mu_cell"+std::to_string(istage),cur_proc,0,getMu());
        ReconstructionInitial();
        updateEdgeValues();
        solveGradient(istage);
        solveSpectralRadius();
        writeRawCellData("spectrum_c"+std::to_string(istage),cur_proc,0,Spectrum_cell_c);
        writeRawCellData("spectrum_v"+std::to_string(istage),cur_proc,0,Spectrum_cell_v);
        updateFlux(istage);
    }
    
    void SolveTime(double** dt, double CFL)
    {
        for(int i = 0;i<d_nmesh;i++)
        {
            auto& curBlk = d_hder->blk2D[i];
            for(int k =0;k<d_hder->nCells(i);k++)
            {
                auto& curCell = curBlk.d_localCells[k];
                dt[i][k] = CFL * curCell.volume()/(Spectrum_cell_c[i][k]+4.0*Spectrum_cell_v[i][k]);
            }
        }
    }

    void postprocessAdvance()
    {
        //Do nothing;
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

    
    void outPutCp(std::string& filename, int mesh_ind);
};
