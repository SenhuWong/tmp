#pragma once

#include "UnstructIntegrator.h"
#include "TopologyHolderStrategy.h"
#include "../toolBox/vector3d.h"

#include "Vankatakrishnan_Limiter.hpp"
#include "FluxStrategy.h"
#include "TimeStrategy.h"


class ViscousFlow2D : public TopologyHolderStrategy
{
private:
    // Const
    const double R = 287.0;
    const double Gamma = 1.4;
    const double PI = 3.141593;

    const double K = 5; // For limiter
    // FreeStream,this should be read from an input file
    double fs_mach = 0;
    double fs_AOA = 0;
    double fs_density = 0;
    double fs_pressure = 0;
    double fs_eigen_Length;
    // Derived freeStream
    double fs_soundSpeed;
    double fs_velocity_magnitude;
    double fs_velocity_components[2];
    double fs_E;
    // Not used yet
    double fs_Temperature;
    double fs_mu_Laminar;
    double fs_Reynolds_number;

    // Nondimensional freeStream Variables,should be calculated right after the above is obtained.
    double fsnd_pressure;
    double fsnd_density;
    double fsnd_soundSpeed;
    double fsnd_velocity_components[3];
    double fsnd_E;

    //Not used yet
    double fsnd_Temperature;
    double fsnd_mu_Laminar;

    // FreeStream cnservatives, only 4 used in this Euler2D
    double fsnd_primVar[4];

    double **commSendBuffer = NULL;
    double **commRecvBuffer = NULL;

    double **U = NULL; // Conservatives
    GeomElements::vector3d<2, double> **gradU = NULL;

    //Not used yet
    GeomElements::vector3d<2, double> **gradT = NULL;
    
    
    double ** mu = NULL;
    double ** Residual = NULL; // Sum of edge's fluxes
    double **Spectrum_cell = NULL;

    // Variables in EDGES  [equ][blk][edge]
    double ** U_edge = NULL;
    double ** mu_edge = NULL;
    double ** UL = NULL;
    double ** UR = NULL;
    double ** Flux_edge = NULL;
    double **Spectrum_edge = NULL;

    // Used for storing Maximum and minimum neighbours for limiter
    // LimiterStrategy:
    LimiterStrategy *d_limiter = NULL;
    FluxStrategy *d_fluxComputer = NULL;

    int **testing_flag = NULL;

private:
    // These should be called in the constructor given inputDatabase
    //  Allocate Data Storage
    void registerTopology();

    // Set freeStream Variable value
    void set_fs_variable(double Ma, double AOA, double density, double pressure, double eigenLen);

    // Update the nondim variables
    void update_nondim_variable();

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

    void AllCellCommunication(double ** value);
    
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

    inline void **getGradientPrimitive() override
    {
        return (void**)gradU;
    }

    inline void **getGradientT() override
    {
        return (void**)gradT;
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

    inline double **getMu() override
    {
        return mu;
    }

    inline double **getMuEdge() override
    {
        return mu_edge;

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
        
        // Get Left right value at edge
        updateEdgeLeftRightValues();
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

    

    void RemoteCellCommunication(double** value);

    void NearCellCommunication(double** value);

    
//These are all about writing
    void writeTiogaFormat(const std::string &filename, int proc, int meshTag, int *icelldata);

    void writeCellData(const std::string &filename, int proc, int meshTag, double **dcelldata);
    void writeCellData(const std::string &filename, int proc, int meshTag, int ***icelldata);

    void outPutNondim(const std::string &filename, int proc)
    {
        std::string localFilename = filename + std::to_string(proc);
        std::ofstream fout;
        fout.open(localFilename);
        fout << "fs_mach " << fs_mach << '\n';
        fout << "fs_AOA " << fs_AOA << '\n';
        fout << "fs_density " << fs_density << '\n';
        fout << "fs_pressure " << fs_pressure << '\n';
        fout << "fs_eigen_Length " << fs_eigen_Length << '\n';
        fout << "fs_soundSpeed " << fs_soundSpeed << '\n';
        fout << "fs_velocity_magnitude " << fs_velocity_magnitude << '\n';
        fout << "fs_velocity_component[0] " << fs_velocity_components[0] << '\n';
        fout << "fs_velocity_component[1] " << fs_velocity_components[1] << '\n';
        fout << "fs_E " << fs_E << '\n';
        fout << "fsnd_velocity_components[0] " << fsnd_velocity_components[0] << '\n';
        fout << "fsnd_velocity_components[1] " << fsnd_velocity_components[1] << '\n';
        fout << "fsnd_pressure " << fsnd_pressure << '\n';
        fout << "fsnd_density " << fsnd_density << '\n';
        fout << "fsnd_E " << fsnd_E << '\n';
        fout << "fsnd_primVar[0] " << fsnd_primVar[0] << '\n';
        fout << "fsnd_primVar[1] " << fsnd_primVar[1] << '\n';
        fout << "fsnd_primVar[2] " << fsnd_primVar[2] << '\n';
        fout << "fsnd_primVar[3] " << fsnd_primVar[3] << '\n';
        fout.close();
    }
    void writeUnacceptableCell(const std::string &filename, int proc, int meshTag, int **dcelldata);
    void updateEdgeConservatives_butWithOutput();
    void writeEdgeValue(const std::string &filename, double ***Wedge)
    {
        for (int i = 0; i < d_hder->d_nmesh; i++)
        {
            std::string localFilename = filename + std::to_string(i) + std::to_string(cur_proc);
            std::ofstream fout;
            fout.open(localFilename);
            for (int j = 0; j < d_hder->nEdges(i); j++)
            {
                auto &curEdge = d_hder->blk2D[i].d_localEdges[j];
                fout << j << " " << curEdge.lCInd() << " " << curEdge.rCInd() << '\n';
                for (int l = 0; l < d_NEQU; l++)
                {
                    fout << Wedge[i][l][j] << " ";
                }
                fout << '\n';
            }
            fout.close();
        }
    }

    void outPutCp(std::string& filename, int mesh_ind);
};
