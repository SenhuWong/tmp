#pragma once
#include "../toolBox/vector3d.h"
#include "TopologyHolderStrategy.h"
#include "Vankatakrishnan_Limiter.h"
#include "FluxStrategy.h"
// I am not ready to decide what's for name.
#include "SAMRAI/tbox/InputDatabase.h"
#include "TimeStrategy.h"


class Euler2D : public TopologyHolderStrategy
{
public:
    int d_NEQU = 4;
    int d_nmesh = -1;

private:
    // Const
    const double R = 287.0;
    const double Gamma = 1.4;
    const double PI = 3.141593;

    const double K = 5; // For limiter
    // Dimension
    const int d_dim = 2;
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

    // double fs_viscos_efficient;
    // double fs_Reynolds_number;

    // Nondimensional freeStream Variables,should be calculated right after the above is obtained.
    double fsnd_soundSpeed;
    double fsnd_density;
    double fsnd_velocity_components[3];
    double fsnd_pressure;
    double fsnd_Temperature;
    double fsnd_E;

    // FreeStream cnservatives, only 4 used in this Euler2D
    double fs_primVar[4];
    // Data storage for MPI communciation
    // Each mesh has only one block as requested by TIOGA;
    // Each block has a sendBUffer for each proc;
    double **commSendBuffer = NULL;
    double **commRecvBuffer = NULL;

    // int buff_size;
    //  Data storage for Cell Variables

    // Variables in CELLS  [blk][][cell]
    double ***U = NULL; // Conservatives
    GeomElements::vector3d<2, double> ***gradU = NULL;
    double ***UL = NULL;
    double ***UR = NULL;
    double ***Residual = NULL; // Sum of edge's fluxes
    double **Spectrum_cell = NULL;

    // Variables in EDGES  [equ][blk][edge]
    double ***U_edge = NULL;
    double ***Flux_edge = NULL;
    double **Spectrum_edge = NULL;

    // Used for storing Maximum and minimum neighbours for limiter
    // LimiterStrategy:
    Vankatakrishnan_Limiter *d_limiter = NULL;
    FluxStrategy *d_fluxComputer = NULL;
    // TimeStrategy *d_timeIntegrator = NULL;

    // double ***Wmax = NULL;
    // double ***Wmin = NULL;

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
    int cur_proc = -1;
    int num_proc = -1;
    UnstructTopologyHolder *getIntegrator()
    {
        return d_hder;
    }

    Euler2D(UnstructTopologyHolder *hder); // I will add an input_db to read the freestream variable value.

    void test_communication();

    void test_partialcomm();
    
    void test_unwantedSweep(int** UnwantedInd,int* nUnwantedInd);

    // Initialize Data to t0
    void initializeData();
    // Perform Communication to fill in the blank of buffer

    void ReconstructionInitial();

    void withinBlockCommunication();
    
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
    inline double ***getU()
    {
        return U;
    }

    inline double ***getUEdge()
    {
        return U_edge;
    }

    inline GeomElements::vector3d<2, double> ***getGradient()
    {
        return gradU;
    }

    inline double ***getUL()
    {
        return UL;
    }

    inline double ***getUR()
    {
        return UR;
    }

    inline double ***getFluxEdge()
    {
        return Flux_edge;
    }

    inline double ***getResidual()
    {
        return Residual;
    }

    
    int getNEQU()
    {
        return d_NEQU;
    }


    // This is intended to be called by TimeStrategy,should be an interface.
    void preprocessAdvance(int istage)
    {
        ReconstructionInitial();
        // Communicate between buffered Cells.
        // writeCellData("BeforeComm"+std::to_string(istage)+"stage", cur_proc, 0, W);
        withinBlockCommunication();
        //writeCellData("AfterComm" + std::to_string(istage) + "stage", cur_proc, 0, U);
        // Get conservative variables at edge

        // Output for debug
        // d_hder->blk2D[0].writeEdgeLeftRight("0000000",0,cur_proc);
        updateEdgeValues();

        solveGradient();
        solveSpectralRadius();
        // writeEdgeValue("AfterCompute"+std::to_string(istage)+"stage",Wi);
        // writeEdgeValue("Conservatives"+std::to_string(istage),Wi);
        //std::string filename = "BoundaryAtStage" + std::to_string(istage);
        //writeBoundaryInfo(filename);
        
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

    void cleaning(double ***W0)
    {
        for (int i = 0; i < d_hder->d_nmesh; i++)
        {
            for (auto iter = d_hder->recvBegin(i); iter != d_hder->recvEnd(i); iter++)
            {
                for (int j = 0; j < d_NEQU; j++)
                {
                    U[i][j][iter->d_localId] = 0;
                    W0[i][j][iter->d_localId] = 0;
                }
            }
        }
    }

    void postprocessAdvance()
    {
    }

    void SolveDiag(double*** de_diag, double** dt )
    {
        double OMEGAN = 5;
        for(int i = 0;i<d_nmesh;i++)
        {
            auto& curBlk = d_hder->blk2D[i];
            for(int k = 0;k<d_hder->nCells(i);k++)
            {
                auto& curCell = curBlk.d_localCells[k];
                double diag = curCell.volume()/dt[i][k] + 0.5*OMEGAN*Spectrum_cell[i][k];
                for(int j = 0;j<d_NEQU;j++)
                {
                    de_diag[i][j][k] = diag;
                }
            }
        }
    }

    void ConservativeParam2ConvectiveFlux(double* W, double* Fc, GeomElements::vector3d<2,double>& norm_vec)
    {
        GeomElements::vector3d<2,double> velocity(W[2]/W[0],W[3]/W[0]);
        double Vn = velocity.dot_product(norm_vec);
        double DatCell = W[0];
        double PatCell = (Gamma-1)*(W[1] - 0.5*DatCell*velocity.L2Square());
        double rhoEatCell = W[1];

        Fc[0] = W[0]*Vn;
        Fc[1] = (W[1] + PatCell)*Vn;
        
        Fc[2] = W[2]*Vn + PatCell*norm_vec[0];
        Fc[3] = W[3]*Vn + PatCell*norm_vec[1];
    }


    
    void SolveDF(int curMesh,int curCell,int curEdge,
                double***W_scratch ,double*** deltaW0,double* DF,
                int ForB,int cellOffset)
    {
        // if(cur_proc==0)

        // std::cout<<"curCell is "<<curCell<<'\n';
        auto& edg = d_hder->blk2D[curMesh].d_localEdges[curEdge];
        double Wp[d_NEQU];
        double W[d_NEQU];
        double W0[d_NEQU];
        double Fcp[d_NEQU];
        double Fc[d_NEQU];
        double DFc[d_NEQU];
        int lC = edg.lCInd();
        int rC = edg.rCInd();
        
        if(curCell + cellOffset==lC)
        {
            if((rC>=0 and rC<lC) && ForB)//If ForB is 1, then calculate when curCell(lC) >the other cell(rC)
            //If ForB is 0, then calculate when curCell(lC) < the other cell(rC)
            {
                if(rC >= d_hder->nCells(curMesh))
                {
                    rC = rC - d_hder->nCells(curMesh);
                }
                for(int j = 0;j<d_NEQU;j++)
                {
                    W0[j] = W_scratch[curMesh][j][rC];
                    Wp[j] =W0[j] + deltaW0[curMesh][j][rC];
                    W[j] = W0[j];
                }
                //Get convectiveFlux
                GeomElements::vector3d<2,double> norm = edg.normal_vector();
                ConservativeParam2ConvectiveFlux(W,Fc,norm);
                ConservativeParam2ConvectiveFlux(Wp,Fcp,norm);
                for(int j = 0;j<d_NEQU;j++)
                {
                    DFc[j] = Fcp[j] - Fc[j];
                }

                for(int j = 0;j<d_NEQU;j++)
                {
                    DF[j] = 0.5*(DFc[j])*edg.area();
                }
            }
            else
            {
                for(int j = 0;j<d_NEQU;j++)
                {
                    DF[j] = 0;
                }
            }
        }
        else if(curCell + cellOffset==rC)
        {
            if((lC<rC)&& ForB)
            {
                if(lC>=d_hder->nCells(curMesh))
                {
                    lC = lC - d_hder->nCells(curMesh);
                }
                for(int j = 0;j<d_NEQU;j++)
                {
                    W0[j] = W_scratch[curMesh][j][lC];
                    Wp[j] = W0[j] + deltaW0[curMesh][j][lC];
                    W[j] = W0[j];
                }
                //Get convectiveFlux
                GeomElements::vector3d<2,double> norm = edg.normal_vector();
                ConservativeParam2ConvectiveFlux(W,Fc,norm);
                ConservativeParam2ConvectiveFlux(Wp,Fcp,norm);
                for(int j = 0;j<d_NEQU;j++)
                {
                    DFc[j] =  Fc[j] - Fcp[j];
                }

                for(int j = 0;j<d_NEQU;j++)
                {
                    DF[j] = 0.5*(DFc[j])*edg.area();
                }
            }
            else
            {
                for(int j = 0;j<d_NEQU;j++)
                {
                    DF[j] = 0;
                }
            }
        }
        else
        {
            if(cur_proc==0)
            {
                std::cout<<lC<<","<<rC<<":"<<curCell<<","<<d_hder->nCells(curMesh)<<","<<cellOffset<<","<<ForB<<'\n';
            }
            throw std::runtime_error("Cur cell is not on either side,impossible\n");
        }
    }



    //The skip point in here are level 1 sendcells and level 1,2 recvcells.
    void SolveLocalForwardSweep(double*** diag,
                                double*** W_scratch,double*** deltaW1,
                                int** skip_point,int* nskip)
    {
        int ForB = 1;
        
        double dFi[d_NEQU];
        double dF[d_NEQU];
        for(int i = 0;i<d_nmesh;i++)
        {
            int cellIdOffset = 0;
            auto& curBlk = d_hder->blk2D[i];
            //Refer to Unstruct
            for(int k = 0;k<nskip[i]-1;k++)
            {
                
                for(int l = skip_point[i][k]+1;l<skip_point[i][k+1];l++)
                {
                    int cellId = l;
                    
                    auto& curCell = curBlk.d_localCells[cellId];
                    for(int j = 0;j<d_NEQU;j++)
                    {
                        dFi[j] = 0;
                    }
                    
                    for(int j = 0;j<curCell.edge_size();j++)
                    {
                        //std::cout<<cellId<<":"<<j<<'\n';
                        SolveDF(i,cellId,curCell.edgeInd(j),W_scratch,deltaW1,dF,ForB,cellIdOffset);
                        for(int m = 0;m<d_NEQU;m++)
                        {
                            dFi[m] += dF[m];
                        }
                    }
                    
                    for(int j = 0;j<d_NEQU;j++)
                    {
                        deltaW1[i][j][cellId] = (-Residual[i][j][cellId] - dFi[j])/diag[i][j][cellId];
                    }
                }
            }
        }
    }
    //The skip point in here are level 1 cells.
    void SolveBoundaryForwardSweep(double*** diag,
                                double*** W_scratch,double*** deltaW1,
                                int** skip_point,int * nskip)
    {
        int ForB = 1;
        double dFi[d_NEQU];
        double dF[d_NEQU];
        for(int i = 0;i<d_nmesh;i++)
        {
            int cellIdOffset = d_hder->nCells(i);
            auto& curBlk = d_hder->blk2D[i];
            //Refer to Unstruct
            for(int k = 0;k<nskip[i];k++)
            {
                int cellId = skip_point[i][k];
                auto& curCell = curBlk.d_localCells[cellId];
                for(int j = 0;j<d_NEQU;j++)
                {
                    dFi[j] = 0;
                }
                for(int j = 0;j<curCell.edge_size();j++)
                {
                    SolveDF(i,cellId,curCell.edgeInd(j),W_scratch,deltaW1,dF,ForB,cellIdOffset);
                    for(int m = 0;m<d_NEQU;m++)
                    {
                        dFi[m] += dF[m];
                    }
                }
                for(int j = 0;j<d_NEQU;j++)
                {
                    deltaW1[i][j][cellId] = (-Residual[i][j][cellId]-dFi[j])/diag[i][j][cellId];
;               }   
            }
        }
    }

    void SolveBoundaryBackwardSweep(double*** diag,
                                    double*** W_scratch,double*** deltaW1,double*** deltaW,
                                    int ** skip_point,int * nskip)
    {
        int ForB = 0;
        double dFi[d_NEQU];
        double dF[d_NEQU];
        for(int i = 0;i< d_nmesh;i++)
        {
            int cellIdOffset = d_hder->nCells(i);
            auto& curBlk = d_hder->blk2D[i];
            for(int k = 0;k<nskip[i];k++)
            {
                int cellId = skip_point[i][k];
                auto& curCell = curBlk.d_localCells[cellId];
                for(int j = 0;j<d_NEQU;j++)
                {
                    dFi[j] = 0;
                }
                for(int j = 0;j<curCell.edge_size();j++)
                {
                    SolveDF(i,cellId,curCell.edgeInd(j),W_scratch,deltaW,dF,ForB,cellIdOffset);
                    for(int m = 0;m<d_NEQU;m++)
                    {
                        dFi[m] += dF[m];
                    }
                }
                for(int j = 0;j < d_NEQU;j++)
                {
                    deltaW[i][j][cellId] = deltaW1[i][j][cellId] - dFi[j]/diag[i][j][cellId];
                }
            }
        }

    }

    void SolveLocalBackwardSweep(double*** diag,
                                double*** W_scratch,double*** deltaW1,double*** deltaW,
                                int** skip_point,int* nskip)
    {
        int ForB = 0;
        double dFi[d_NEQU];
        double dF[d_NEQU];
        for(int i = 0;i< d_nmesh;i++)
        {
            int cellIdOffset = 0;
            auto& curBlk = d_hder->blk2D[i];
            for(int k = 0;k<nskip[i]-1;k++)
            {
                for(int l = skip_point[i][k]+1;l<skip_point[i][k+1];l++)
                {
                    int cellId = l;
                    auto& curCell = curBlk.d_localCells[cellId];
                    for(int j =0;j<d_NEQU;j++)
                    {
                        dFi[j] = 0;
                    }
                    for(int j = 0;j<curCell.edge_size();j++)
                    {
                        SolveDF(i,cellId,curCell.edgeInd(j),W_scratch,deltaW,dF,ForB,cellIdOffset);
                        for(int m = 0;m<d_NEQU;m++)
                        {
                            dFi[m] += dF[m];
                        }
                    }
                    for(int j = 0;j<d_NEQU;j++)
                    {
                        deltaW[i][j][cellId] = deltaW1[i][j][cellId] - dFi[j]/diag[i][j][cellId];
                    }
                }
            }
        }
    }

    void SolveForwardSweep(double*** diag,
                        double*** W_scratch,double*** deltaW1)
    {
        int ForB = 1;
        double dFi[d_NEQU];
        double dF[d_NEQU];
        for(int i = 0;i<d_nmesh;i++)
        {
            auto& curBlk = d_hder->blk2D[i];
            for(int k = 0;k<d_hder->nCells(i);k++)
            {
                auto& curCell = curBlk.d_localCells[k];
                for(int j = 0;j<d_NEQU;j++)
                {
                    dFi[j] = 0;
                }
                for(int j = 0 ;j<curCell.edge_size();j++)
                {
                    SolveDF(i,k,curCell.edgeInd(j),W_scratch,deltaW1,dF,ForB,0);
                    for(int l = 0;l<d_NEQU;l++)
                    {
                        dFi[l] += dF[l];
                    }
                }
                for(int j = 0;j<d_NEQU;j++)
                {
                    deltaW1[i][j][k] = (-Residual[i][j][k]-dFi[j])/diag[i][j][k];
                }
            }
        }
    }

    void SolveBackwardSweep(double*** diag,
                            double*** W_scratch,double*** deltaW1,double*** deltaW)
    {
        int ForB = 0;
        double dFi[d_NEQU];
        double dF[d_NEQU];
        for(int i = 0;i<d_nmesh;i++)
        {
            auto& curBlk = d_hder->blk2D[i];
            for(int k = 0;k<d_hder->nCells(i);k++)
            {
                auto& curCell = curBlk.d_localCells[k];
                for(int j = 0;j<d_NEQU;j++)
                {
                    dFi[j] = 0;
                }
                for(int j = 0;j<curCell.edge_size();j++)
                {
                    SolveDF(i,k,curCell.edgeInd(j),W_scratch,deltaW,dF,ForB,0);
                    for(int l = 0;l<d_NEQU;l++)
                    {
                        dFi[l] += dF[l];
                    }
                }
                for(int j = 0;j<d_NEQU;j++)
                {
                    deltaW[i][j][k] = deltaW1[i][j][k] - dFi[j]/diag[i][j][k];
                }
            }
        }
    }

    void RemoteCellCommunication(double*** value);

    void NearCellCommunication(double*** value);

    
//These are all about writing
    void writeTiogaFormat(const std::string &filename, int proc, int meshTag, int *icelldata);

    void writeCellData(const std::string &filename, int proc, int meshTag, double ***dcelldata);
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
        fout << "fs_primVar[0] " << fs_primVar[0] << '\n';
        fout << "fs_primVar[1] " << fs_primVar[1] << '\n';
        fout << "fs_primVar[2] " << fs_primVar[2] << '\n';
        fout << "fs_primVar[3] " << fs_primVar[3] << '\n';
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
    void writeBoundaryInfo(std::string &filename)
    {

        for (int i = 0; i < d_nmesh; i++)
        {
            std::string local_filename = filename + std::to_string(i) + std::to_string(cur_proc);
            std::ofstream fout;
            fout.open(local_filename);
            auto &curBlk = d_hder->blk2D[i];
            for (int j = 0; j < d_hder->nEdges(i); j++)
            {
                auto &curEdge = curBlk.d_localEdges[j];
                int rC = curEdge.rCInd();
                int lC = curEdge.lCInd();
                if (rC == GeomElements::edge3d<2>::BoundaryType::WALL)
                {
                    fout << "Edge " << j << " Right is Wall:\n";
                    fout << "Edge's normal is " << curEdge.normal_vector()[0] << "," << curEdge.normal_vector()[1] << '\n';
                    fout << "Edge's center is " << curEdge.center()[0] << "," << curEdge.center()[1] << '\n';
                    fout << "Edge's area is " << curEdge.area() << '\n';
                    fout << U_edge[i][0][j] << "," << U_edge[i][1][j] << "," << U_edge[i][2][j] << "," << U_edge[i][3][j] << "\n";
                    fout << "\t left Cell is:\n";
                    fout << U[i][0][lC] << "," << U[i][1][lC] << "," << U[i][2][lC] << "," << U[i][3][lC] << "\n";
                }
            }
            for (int j = 0; j < d_hder->nEdges(i); j++)
            {
                auto &curEdge = curBlk.d_localEdges[j];
                int lC = curEdge.lCInd();
                int rC = curEdge.rCInd();
                if (rC == GeomElements::edge3d<2>::BoundaryType::FARFIELD)
                {
                    fout << "Edge " << j << " Right is FARFIELD:\n";
                    fout << "Edge's normal is " << curEdge.normal_vector()[0] << "," << curEdge.normal_vector()[1] << '\n';
                    fout << "Edge's center is " << curEdge.center()[0] << "," << curEdge.center()[1] << '\n';
                    fout << "Edge's area is " << curEdge.area() << '\n';
                    fout << U_edge[i][0][j] << "," << U_edge[i][1][j] << "," << U_edge[i][2][j] << "," << U_edge[i][3][j] << "\n";
                    fout << "\t left Cell is:\n";
                    fout << U[i][0][lC] << "," << U[i][1][lC] << "," << U[i][2][lC] << "," << U[i][3][lC] << "\n";
                    fout << "Far field var is:\n";
                    fout << fsnd_density << "," << fsnd_pressure << "," << fsnd_velocity_components[0] << "," << fsnd_velocity_components[1] << '\n';
                }
            }
            fout.close();
        }
    }

    void outPutCp(std::string& filename, int mesh_ind);
};
