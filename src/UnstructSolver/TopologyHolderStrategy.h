#pragma once
//This should be an abstract base class for Euler Equation or N-S Equation advancing.
//We need to specify the following routines to perform data management on TopologyHolder:
//1.Allocate data storage for Conservative Variables' Data
//2.Allocate data storage for Flux Variables' Data
//3.Estimate and allocate buffer storage for cross-process communication.
//4.Register subStrategies for Flux computation(by passing left right values or something else) and so on.
//That's all I can plan for now
#include <string>
#include "toolBox/vector3d.h"
template<int ndim>
class SSTkomegaModel;

class UnstructTopologyHolder;
class TopologyHolderStrategy
{
public:
    int d_NEQU = -1;
    int d_nmesh = -1;
    int d_dim = -1;
    UnstructTopologyHolder* d_hder = nullptr;

//Below are all the freeStream variables.
    // Const
    const double R = 287.0;
    const double Gamma = 1.4;
    const double PI = 3.141593;    
    //Input
    double fs_mach = 0;
    double fs_AOA = 0;
    double fs_density = 0;
    double fs_pressure = 0;
    double fs_eigen_Length = 0;
    // Derived freeStream
    double fs_soundSpeed;
    double fs_velocity_magnitude;
    double fs_velocity_components[3];
    double fs_Temperature;

    double fs_Pr = 0.72;
    double fs_Prt = 0.9;
    double fs_mu;
    double fs_Re;

    // Nondimensional freeStream Variables,should be calculated right after the above is obtained.
    double fsnd_density;
    double fsnd_pressure;
    double fsnd_soundSpeed;
    double fsnd_velocity_components[3];
    double fsnd_Temperature;
    double fs_primVar[5];

    //It's requested by both Flow and TurbulenceModel.
    double ** d_stableDt = NULL;

    

    // Data storage for MPI communciation
    // Each mesh has only one block as requested by TIOGA;
    // Each block has a sendBUffer for each proc;
    double **commSendBuffer = NULL;
    double **commRecvBuffer = NULL;

    // For debug.
    int **testing_flag = NULL;
    // For Tioga's iblank and q:
    int **d_iblank_cell = NULL;
    
friend class BoundaryStrategy;
friend class LUSGSIntegrator;
friend class LaminarFluxStrategy;

friend class SSTkomegaModel<2>;
    
public:
    int cur_proc = -1;
    int num_proc = -1;
public:
    TopologyHolderStrategy()
    {

    }
    void initializeStrategy(UnstructTopologyHolder* hder);
    virtual void initializeSubStrategy()
    {
        
    }
    TopologyHolderStrategy(UnstructTopologyHolder* hder, int nequ);
    ~TopologyHolderStrategy();

    void setNEquation(int nequ)
    {
        d_NEQU = nequ;
    }
    virtual int getNEquation()
    {
        return d_NEQU;
    }
    void setNMesh(int nmesh)
    {
        d_nmesh = nmesh;
    }
    int getNMesh()
    {
        return d_nmesh;
    }
    void setNDim(int ndim)
    {
        d_dim = ndim;
    }
    int getNDim()
    {
        return d_dim;
    }

    // Set freeStream Variable value
    void set_fs_variable(double Ma, double AOA, double density, double pressure, double eigenLen);
    void set_fs_variable2(double Ma, double AOA, double density, double pressure, double Re);
    void set_fs_variable(TopologyHolderStrategy* another)
    {

        std::cout<<"Another fs_mach"<<another->fs_mach<<'\n';
        fs_mach = another->fs_mach;
        std::cout<<"Now my mach is "<<fs_mach<<'\n';


        fs_AOA = another->fs_AOA;
        
        std::cout<<"Another fs_density"<<another->fs_density<<'\n';
        fs_density = another->fs_density;
        std::cout<<"Now my density is "<<fs_density<<'\n';
        fs_pressure = another->fs_pressure;

        fs_eigen_Length = another->fs_eigen_Length;
        fs_soundSpeed = another->fs_soundSpeed;
        fs_velocity_magnitude = another->fs_velocity_magnitude;
        for(int i = 0;i<3;i++)
        {
            fs_velocity_components[i] = another->fs_velocity_components[i];
            fsnd_velocity_components[i] = another->fsnd_velocity_components[i];
        }
        fs_Temperature = another->fs_Temperature;      
        fs_Pr = 0.72;
        fs_Prt = 0.9;
        fs_mu = another->fs_mu;
        fs_Re = another->fs_Re;
        fsnd_density = another->fsnd_density;
        fsnd_pressure = another->fsnd_pressure;
        fsnd_soundSpeed = another->fsnd_soundSpeed;
        
        fsnd_Temperature = another->fsnd_Temperature;
        for(int i = 0;i<5;i++)
        {
            fs_primVar[i] = another->fs_primVar[i];
        }
    }
    void set_fs_variable3(double Ma, double AOA, double temperature, double pressure, double Re)
    {
        fs_mach = Ma;
        std::cout<<"fs_mach "<<fs_mach<<'\n';
        fs_AOA = PI * AOA/180;
        std::cout<<"fs_AOA "<<fs_AOA<<'\n';
        fs_Temperature=temperature;
        std::cout<<"fs_Temperature "<<fs_Temperature<<'\n';
        //fs_pressure = pressure;
        //std::cout<<"fs_pressure "<<fs_pressure<<'\n';
        fs_Re = Re;
        std::cout<<"fs_Re "<<fs_Re<<'\n';

        fs_soundSpeed = sqrtf64(Gamma*R*fs_Temperature);
        std::cout<<"fs_soundSpeed "<<fs_soundSpeed<<'\n';
        fs_velocity_magnitude = fs_soundSpeed*fs_mach;
        std::cout<<"fs_velocity_magnitude "<<fs_velocity_magnitude<<'\n';
        fs_mu = 1.45e-6*pow(fs_Temperature,1.5)/(fs_Temperature+110);
        std::cout<<"fs_mu "<<fs_mu<<'\n';
        fs_density=fs_mu*fs_Re/(fs_velocity_magnitude*1.0);
        std::cout<<"fs_density "<<fs_density<<'\n';
        fs_pressure=fs_density*R*fs_Temperature;
        std::cout<<"fs_pressure "<<fs_pressure<<'\n';



        fs_velocity_components[0] = fs_velocity_magnitude * cosf64(fs_AOA);
        fs_velocity_components[1] = fs_velocity_magnitude * sinf64(fs_AOA);
        std::cout<<fs_velocity_components[0]<<" "<<fs_velocity_components[1]<<'\n';
        fs_eigen_Length = fs_Re*fs_mu/(fs_density*fs_velocity_magnitude);
        std::cout<<"fs_eigen_Length "<<fs_eigen_Length<<'\n';
        update_nondim_variable();
    }
    // Update the nondim variables
    void update_nondim_variable();

    void AllCellCommunication(double** value);

    void RemoteCellCommunication(double** value);

    void NearCellCommunication(double** value);
    
    void getFreeStreamVars(double* out)
    {
        for(int i = 0;i <d_NEQU;i++)
        {
            out[i] = fs_primVar[i];
        }
        out[d_NEQU] = fsnd_soundSpeed;
    }
    
    //Optional Interface for NS Equations.
    virtual void **getGradientT() 
    {
        return nullptr;
    }

    virtual double **getSpectrumViscous()
    {
        return nullptr;
    }

    virtual double getMuOverPrCell(int curMesh,int curCell)
    {
        return 0;
    }

    virtual double getMuOverPrEdge(int curMesh,int curEdge)
    {
        return 0;
    }

    virtual double getMuCell(int curMesh,int curCell)
    {
        return 0;
    }
    
    virtual double getMuEdge(int curMesh,int curEdge)
    {
        return 0;
    }

public:

    // Interface for both Euler and NS
    virtual double **getU() = 0;

    inline int** getIblankCell()
    {
        return d_iblank_cell;
    }

    // inline int** setIblankCell(int** iblank_cell_in)
    // {
    //     d_iblank_cell = iblank_cell_in;
    // }

    virtual double **getUEdge() = 0;

    virtual void** getGradientPrimitive() = 0;

    virtual double **getResidual() = 0;

    virtual double **getSpectrumConvective() = 0;

    virtual void initializeData() = 0;

    //Used by LUSGS to update Variables.
    virtual void initializeConsData(double** W);
    //Used by RungeKutta to update Variables.
    virtual void initializePrimData(double** W);

    virtual void UpdatePrimitiveVariables(double** W);

    virtual void UpdateConservativeForRecvCells(double** W);

    virtual void ConservativeParam2ConvectiveFlux(double* W, double* Fc, GeomElements::vector3d<2,double>& norm_vec);
    
    virtual void SolveViscousFlux(int curMesh,int curCell,int curEdge,int anotherCell,double* detlaW,double* Fv);

    //Used by LUSGS to perform Sweep

    virtual void SolveDiag(double** de_diag);

    virtual void SolveDF(int curMesh,int curCell,int curEdge,
                double** W_scratch ,double** deltaW0,double* DF,
                int ForB,int cellOffset);
    // virtual double getDiag(int curMesh,int curCell)
    // {
    //     return 0.0;
    // }

    // virtual double getRa(int curMesh,int curCell)
    // {
    //     return 0.0;
    // }

    virtual void preprocessAdvance(int istage) = 0;

    virtual void preprocessAdvanceSerial(int istage) = 0;

    virtual void SolveTime(double CFL) = 0;

    virtual double** getStableTime()
    {
        return d_stableDt;
    }

    virtual void postprocessAdvance(int istage) = 0;

    virtual void postprocessAdvanceSerial(int istage) = 0;

    virtual void cleaning(double **W0) = 0;

    //Especially Requested by BoundaryStrategy.
    virtual bool isInvicid() = 0;
    virtual bool isLaminar()
    {
        return true;
    }

public:

    //These are output not necessary but helps to debug code.
    void test_communication();

    void test_partialcomm();
    
    void test_unwantedSweep(int** UnwantedInd,int* nUnwantedInd);

    void writeTiogaFormat(const std::string &filename, int proc, int meshTag, int *icelldata);

        
    virtual void writeCellData(const std::string &filename, int proc, int meshTag, double **dcelldata);
    virtual void writeNOdeData(const std::string &filename, int proc, int meshTag, double **dnodedata);

    void writeRawCellData(const std::string &filename, int proc, int meshTag, double **dcelldata);

    void outPutCp(std::string& filename, int mesh_ind, int step);
};