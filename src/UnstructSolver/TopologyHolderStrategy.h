#pragma once
//This should be an abstract base class for Euler Equation or N-S Equation advancing.
//We need to specify the following routines to perform data management on TopologyHolder:
//1.Allocate data storage for Conservative Variables' Data
//2.Allocate data storage for Flux Variables' Data
//3.Estimate and allocate buffer storage for cross-process communication.
//4.Register subStrategies for Flux computation(by passing left right values or something else) and so on.
//That's all I can plan for now
#include <string>
class UnstructTopologyHolder;
class TopologyHolderStrategy
{
protected:
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
    double fs_velocity_components[2];
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
    double fs_primVar[4];

    

    // Data storage for MPI communciation
    // Each mesh has only one block as requested by TIOGA;
    // Each block has a sendBUffer for each proc;
    double **commSendBuffer = NULL;
    double **commRecvBuffer = NULL;

    // For debug.
    int **testing_flag = NULL;
friend class BoundaryStrategy;
friend class LUSGSStrategy;
friend class LaminarFluxStrategy;
    
public:
    int cur_proc = -1;
    int num_proc = -1;
public:
    TopologyHolderStrategy(UnstructTopologyHolder* hder, int nequ);
    ~TopologyHolderStrategy();

    void setNEquation(int nequ)
    {
        d_NEQU = nequ;
    }
    int getNEquation()
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

    virtual double **getMu()
    {
        return nullptr;
    }

    virtual double **getMuEdge()
    {
        return nullptr;
    }

public:

    // Interface for both Euler and NS
    virtual double **getU() = 0;

    virtual double **getUEdge() = 0;

    virtual void** getGradientPrimitive() = 0;

    virtual double **getResidual() = 0;

    virtual double **getSpectrumConvective() = 0;

    virtual void initializeData() = 0;

    virtual void preprocessAdvance(int istage) = 0;

    virtual void preprocessAdvanceSerial(int istage) = 0;

    virtual void SolveTime(double** dt, double CFL) = 0;

    virtual void postprocessAdvance() = 0;

    virtual void cleaning(double **W0) = 0;

    //Especially Requested by BoundaryStrategy.
    virtual bool isInvicid() = 0;

public:

    //These are not useful but helps to debug code.
    void test_communication();

    void test_partialcomm();
    
    void test_unwantedSweep(int** UnwantedInd,int* nUnwantedInd);

    void writeTiogaFormat(const std::string &filename, int proc, int meshTag, int *icelldata);

        
    void writeCellData(const std::string &filename, int proc, int meshTag, double **dcelldata);

    void writeRawCellData(const std::string &filename, int proc, int meshTag, double **dcelldata);


};