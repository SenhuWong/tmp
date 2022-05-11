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
public:
    int cur_proc = -1;
    int num_proc = -1;
public:
    TopologyHolderStrategy(UnstructTopologyHolder* hder);
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

    virtual void registerTopology() = 0;

public:

    // Interface for Euler and NS
    virtual double **getU() = 0;

    virtual void** getGradient() = 0;

    virtual double **getUL() = 0;

    virtual double **getUR() = 0;

    virtual double **getFluxEdge() = 0;

    virtual double **getResidual() = 0;

    virtual double **getSpectrum() = 0;
    
    virtual void initializeData() = 0;
    virtual void withinBlockCommunication() = 0;

    virtual void RemoteCellCommunication(double** value) = 0;

    virtual void NearCellCommunication(double** value) = 0;

    virtual void preprocessAdvance(int istage) = 0;

    virtual void writeCellData(const std::string &filename, int proc, int meshTag, double **dcelldata) = 0;

    virtual void preprocessAdvanceSerial(int istage) = 0;

    virtual void SolveTime(double** dt, double CFL) = 0;
    
    virtual void test_unwantedSweep(int** UnwantedInd,int* nUnwantedInd) = 0;

    virtual void postprocessAdvance() = 0;

    virtual void cleaning(double **W0) = 0;
};