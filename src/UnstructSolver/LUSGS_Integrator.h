#pragma once
#include "TimeIntegrator.h"
#include <set>
#include "TopologyHolderStrategy.h"
#include "toolBox/vector3d.h"
#include "toolBox/edge3d_int.h"
#include "UnstructIntegrator.h"
class LUSGSIntegrator : public TimeIntegrator
{
public:

    //std::set<int,std::less<int>>* ForwardSweep_excluded_cells = NULL;
    int ** ForwardSweep_exceptions = NULL;
    int * nForwardSweep_exceptions = NULL;
    int ** ForwardSweep_Only = NULL;
    int * nForwardSweep_Only = NULL;

    double d_CFL_number = 400;
    double ** d_diagOperator = NULL;

    double ** W_scratch = NULL;
    double ** d_deltaW = NULL;
    double ** d_deltaW1 = NULL;
    double ** dFi = NULL;

public:
    LUSGSIntegrator();
    LUSGSIntegrator(UnstructTopologyHolder *hder, TopologyHolderStrategy* hder_strategy);
    virtual void initializeData()override;
    void singleStep(int curStage)override;
    void singleStepSerial(int curStep);
    bool converged()override;
    void postprocessUpdate();
private:
    void initializeSubStrategy() override;
    void initializeParallelSweep();
    void preprocessUpdate();
    void Update();
    void UpdateSerial();
    
    //Uppering and lowering id.
    void preprocessLUSGS();
    void postprocessLUSGS();

    //Forward and Backward Sweep.
    void SolveForwardSweep();

    void SolveBackwardSweep();

    void SolveForwardSweepSerial();

    void SolveBackwardSweepSerial();

    //The skip point in here are level 1 sendcells and level 1,2 recvcells.
    void SolveLocalForwardSweep(double** diag,
                                double** W_scratch,double** deltaW1,
                                int** skip_point,int* nskip);

    //The skip point in here are level 1 cells.
    void SolveBoundaryForwardSweep(double** diag,
                                double** W_scratch,double** deltaW1,
                                int** skip_point,int * nskip);

    void SolveBoundaryBackwardSweep(double** diag,
                                    double** W_scratch,double** deltaW1,double** deltaW,
                                    int ** skip_point,int * nskip);

    void SolveLocalBackwardSweep(double** diag,
                                double** W_scratch,double** deltaW1,double** deltaW,
                                int** skip_point,int* nskip);

    void SolveForwardSweep(double** diag,
                        double** W_scratch,double** deltaW1);

    void SolveBackwardSweep(double** diag,
                            double** W_scratch,double** deltaW1,double** deltaW);

    //PostprocessUpdate
    void UpdatePrimitiveVariable();

    
};