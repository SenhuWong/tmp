#include "TimeStrategy.h"
#include <set>
class LUSGSStrategy : public TimeStrategy
{
private:

    //std::set<int,std::less<int>>* ForwardSweep_excluded_cells = NULL;
    int ** ForwardSweep_exceptions = NULL;
    int * nForwardSweep_exceptions = NULL;
    int ** ForwardSweep_Only = NULL;
    int * nForwardSweep_Only = NULL;

    double d_CFL_number = 400;
    double ** d_diagOperator = NULL;
    double ** d_stableDt = NULL;
    double ** W_scratch = NULL;
    double ** d_deltaW = NULL;
    double ** d_deltaW1 = NULL;
    double ** dFi = NULL;

public:
    LUSGSStrategy(UnstructTopologyHolder *hder, Euler2D* hder_strategy);
    void initialize()override;
    void singleStep(int curStage)override;
    void singleStepSerial(int curStep);
    bool converged()override;
private:
    void preprocessUpdate();
    void Update();
    void UpdateSerial();
    void postprocessUpdate();

    void preprocessLUSGS();
    void postprocessLUSGS();

    void SolveDiagOperator();

    void SolveForwardSweep();

    void SolveBackwardSweep();

    void SolveForwardSweepSerial();

    void SolveBackwardSweepSerial();

    void UpdatePrimitiveVariable();

    void UpdateConservativeForRecvCells();


};