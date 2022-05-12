#include "TimeStrategy.h"
#include <set>
#include "TopologyHolderStrategy.h"
#include "toolBox/vector3d.h"
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
    LUSGSStrategy(UnstructTopologyHolder *hder, TopologyHolderStrategy* hder_strategy);
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

    void ConservativeParam2ConvectiveFlux(double* W, double* Fc, GeomElements::vector3d<2,double>& norm_vec)
    {
        double Gamma = 1.4;
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

    void SolveDiag(double** de_diag, double** dt);

    
    void SolveDF(int curMesh,int curCell,int curEdge,
                double** W_scratch ,double** deltaW0,double* DF,
                int ForB,int cellOffset);

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

};