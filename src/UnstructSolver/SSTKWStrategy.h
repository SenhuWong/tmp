#include "TopologyHolderStrategy.h"
#include "UnstructIntegrator.h"
#include "toolBox/edge3d_int.h"
#include "toolBox/vector3d.h"
#include "LaminarFluxStrategy.h"
#include <math.h>

template<int ndim>
class SSTkomegaModel:public TopologyHolderStrategy
{
private:
    //Constants used by SSTModel.
    const double C1 = 10.0;
    const double C2 = 3.0;
    const double L = 1.0;

    const double beta1 = 0.075;
    const double beta2 = 0.0828;
    const double sigma_k1 = 0.85;
    const double sigma_k2 = 1.0;
    const double sigma_w1 = 0.5;
    const double sigma_w2 = 0.856;
    const double Cw1 = 0.533;
    const double Cw2 = 0.440;

    const double a1 = 0.31;
    const double _beta = 0.09;
    const double Cepsilion = 0.61;
    const double Comega = 0.78;
    //Arrays used by SSTModel.
    double** cellF1 = NULL;
    double** cellF2 = NULL;
    double** vCellMut = NULL;
    double** vFaceMut = NULL;
    //FreeStream Variables.
    double fs_mut = 0.0;
    double fsnd_mut = 0.0;
    double fsnd_mu = 0.0;

    double Gamma = 1.4;
    double fs_turbulenceValue[2];
    double fsnd_turbulenceValue[2];

    TopologyHolderStrategy* d_flow_strategy;
    LaminarFluxStrategy* d_laminar_strategy;

    double** turbulenceCellValue = NULL;//These are  now actually U
    GeomElements::vector3d<ndim,double>** turbulenceCellGrad = NULL;//These are now actually GradU.



    double** turbulenceFaceValue = NULL;
    double** FaceConvectiveTFlux = NULL;
    double** FaceViscousTFlux = NULL;
    double** CellTotalTFlux = NULL;//These are now actually Residual

    
    int nequ = 0;// Flow Strategy's d_NEQU;
    
    //Get from TIOGA;
    double** vCellSolidDistance = NULL;

    int getNEquation() override
    {
        return 2;
    }


public:
    double** extern_W_scratch = NULL;
    double** extern_deltaW = NULL;
    double** extern_deltaW1 = NULL;
    //Interface for TopologyHolderStrategy
    double** getU() override
    {
        return turbulenceCellValue;

    }
    double** getUEdge() override
    {
        return turbulenceFaceValue;

    }
    void** getGradientPrimitive() override
    {
        return (void**)(turbulenceCellGrad);
    }

    double** getResidual() override
    {
        return CellTotalTFlux;

    }

    double** getSpectrumConvective() override
    {
        return nullptr;
    }

    void initializeConsData(double** W) override
    {
        double** U = d_flow_strategy->getU();
        for(int i = 0;i<d_nmesh;i++)
        {
            for(int k = 0;k<d_hder->nCells(i);k++)
            {
                double rho = U[i][nequ*k];
                W[i][2*k] = turbulenceCellValue[i][2*k]*rho;
                W[i][2*k+1] = turbulenceCellValue[i][2*k+1]*rho;                
            }
        }
    }

    void initializePrimData(double** W) override
    {
        double **U = d_flow_strategy->getU();
        for(int i = 0;i<d_nmesh;i++)
        {
            for(int k = 0;k<d_hder->nCells(i);k++)
            {
                double rho = U[i][nequ*k];
                turbulenceCellValue[i][2*k] = W[i][2*k] / rho;
                turbulenceCellValue[i][2*k+1] = W[i][2*k+1]/rho;     
            }
        }
    }

    void initializeData() override
    {
        InitializeturbulenceVariable();
    }

    void setStableTime()
    {
        d_stableDt = d_flow_strategy->getStableTime();
    }

    void preprocessAdvance(int istage)
    {
        //Compute of Residual should have been done here.
    }

    void preprocessAdvanceSerial(int istage)
    {

    }

    void UpdateConservativeForRecvCells(double** W)
    {
        double Gamma = 1.4;
        double** U = d_flow_strategy->getU();
        for(int i = 0;i<d_nmesh;i++)
        {
            for(auto iter = d_hder->recvBegin(i);iter!=d_hder->recvEnd(i);iter++)
            {
                int cellId = iter->d_localId;
                double rho = U[i][nequ*cellId];
                W[i][2*cellId] = rho*turbulenceCellValue[i][2*cellId];
                W[i][2*cellId+1] = rho*turbulenceCellValue[i][2*cellId+1];
            }
        }
    }

    void SolveTime(double CFL)
    {
        //Do nothing
        return;
    }
    //These 3 for LUSGS must be overridden.
    void SolveDiag(double** de_diag) override
    {
        double** spectrum_convec = d_flow_strategy->getSpectrumConvective();
        double** spectrum_viscous = d_flow_strategy->getSpectrumViscous();
        double OMEGAN =1.5;
        for(int i = 0;i<d_nmesh;i++)
        {
            auto& curBlk = d_hder->blk2D[i];
            for(int k = 0;k < d_hder->nCells(i);k++)
            {

                auto& curCell = curBlk.d_localCells[k];
                double diag = curCell.volume()/d_stableDt[i][k] + 0.5*OMEGAN*spectrum_convec[i][k] + spectrum_viscous[i][k];
                de_diag[i][2*k] =  diag+2.0* _beta * turbulenceCellValue[i][2*k+1]*curCell.volume();

                double gradVal = turbulenceCellGrad[i][2*k].dot_product(turbulenceCellGrad[i][2*k+1]);
                double beta = cellF1[i][k]*beta1 + (1-cellF1[i][k])*beta2;
                double w = turbulenceCellValue[i][2*k+1];
                de_diag[i][2*k+1] = diag+curCell.volume()*fabsf64(2.0*beta*w+2.0*sigma_w2*(1-cellF1[i][k])*gradVal/(w*w));
            }
        }

    }

    void UpdatePrimitiveVariables(double** W) override
    {
        //This has to be done only after FLowStrategy's UpdatePrimitiveVariables
        double** U = d_flow_strategy->getU();
        GeomElements::vector3d<ndim,double>** gradU = (GeomElements::vector3d<ndim,double>**)d_flow_strategy->getGradientPrimitive();
        double rho;
        double curl;
        for(int i = 0;i<d_nmesh;i++)
        {
            for(int k = 0;k<d_hder->nCells(i);k++)
            {
                rho = U[i][nequ*k];
                turbulenceCellValue[i][2*k] = W[i][2*k]/rho;
                turbulenceCellValue[i][2*k+1] = W[i][2*k+1]/rho;
                // if(turbulenceCellValue[i][2*k]<fs_turbulenceValue[0]*1.0e-4)
                // {
                //     turbulenceCellValue[i][2*k] = fs_turbulenceValue[0]*1.0e-4;
                //     W[i][2*k] = rho * turbulenceCellValue[i][2*k];
                // }
                // if(turbulenceCellValue[i][2*k+1]<fs_turbulenceValue[1]*1.0e-4)
                // {
                //     turbulenceCellValue[i][2*k+1] = fs_turbulenceValue[1]*1.0e-4;
                //     W[i][2*k+1] = rho * turbulenceCellValue[i][2*k+1];
                // }
                curl = fabsf64(gradU[i][3+nequ*k][0] - gradU[i][2+nequ*k][1]);
                double h=std::max(a1*W[i][2*k+1]/rho,cellF2[i][k]*curl);
                vCellMut[i][k]=a1*W[i][2*k]/h;
            }
        }
    }

    //THese are interfaces need later rejust;
    virtual void postprocessAdvance(int iStage)
    {
        std::cout<<"Do nothing and return\n";
        return;

    }
    virtual void postprocessAdvanceSerial(int iStage)
    {
        std::cout<<"Do nothing and return\n";
        return;
    }
    virtual void cleaning(double **W0)
    {
        return;

    }

    //Especially Requested by BoundaryStrategy.
    virtual bool isInvicid()
    {
        return false;
    }

    //End of Interface

    void SolveDF(int curMesh,int curCell,int curEdge,
            double** W_scratch ,double** deltaW0,double* DF,
            int ForB,int cellOffset) override
    {
        double** W_scratch_extern = extern_W_scratch;
        double** deltaW0_extern = ForB == 1 ? extern_deltaW1 : extern_deltaW;
        const double OMEGAN  =1.5;
        auto& edg = d_hder->blk2D[curMesh].d_localEdges[curEdge];
        double Wp[nequ+d_NEQU];
        double W[nequ+d_NEQU];
        double dW[nequ+d_NEQU];
        double W0[nequ+d_NEQU];
        double Fcp[d_NEQU];
        double Fc[d_NEQU];

        double rA;
        double Fv[d_NEQU];
        double DFc[d_NEQU];

        double** U = d_flow_strategy->getU();
        int lC = edg.lCInd();
        int rC = edg.rCInd();
        if(curCell + cellOffset==lC)//Current cell at left side
        {
            bool dothis = false;
            if(rC>=0 and ((int(rC<lC))==(ForB)))
            {
                //If ForB is 1, then calculate when curCell(lC) >the other cell(rC)
                //If ForB is 0, then calculate when curCell(lC) < the other cell(rC)
                if(rC >= d_hder->nCells(curMesh))
                {
                    rC = rC - d_hder->nCells(curMesh);
                }
                for(int j = 0;j<nequ;j++)
                {
                    W0[j] = W_scratch_extern[curMesh][j+nequ*rC];
                    dW[j] = deltaW0_extern[curMesh][j+nequ*rC];
                    Wp[j] =  W0[j] + dW[j];
                    W[j] = W0[j];
                }
                for(int j = nequ;j<nequ+d_NEQU;j++)
                {
                    W0[j] = W_scratch[curMesh][j+d_NEQU*rC];
                    dW[j] = deltaW0[curMesh][j+d_NEQU*rC];
                    Wp[j] = W0[j] + dW[j];
                    W[j] = W0[j];
                }
                


                //Get convectiveFlux
                // GeomElements::vector3d<ndim,double> curVelocity(&U[curMesh][2+nequ*rC]);
                GeomElements::vector3d<2,double> norm = edg.normal_vector();
                ConservativeParam2ConvectiveFlux(W,Fc,norm);
                ConservativeParam2ConvectiveFlux(Wp,Fcp,norm);
                SolveViscousFlux(curMesh,curCell,curEdge,rC,dW,Fv);
                
            //     GeomElements::vector3d<ndim,double> norm = edg.normal_vector();
            //     double Vn = curVelocity.dot_product(norm);
            //     double c = sqrtf64(Gamma*U[curMesh][1+nequ*rC]/U[curMesh][nequ*rC]);
            //     GeomElements::vector3d<ndim,double> centerLine;
         
            //     centerLine = d_hder->blk2D[curMesh].d_localCells[rC].center() - d_hder->blk2D[curMesh].d_localCells[curCell].center();

            //     double rA = OMEGAN*(fabsf64(Vn)+c) + std::max<double>(4/(3*U[curMesh][0+nequ*rC]),Gamma/U[curMesh][nequ*rC])*
            //                 (d_flow_strategy->getMuOverPrCell(curMesh,curCell))/sqrtf64(centerLine.L2Square());
            //      //std::cout<<"Vn and rA"<<Vn<<" "<<rA<<'\n';
            //      double efficient;
                 
            //      //if(rC == GeomElements::edge3d<2>::BoundaryType::WALL)
            //    // {efficient = (((-Vn) - rA))*0.5*edg.area();}
            //     //else{
            //         efficient = (((Vn) - rA))*0.5*edg.area();

                
                // for(int j = 0;j<d_NEQU;j++)
                // {
                //     DF[j] = deltaW0[curMesh][j+d_NEQU*rC]*efficient;
                // }
                for(int j = 0;j<d_NEQU;j++)
                {
                    DFc[j] = Fcp[j] - Fc[j];

                }
                for(int j = 0;j<d_NEQU;j++)
                {
                    DF[j] = 0.5*(DFc[j] - Fv[j])*edg.area();
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
            if((int(lC<rC)) == (ForB))
            {
                if(lC>=d_hder->nCells(curMesh))
                {
                    lC = lC - d_hder->nCells(curMesh);
                }

                for(int j = 0;j<nequ;j++)
                {
                    W0[j] = W_scratch_extern[curMesh][j+nequ*lC];
                    dW[j] = deltaW0_extern[curMesh][j+nequ*lC];
                    Wp[j] =  W0[j] + dW[j];
                    W[j] = W0[j];
                }
                for(int j = nequ;j<nequ+d_NEQU;j++)
                {
                    W0[j] = W_scratch[curMesh][j+d_NEQU*lC];
                    dW[j] = deltaW0[curMesh][j+d_NEQU*lC];
                    Wp[j] = W0[j] + dW[j];
                    W[j] = W0[j];
                }
                //Get convectiveFlux
                //Get convectiveFlux
                // GeomElements::vector3d<ndim,double> curVelocity(&U[curMesh][2+nequ*rC]);
                GeomElements::vector3d<2,double> norm = edg.normal_vector();
                ConservativeParam2ConvectiveFlux(W,Fc,norm);
                ConservativeParam2ConvectiveFlux(Wp,Fcp,norm);
                SolveViscousFlux(curMesh,curCell,curEdge,rC,dW,Fv);

                // GeomElements::vector3d<ndim,double> curVelocity(&U[curMesh][2+nequ*lC]);
                // GeomElements::vector3d<ndim,double> norm = edg.normal_vector();
                // double Vn = curVelocity.dot_product(norm);
                // double c = sqrtf64(Gamma*U[curMesh][1+nequ*lC]/U[curMesh][nequ*lC]);
                // GeomElements::vector3d<ndim,double> centerLine;
                
                // centerLine = d_hder->blk2D[curMesh].d_localCells[lC].center() - d_hder->blk2D[curMesh].d_localCells[curCell].center();

                // double rA = OMEGAN*(fabsf64(Vn)+c) + std::max<double>(4/(3*U[curMesh][0+nequ*lC]),Gamma/U[curMesh][nequ*lC])*
                //             (d_flow_strategy->getMuOverPrCell(curMesh,curCell))/sqrtf64(centerLine.L2Square());
                // double efficient= (-Vn - rA)*0.5*edg.area();
                // for(int j = 0;j<d_NEQU;j++)
                // {
                //     DF[j] = deltaW0[curMesh][j+d_NEQU*lC]*efficient;
                // }
                for(int j = 0;j<d_NEQU;j++)
                {
                    DFc[j] = Fc[j] - Fcp[j];

                }
                for(int j = 0;j<d_NEQU;j++)
                {
                    DF[j] = 0.5*(DFc[j] - Fv[j])*edg.area();
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
            std::cout<<curCell<<" with offset "<< cellOffset <<" not found in "<<lC<<" and "<<rC<<'\n';
            throw std::runtime_error("Cur cell is not on either side,impossible\n");
        }
    }
    void ConservativeParam2ConvectiveFlux(double* W, double* Fc, GeomElements::vector3d<2,double>& norm_vec) override
    {
        double Gamma = 1.4;
        GeomElements::vector3d<2,double> velocity(W[2]/W[0],W[3]/W[0]);
        double Vn = velocity.dot_product(norm_vec);
        double DatCell = W[0];
        double PatCell = (Gamma-1)*(W[1] - 0.5*DatCell*velocity.L2Square());
        double rhoEatCell = W[1];
        Fc[0] = W[nequ+0]*Vn;
        Fc[1] = W[nequ+1]*Vn;
    }
    void SolveViscousFlux(int curMesh,int curCell,int curEdge,int anotherCell,double* detlaW,double* Fv) override
    {
        const double Gamma = 1.4;
        const double OMEGAN = 1.5;
        auto& curBlk = d_hder->blk2D[curMesh];
        double rA;
        double** U_extern =d_flow_strategy->getU();
        auto& curE = curBlk.d_localEdges[curEdge];
        auto& curC = curBlk.d_localCells[curCell];
        double edgeAreaSquare = curE.area()*curE.area();
        GeomElements::vector3d<2,double> centerLine;
        //It is impossible to have anotherCell at first negative to be raised amid 0 and nCells 
        if(anotherCell>= d_hder->nCells(curMesh))
        {
            throw std::runtime_error("Oh Look what you have done\n");
        }
        // std::cout<<"Viscous Flux computed\n";
        if(anotherCell >=0)
        {
            auto& anotherC = curBlk.d_localCells[anotherCell];
            centerLine = anotherC.center() - curC.center();
        }
        else
        {
            centerLine = curE.center() - curC.center();
        }
        if(anotherCell<0)
        {
            throw std::runtime_error("Another Cell is below 0\n");
        }
        GeomElements::vector3d<2,double> curVelocity(&U_extern[curMesh][2+nequ*anotherCell]);
        GeomElements::vector3d<2,double> norm = curE.normal_vector();
        double Vn = curVelocity.dot_product(norm);
        double c = sqrtf64(Gamma * U_extern[curMesh][1+nequ*anotherCell] / U_extern[curMesh][0+nequ*anotherCell]);

        rA  = OMEGAN*(fabsf64(Vn)+c);
        if(!isInvicid())
        {
            rA += std::max<double>(4/(3*U_extern[curMesh][0+nequ*anotherCell]),Gamma/U_extern[curMesh][0+nequ*anotherCell])*
            (d_flow_strategy->getMuOverPrCell(curMesh,curCell))/sqrtf64(centerLine.L2Square());
        }
        for(int j = 0;j<d_NEQU;j++)
        {
            Fv[j] = rA*detlaW[nequ+j];
        }
    }

    void ComputeFaceTValue();
    void ComputeCellTGradient();

    void ComputeSSTf1f2();

    void ComputeTurbulenceViscousCoefficient();

    void SolveConvectiveTFlux();

    void SolveViscousTFlux();

    double rans_P_Solver(int meshIndex,int cellIndex);

    void SolveTotalTFlux();

    SSTkomegaModel(UnstructTopologyHolder* hder,TopologyHolderStrategy* flower,LaminarFluxStrategy* laminar);

    double getMuOverPrCell(int curMesh,int curCell)
    {
        return vCellMut[curMesh][curCell]/fs_Prt;
    }

    double getMuOverPrEdge(int curMesh,int curEdge)
    {
        return vFaceMut[curMesh][curEdge]/fs_Prt;
    }

    virtual double getMuCell(int curMesh,int curCell)
    {
        return vCellMut[curMesh][curCell];
    }
    
    virtual double getMuEdge(int curMesh,int curEdge)
    {
        return vFaceMut[curMesh][curEdge];
    }

    void InitializeturbulenceVariable();

    void SolveTurbulenceEquation()
    {
        ComputeFaceTValue();
        ComputeCellTGradient();
        ComputeSSTf1f2();
        ComputeTurbulenceViscousCoefficient();
        SolveConvectiveTFlux();
        SolveViscousTFlux();
        SolveTotalTFlux();
    }

    void LUSGSAdvance();

    //The last 3 looks like an AMR routine to me,I will leave it here.
    void writeCellData(const std::string &filename, int proc, int meshTag, double **dcelldata);
};


template<int ndim>
SSTkomegaModel<ndim>::SSTkomegaModel(UnstructTopologyHolder* hder,TopologyHolderStrategy* flower,LaminarFluxStrategy* laminar)
:TopologyHolderStrategy(hder,2)
{
    turbulenceCellValue = new double*[d_nmesh];
    //turbulenceConsCellValue = new double*[d_nmesh];
    cellF1 = new double*[d_nmesh];
    cellF2 = new double*[d_nmesh];
    turbulenceCellGrad = new GeomElements::vector3d<ndim,double>*[d_nmesh];

    vCellMut = new double*[d_nmesh];
    vFaceMut = new double*[d_nmesh];

    turbulenceFaceValue = new double*[d_nmesh];
    FaceConvectiveTFlux = new double*[d_nmesh];
    FaceViscousTFlux = new double*[d_nmesh];
    CellTotalTFlux = new double*[d_nmesh];

    for(int i = 0;i<d_nmesh;i++)
    {
        turbulenceCellValue[i] = new double[d_NEQU*d_hder->nCells(i)];

        cellF1[i] = new double[d_hder->nCells(i)];
        cellF2[i] = new double[d_hder->nCells(i)];
        turbulenceCellGrad[i] = new GeomElements::vector3d<ndim,double>[d_NEQU*d_hder->nCells(i)];

        vCellMut[i] = new double[d_hder->nCells(i)];
        vFaceMut[i] = new double[d_hder->nEdges(i)];
        turbulenceFaceValue[i] = new double[d_NEQU*d_hder->nEdges(i)];

        FaceConvectiveTFlux[i] = new double[d_NEQU*d_hder->nCells(i)];
        FaceViscousTFlux[i] = new double[d_NEQU*d_hder->nCells(i)];
        CellTotalTFlux[i] = new double[d_NEQU*d_hder->nCells(i)];
    }
    d_flow_strategy = flower;
    d_laminar_strategy = laminar;
    std::cout<<"End of SSTKomega Model\n";
    std::cout<<"Before :"<<fs_density<<'\t';
    TopologyHolderStrategy::set_fs_variable(d_flow_strategy);

    std::cout<<"After :"<<TopologyHolderStrategy::fs_density<<'\n';
    std::cout<<"Holder's density: "<<d_flow_strategy->fs_density<<'\n';

    
    nequ = d_flow_strategy->getNEquation();
    setStableTime();
}

template<int ndim>
void SSTkomegaModel<ndim>::InitializeturbulenceVariable()
{
    std::cout<<"p and d "<<fs_pressure<<" "<<fs_density<<'\n';
    fs_mut = fs_mu*0.035;
    fsnd_mu = sqrtf64(Gamma)*fs_mach/fs_Re;
    fsnd_mut = fs_mut/fs_mu*sqrtf64(Gamma)*fs_mach/fs_Re;
    std::cout<<"fsnd_mu is "<<fsnd_mu<<'\n';
    std::cout<<"fsnd_mut is "<<fsnd_mut<<'\n';
    
    double fsnd_velocity_magnitude = GeomElements::vector3d<ndim,double>(&fs_primVar[2]).normalize();


    std::cout<<"fs_velocity_magnitude is "<<fsnd_velocity_magnitude<<'\n';

    fs_turbulenceValue[0] = 1.5*fs_mach*fs_soundSpeed*0.01*fs_mach*fs_soundSpeed*0.01;
    fs_turbulenceValue[1] = fs_density * fs_turbulenceValue[0] /(0.035*fs_mu);

    fsnd_turbulenceValue[0]=Gamma*fs_turbulenceValue[0]/(fs_soundSpeed*fs_soundSpeed);
    fsnd_turbulenceValue[1]=sqrt(Gamma)*L*fs_turbulenceValue[1]/fs_soundSpeed;

    for(int i = 0;i<d_nmesh;i++)
    {
        for(int k = 0;k<d_hder->nCells(i);k++)
        {
            int offset= 2*k;
            turbulenceCellValue[i][offset] = fsnd_turbulenceValue[0];
            turbulenceCellValue[i][offset+1] = fsnd_turbulenceValue[1];
            vCellMut[i][k] = fsnd_mut;
        }
    }
    std::cout<<"turbulenceValue: "<<fsnd_turbulenceValue[0]<<" "<<fsnd_turbulenceValue[1]<<'\n';
}

template<int ndim>
void SSTkomegaModel<ndim>::ComputeFaceTValue()
{
    double** vFaceMul = d_laminar_strategy->getMuEdge();
    double** U_edge = d_flow_strategy->getUEdge();
    for(int i = 0;i<d_nmesh;i++)
    {
        auto& curBlk = d_hder->blk2D[i];
        double* distances = curBlk.wallCellDistance;
        for(int k =0;k<d_hder->nEdges(i);k++)
        {
            auto& curEdge = curBlk.d_localEdges[k];
            int lC = curEdge.lCInd();
            int rC = curEdge.rCInd();
            if(rC>=0)//Inner
            {
                turbulenceFaceValue[i][2*k] = 0.5*(turbulenceCellValue[i][2*lC] + turbulenceCellValue[i][2*rC]);
                turbulenceFaceValue[i][2*k+1] = 0.5*(turbulenceCellValue[i][2*lC+1] + turbulenceCellValue[i][2*rC+1]);
            }
            else if(rC == GeomElements::edge3d<2>::BoundaryType::WALL)
            {
                turbulenceFaceValue[i][2*k] = 0.0;
            
                double dis = distances[lC];
                turbulenceFaceValue[i][2*k+1] = 10.0 * 6.0 * vFaceMul[i][k] / (U_edge[i][nequ*k]*beta1*dis*dis);
            }
            else if(rC == GeomElements::edge3d<2>::BoundaryType::FARFIELD)
            {
                GeomElements::vector3d<ndim,double> VelocityEdge(&U_edge[i][nequ*k+2]);
                GeomElements::vector3d<ndim,double> norm = curEdge.normal_vector();
                double Vn = VelocityEdge.dot_product(norm);
                if(Vn > 0.0)
                {
                    turbulenceFaceValue[i][2*k] = turbulenceCellValue[i][2*lC];
                    turbulenceFaceValue[i][2*k+1] = turbulenceCellValue[i][2*lC+1];
                }
                else
                {
                    turbulenceFaceValue[i][2*k] = fsnd_turbulenceValue[0];
                    turbulenceFaceValue[i][2*k+1] = fsnd_turbulenceValue[1];
                }
            }
        }
    }
}

template<int ndim>
void SSTkomegaModel<ndim>::ComputeCellTGradient()
{
    for(int i = 0;i<d_nmesh;i++)
    {
        auto& curBlk = d_hder->blk2D[i];
        for(int k = 0;k<d_hder->nCells(i);k++)
        {
            turbulenceCellGrad[i][2*k].reset();
            turbulenceCellGrad[i][2*k+1].reset();
        }
        for(int k = 0;k<d_hder->nEdges(i);k++)
        {
            auto& curEdge = curBlk.d_localEdges[k];
            int lC = curEdge.lCInd();
            int rC = curEdge.rCInd();
            turbulenceCellGrad[i][2*lC] = turbulenceCellGrad[i][2*lC] + curEdge.normal_vector()*(curEdge.area()*turbulenceFaceValue[i][2*k]);
            turbulenceCellGrad[i][2*lC+1] = turbulenceCellGrad[i][2*lC+1] + curEdge.normal_vector()*(curEdge.area()*turbulenceFaceValue[i][2*k+1]);
            if(rC>=0)
            {
                turbulenceCellGrad[i][2*rC] = turbulenceCellGrad[i][2*rC] - curEdge.normal_vector()*(curEdge.area()*turbulenceFaceValue[i][2*k]);
                turbulenceCellGrad[i][2*rC+1] = turbulenceCellGrad[i][2*rC+1] - curEdge.normal_vector()*(curEdge.area()*turbulenceFaceValue[i][2*k+1]);
            }
        }
        for(int k = 0;k<d_hder->nCells(i);k++)
        {
            auto& curCell = curBlk.d_localCells[k];
            turbulenceCellGrad[i][2*k] = turbulenceCellGrad[i][2*k] / curCell.volume();
            turbulenceCellGrad[i][2*k+1] = turbulenceCellGrad[i][2*k+1] / curCell.volume(); 
        }
    }
}

template<int ndim>
void SSTkomegaModel<ndim>::ComputeSSTf1f2()
{
    double** U = d_flow_strategy->getU();
    double** mu = d_laminar_strategy->getMu();
    double CD_kw, arg1, arg2, f1,f2;
    for(int i = 0;i<d_nmesh;i++)
    {
        auto& curBlk = d_hder->blk2D[i];
        double* distances = curBlk.wallCellDistance;
        for(int k = 0;k<d_hder->nCells(i);k++)
        {
            double gradVal = turbulenceCellGrad[i][2*k].dot_product(turbulenceCellGrad[i][2*k+1]);
            double rho = U[i][nequ*k];
            double dis = distances[k];
            double K = turbulenceCellValue[i][2*k];
            double w = turbulenceCellValue[i][2*k+1];
            double muL = mu[i][k];
            CD_kw = std::max<double>(2.0*rho*sigma_w2 / w * gradVal,1.0e-20);
            arg1 = std::max<double>(sqrtf64(K)/(0.09*w*dis),500.0*muL/(rho*w*dis*dis));
            arg1 = std::min<double>(arg1,4.0*rho*sigma_w2*K/(CD_kw*dis*dis));
            arg2 = std::max(2.0*sqrtf64(K)/(0.09*w*dis),500.0*muL/(rho*w*dis*dis));
            f1 = tanhf64(arg1*arg1*arg1*arg1);
            f2 = tanhf64(arg2*arg2);
            cellF1[i][k] = f1;
            cellF2[i][k] = f2;
        }
    }
}

template<int ndim>
void SSTkomegaModel<ndim>::ComputeTurbulenceViscousCoefficient()
{
    //First compute curlV
    double curl;
    GeomElements::vector3d<ndim,double>** gradU = (GeomElements::vector3d<ndim,double>**)d_flow_strategy->getGradientPrimitive();
    double** U = d_flow_strategy->getU();
    double** U_edge = d_flow_strategy->getUEdge();
    if(ndim==2)
    {
        // for(int i = 0;i<d_nmesh;i++)
        // {
        //     for(int k = 0;k<d_hder->nCells(i);k++)
        //     {
        //         curl = fabsf64(gradU[i][3+nequ*k][0] - gradU[i][2+nequ*k][1]);
        //         double rho = U[i][nequ*k];
        //         vCellMut[i][k]= a1 *rho* turbulenceCellValue[i][2*k] / std::max<double>(a1 * turbulenceCellValue[i][2*k+1], cellF2[i][k]* curl);
        //     }
        // }
    }
    else if(ndim==3)
    {
        throw std::runtime_error("Not implemented yet\n");
    }
                
    for(int i = 0;i<d_nmesh;i++)
    {
        auto& curBlk = d_hder->blk2D[i];
        for(int k = 0;k<d_hder->nEdges(i);k++)
        {
            auto& curEdge =curBlk.d_localEdges[k];
            int lC = curEdge.lCInd();
            int rC = curEdge.rCInd();
            if(rC>=0)
            {
                vFaceMut[i][k] = 0.5*(vCellMut[i][lC] + vCellMut[i][rC]);
            }
            else if(rC==GeomElements::edge3d<2>::BoundaryType::WALL)
            {
                vFaceMut[i][k] = 0.0;
            }
            else if(rC==GeomElements::edge3d<2>::BoundaryType::FARFIELD)
            {
                GeomElements::vector3d<ndim,double> velo(&(U_edge[i][2+nequ*k]));
                GeomElements::vector3d<ndim,double> norm = curEdge.normal_vector();
                double Vn = velo.dot_product(norm);
                vFaceMut[i][k]  = Vn > 0.0 ? vCellMut[i][lC] : fsnd_mut;
            }
            else
            {
                throw std::runtime_error("Undefined boundary type\n");
            }

        }
    }
    
}

template<int ndim>
void SSTkomegaModel<ndim>::SolveConvectiveTFlux()
{
    double** U =d_flow_strategy->getU();
    double** U_edge = d_flow_strategy->getUEdge();
    for(int i = 0;i<d_nmesh;i++)
    {
        auto& curBlk = d_hder->blk2D[i];
        for(int k = 0;k<d_hder->nEdges(i);k++)
        {
            auto& curEdge = curBlk.d_localEdges[k];
            int lC = curEdge.lCInd();
            int rC = curEdge.rCInd();
            GeomElements::vector3d<ndim,double> velo(&(U_edge[i][2+nequ*k]));
            GeomElements::vector3d<ndim,double> norm = curEdge.normal_vector();
            double Vn = velo.dot_product(norm);
            double VolumePassing = (Vn * curEdge.area());
            if(rC>=0)
            {
                int cellIndex = Vn > 0.0 ? lC:rC;
                double Vp=0.5*(Vn+fabsf64(Vn));
                double Vs=0.5*(Vn-fabsf64(Vn));

                FaceConvectiveTFlux[i][2*k] = turbulenceCellValue[i][2*lC]* U[i][nequ*lC]*Vp*curEdge.area()
                 +turbulenceCellValue[i][2*rC]* U[i][nequ*rC]*Vs*curEdge.area();
                FaceConvectiveTFlux[i][2*k+1] = turbulenceCellValue[i][2*lC+1]* U[i][nequ*lC]*Vp*curEdge.area()
                 +turbulenceCellValue[i][2*rC+1]* U[i][nequ*rC]*Vs*curEdge.area();
            }
            else 
            {
                FaceConvectiveTFlux[i][2*k] = Vn*U_edge[i][nequ*k]*turbulenceFaceValue[i][2*k]*curEdge.area() ;
                FaceConvectiveTFlux[i][2*k+1] = Vn*U_edge[i][nequ*k]*turbulenceFaceValue[i][2*k+1]*curEdge.area();
            }
            // if(rC>=0)
            // {
            //     int cellIndex = Vn > 0.0 ? lC:rC;
            //     double massPassing = VolumePassing * U[i][nequ*cellIndex];
            //     FaceConvectiveTFlux[i][2*k] = turbulenceCellValue[i][2*cellIndex] * massPassing;
            //     FaceConvectiveTFlux[i][2*k+1] = turbulenceCellValue[i][2*cellIndex+1] * massPassing; 
            // }
            // else
            // {
            //     double massPassing = VolumePassing * U_edge[i][nequ*k];
            //     FaceConvectiveTFlux[i][2*k] = turbulenceFaceValue[i][2*k] * massPassing;
            //     FaceConvectiveTFlux[i][2*k+1] = turbulenceFaceValue[i][2*k+1] * massPassing;
            // }
        }
    }

}

template<int ndim>
void SSTkomegaModel<ndim>::SolveViscousTFlux()
{
    double** vFaceMul = d_laminar_strategy->getMuEdge();
    double txx_k, tyy_k, txx_w, tyy_w;
    double tzz_k, tzz_w;
    double left_f1, right_f1, face_f1, sigma_k,sigma_w, faceMul, faceMut;
    double left_Tvalue[2], right_Tvalue[2];
    GeomElements::vector3d<ndim,double> faceTGrad[2];
    GeomElements::vector3d<ndim,double> left2Right;
    for(int i = 0;i<d_nmesh;i++)
    {
        auto& curBlk = d_hder->blk2D[i];
        for(int k = 0;k<d_hder->nEdges(i);k++)
        {
            auto& curEdge = curBlk.d_localEdges[k];
            auto norm = curEdge.normal_vector();
            int lC = curEdge.lCInd();
            int rC = curEdge.rCInd();
            if(rC>=0)
            {
                left2Right = curBlk.d_localCells[rC].center() - curBlk.d_localCells[lC].center();
                
                left_f1 = cellF1[i][lC];
                left_Tvalue[0] = turbulenceCellValue[i][2*lC];
                left_Tvalue[1] = turbulenceCellValue[i][2*lC+1];

                right_f1 = cellF1[i][rC];
                right_Tvalue[0] = turbulenceCellValue[i][2*rC];
                right_Tvalue[1] = turbulenceCellValue[i][2*rC+1];

                face_f1 = 0.5*(left_f1+right_f1);
                faceTGrad[0] = (turbulenceCellGrad[i][2*lC] + turbulenceCellGrad[i][2*rC])*0.5;
                faceTGrad[1] = (turbulenceCellGrad[i][2*lC+1] + turbulenceCellGrad[i][2*rC+1])*0.5;

                faceTGrad[0][0]=faceTGrad[0][0]+ ((right_Tvalue[0] - left_Tvalue[0])-(faceTGrad[0][0]*left2Right[0]+faceTGrad[0][1]*left2Right[1]))*left2Right[0]/left2Right.L2Square();
                faceTGrad[0][1]=faceTGrad[0][1]+ ((right_Tvalue[0] - left_Tvalue[0])-(faceTGrad[0][0]*left2Right[0]+faceTGrad[0][1]*left2Right[1]))*left2Right[1]/left2Right.L2Square();
           
                faceTGrad[1][0]=faceTGrad[1][0]+ ((right_Tvalue[1] - left_Tvalue[1])-(faceTGrad[1][0]*left2Right[0]+faceTGrad[1][1]*left2Right[1]))*left2Right[0]/left2Right.L2Square();
                faceTGrad[1][1]=faceTGrad[1][1]+ ((right_Tvalue[1] - left_Tvalue[1])-(faceTGrad[1][0]*left2Right[0]+faceTGrad[1][1]*left2Right[1]))*left2Right[1]/left2Right.L2Square();
            }
            else
            {
                left2Right = curEdge.center() - curBlk.d_localCells[lC].center();
                
                
                left_Tvalue[0] = turbulenceCellValue[i][2*lC];
                left_Tvalue[1] = turbulenceCellValue[i][2*lC+1];
                
                right_Tvalue[0] = turbulenceFaceValue[i][2*k];
                right_Tvalue[1] = turbulenceFaceValue[i][2*k+1];

                face_f1 = cellF1[i][lC];
                faceTGrad[0] = turbulenceCellGrad[i][2*lC];
                faceTGrad[1] = turbulenceCellGrad[i][2*lC+1];

            }
            // faceTGrad[0] = faceTGrad[0] + left2Right*((right_Tvalue[0] - left_Tvalue[0]-faceTGrad[0].dot_product(left2Right))/left2Right.L2Square());
            // faceTGrad[1] = faceTGrad[1] + left2Right*((right_Tvalue[1] - left_Tvalue[1]-faceTGrad[1].dot_product(left2Right))/left2Right.L2Square());

            sigma_k = face_f1 * sigma_k1 + (1-face_f1)*sigma_k2;
            sigma_w = face_f1 * sigma_w1 + (1-face_f1)*sigma_w2;

            faceMul = vFaceMul[i][k];
            faceMut = vFaceMut[i][k];

            txx_k = (faceMul + sigma_k * faceMut) * faceTGrad[0][0];
            tyy_k = (faceMul + sigma_k * faceMut) * faceTGrad[0][1];
            txx_w = (faceMul + sigma_w * faceMut) * faceTGrad[1][0];
            tyy_w = (faceMul + sigma_w * faceMut) * faceTGrad[1][1];
            if(ndim==3)
            {
                tzz_k = (faceMul + sigma_k * faceMut) * faceTGrad[0][2];
                tzz_w = (faceMul + sigma_k * faceMut) * faceTGrad[1][2];
            }
            FaceViscousTFlux[i][2*k] = txx_k * norm[0] + tyy_k * norm[1];
            FaceViscousTFlux[i][2*k+1] = txx_w * norm[0] + tyy_w * norm[1];
            if(ndim==3)
            {
                FaceViscousTFlux[i][2*k] += tzz_k * norm[2];
                FaceViscousTFlux[i][2*k+1] += tzz_w * norm[2];
            }
            FaceViscousTFlux[i][2*k] *= curEdge.area();
            FaceViscousTFlux[i][2*k+1] *= curEdge.area();
        }
    }

}

template<int ndim>
double SSTkomegaModel<ndim>::rans_P_Solver(int meshind, int cellIndex)
{
    GeomElements::vector3d<ndim,double>** gradU = (GeomElements::vector3d<ndim,double>**)d_flow_strategy->getGradientPrimitive();
    double ransP = 0.0;
    auto det_ij = [=](int i ,int j)->int
    {
        if(i ==j)
        {
            return 1;
        }
        else
        {
            return 0;
        }
    };

    double div = gradU[meshind][2+nequ*cellIndex][0] + gradU[meshind][3+nequ*cellIndex][1];
    if(ndim==3)
    {
        div += gradU[meshind][4+nequ*cellIndex][2];
    }
    double** U = d_flow_strategy->getU();
    double rho = U[meshind][nequ*cellIndex];

    double ddls0,ddls1,ddls2;
    ddls0=2.0/3.0*((gradU[meshind][2+nequ*cellIndex][0]-gradU[meshind][3+nequ*cellIndex][1])*(gradU[meshind][2+nequ*cellIndex][0]-gradU[meshind][3+nequ*cellIndex][1]));		
	ddls1=(gradU[meshind][2+nequ*cellIndex][1]+gradU[meshind][3+nequ*cellIndex][0])*(gradU[meshind][2+nequ*cellIndex][1]+gradU[meshind][3+nequ*cellIndex][0]);		
	ddls2=2.0/3.0*rho*turbulenceCellValue[meshind][2*cellIndex]*(gradU[meshind][2+nequ*cellIndex][0]+gradU[meshind][3+nequ*cellIndex][1]);	
    ransP=(ddls0+ddls1)*vCellMut[meshind][cellIndex] -ddls2;
    // for(int i = 0;i<ndim;i++)
    // {
    //     for(int j = 0;j<ndim;j++)
    //     {
    //         double part_P = vCellMut[meshind][cellIndex] * 
    //         (gradU[meshind][i+2+nequ*cellIndex][j]+gradU[meshind][j+2+nequ*cellIndex][i]
    //         - 2.0 / 3.0 *det_ij(i,j)*div) - 2.0 /3.0* det_ij(i,j)*rho*turbulenceCellValue[meshind][2*cellIndex];
    //         ransP += part_P * gradU[meshind][i+2+nequ*cellIndex][j];
    //     }
    // }

    return ransP;


}

template<int ndim>
void SSTkomegaModel<ndim>::SolveTotalTFlux()
{
    double Qflux[2];
    double** U  =d_flow_strategy->getU();
    for(int i = 0;i<d_nmesh;i++)
    {
        auto& curBlk = d_hder->blk2D[i];
        for(int k = 0;k<d_hder->nCells(i);k++)
        {
            auto& curCell = curBlk.d_localCells[k];
            double f1 = cellF1[i][k];
            double f2 = cellF2[i][k];
            double rdiv = turbulenceCellGrad[i][2*k].dot_product(turbulenceCellGrad[i][2*k+1]);
            double Cw = f1 * Cw1 + (1.0 - f1)* Cw2;
            double beta = f1 * beta1 + (1.0 - f1)* beta2;
            double rho = U[i][nequ*k];

            double rans_P = rans_P_Solver(i,k);

            Qflux[0] = rans_P - _beta * rho * turbulenceCellValue[i][2*k] * turbulenceCellValue[i][2*k+1];
            Qflux[1] = Cw * rho * rans_P / vCellMut[i][k] - beta * rho 
            * turbulenceCellValue[i][2*k+1] * turbulenceCellValue[i][2*k+1] + 2 *(1-f1) *rho * sigma_w2/turbulenceCellValue[i][2*k+1]*rdiv;
            // std::cout<<Qflux[1]<<" = "<< Cw << " * " << rho << " * " << rans_P <<" / "<<vCellMut[i][k] <<" - "
            // <<beta<<" * "<<rho << " * "<<turbulenceCellValue[i][2*k+1] << " * "<<turbulenceCellValue[i][2*k+1]<<" + "
            // <<" 2 * (1 - "<<f1<<" ) * "<< rho <<" * "<<sigma_w2<<" / "<<turbulenceCellValue[i][2*k+1]<<" * "<< rdiv<<'\n';
            CellTotalTFlux[i][2*k] =  Qflux[0]*curCell.volume();
            CellTotalTFlux[i][2*k+1] =  Qflux[1]*curCell.volume();
        }
        for(int k = 0;k<d_hder->nEdges(i);k++)
        {
            auto& curEdge = curBlk.d_localEdges[k];
            int lC = curEdge.lCInd();
            int rC = curEdge.rCInd();

            CellTotalTFlux[i][2*lC] += FaceConvectiveTFlux[i][2*k] - FaceViscousTFlux[i][2*k];
            CellTotalTFlux[i][2*lC+1] += FaceConvectiveTFlux[i][2*k+1] - FaceViscousTFlux[i][2*k+1];
            if(rC>=0)
            {
                CellTotalTFlux[i][2*rC] -= FaceConvectiveTFlux[i][2*k] - FaceViscousTFlux[i][2*k];
                CellTotalTFlux[i][2*rC+1] -= FaceConvectiveTFlux[i][2*k+1] - FaceViscousTFlux[i][2*k+1];
            }
        }
    }

}


template<int ndim>
void SSTkomegaModel<ndim>::writeCellData(const std::string &filename, int proc, int meshTag, double **dcelldata)
{
    
    std::string local_filename = filename + "_"+std::to_string(meshTag+1)+"_" + std::to_string(proc+1);
    std::ofstream fout;
    fout.open(local_filename,std::ios::out);
    
    
    if (!fout.is_open())
    {
        std::cout << "Open file failure\n";
        return;
    }
    
    fout << "TITLE =\"Tioga output\"\n";
    fout << "VARIABLES=\"X\",\"Y\",";
    if(d_dim==3)
    {
        fout << "\"Z\",";
    }
    fout << "\"K\",\"w\"\n";
    int nnodes = d_hder->nPoints(meshTag);
    int ncells = d_hder->nCells(meshTag);
    fout << "ZONE T=\"VOL_MIXED\",N=" << nnodes << " E=" << ncells;
    fout << " ET = QUADRILATERAL, F = FEBLOCK\n";
    fout << "VARLOCATION =  (";
    
    for (int i = 0; i < d_dim; i++)
    {
        fout << i + 1 << "=NODAL,";
    }
    for (int i = 0; i < d_NEQU -1; i++)
    {
        fout << d_dim + 1 + i << "=CELLCENTERED,";
    }
    fout <<d_dim + d_NEQU <<"=CELLCENTERED)\n";
    for (int j = 0; j < d_dim; j++)
    {
        for (int i = 0; i < nnodes; i++)
        {
            fout << d_hder->blk2D[meshTag].d_localPoints[i][j] << '\n';
        }
    }
    for (int i = d_dim; i < d_dim + d_NEQU; i++)
    {
        for (int j = 0; j < d_hder->nCells(meshTag); j++)
        {
            fout << dcelldata[meshTag][i - d_dim+d_NEQU*j] << '\n';
        }
    }
    double fs_velocityMag = (fsnd_velocity_components[0]*fsnd_velocity_components[0]+fsnd_velocity_components[1]*fsnd_velocity_components[1]);

    for (int i = 0; i < ncells; i++)
    {
        for (int j = 0; j < d_hder->blk2D[meshTag].d_localCells[i].size(); j++)
        {
            fout << d_hder->blk2D[meshTag].d_localCells[i].pointInd(j) + 1 << '\t';
        }
        for (int j = 0; j < 4 - d_hder->blk2D[meshTag].d_localCells[i].size(); j++)
        {
            fout << d_hder->blk2D[meshTag].d_localCells[i].pointInd(d_hder->blk2D[meshTag].d_localCells[i].size() - 1) + 1 << '\t';
        }
        fout << '\n';
    }
    fout.close();
}