#include "LU_SGS_Strategy.h"
#include "UnstructIntegrator.h"
#include "TopologyHolderStrategy.h"
#include <algorithm>
#include <iostream>
#include <mpi/mpi.h>
#include "toolBox/edge3d_int.h"
LUSGSStrategy::LUSGSStrategy(UnstructTopologyHolder *hder, TopologyHolderStrategy* hder_strategy)
    :TimeStrategy(hder,hder_strategy)
{
    fs_Pr = d_hder_strategy->fs_Pr;
    fs_Prt = d_hder_strategy->fs_Prt;

    ForwardSweep_exceptions = new int*[d_nmesh];
    nForwardSweep_exceptions = new int[d_nmesh];
    ForwardSweep_Only = new int*[d_nmesh];
    nForwardSweep_Only = new int[d_nmesh];

    auto ForwardSweep_excluded_cells = new std::set<int,std::less<int>>[d_nmesh];
    for (int i = 0;i<d_nmesh;i++)
    {
        auto& curBlk = d_hder->blk2D[i];
        ForwardSweep_excluded_cells[i].clear();
        for(auto iter = d_hder->sendBegin(i);iter!=d_hder->sendEnd(i);iter++)
        {
            if(iter->d_layer_level<=1)
            {
                ForwardSweep_excluded_cells[i].insert(iter->d_localId);
            }
        }
        for(auto iter = d_hder->recvBegin(i);iter!=d_hder->recvEnd(i);iter++)
        {
            if(iter->d_layer_level<=1)
            {
                ForwardSweep_excluded_cells[i].insert(iter->d_localId);
            }
        }
        nForwardSweep_Only[i] = ForwardSweep_excluded_cells[i].size();//Add -1 and nCell
        ForwardSweep_Only[i] =  new int[nForwardSweep_Only[i]];
        int temp_int = 0;
        
        for(auto iter = ForwardSweep_excluded_cells[i].begin();iter != ForwardSweep_excluded_cells[i].end();iter++)
        {
            ForwardSweep_Only[i][temp_int++] = *iter;
        }
        
        for(auto iter = d_hder->recvBegin(i);iter!=d_hder->recvEnd(i);iter++)
        {
            if(iter->d_layer_level>1)
            {
                ForwardSweep_excluded_cells[i].insert(iter->d_localId);
            }
        }
        nForwardSweep_exceptions[i] = ForwardSweep_excluded_cells[i].size()+2;
        ForwardSweep_exceptions[i] = new int[ForwardSweep_excluded_cells[i].size()+2];
        temp_int = 0;
        ForwardSweep_exceptions[i][temp_int++] = -1;
        for(auto iter = ForwardSweep_excluded_cells[i].begin();iter != ForwardSweep_excluded_cells[i].end();iter++)
        {
            ForwardSweep_exceptions[i][temp_int++] = *iter;
        }
        ForwardSweep_exceptions[i][temp_int++] = d_hder->nCells(i);
    }
    delete[] ForwardSweep_excluded_cells;
    
    d_hder_strategy->test_unwantedSweep(ForwardSweep_Only,nForwardSweep_Only);

    d_diagOperator = new double*[d_nmesh];
    d_deltaW = new double*[d_nmesh];
    d_deltaW1 = new double*[d_nmesh];
    W_scratch = new double*[d_nmesh];

    d_stableDt = new double*[d_nmesh];

    for(int i = 0 ;i<d_nmesh;i++)
    {
        d_diagOperator[i] = new double[d_NEQU*d_hder->nCells(i)];
        d_deltaW[i] = new double[d_NEQU*d_hder->nCells(i)];
        d_deltaW1[i] = new double[d_NEQU*d_hder->nCells(i)];
        W_scratch[i] = new double[d_NEQU*d_hder->nCells(i)];

        d_stableDt[i] = new double[d_hder->nCells(i)];
    }
}
void LUSGSStrategy::initialize()
{
    d_hder_strategy->initializeData();
    
    const double Gamma = 1.4;
    double **U = d_hder_strategy->getU();
    for(int i = 0;i<d_nmesh;i++)
    {
        for(int k = 0; k< d_hder->nCells(i); k++)
        {
            double VMSquare = 0;
            for (int l = 0; l < d_dim; l++)
            {
                VMSquare += U[i][l + 2+d_NEQU*k] * U[i][l + 2+d_NEQU*k];
                W_scratch[i][l + 2+d_NEQU*k] = U[i][0+d_NEQU*k] * U[i][l + 2+d_NEQU*k];
            }
            W_scratch[i][0+d_NEQU*k] = U[i][0+d_NEQU*k];
            W_scratch[i][1+d_NEQU*k] = U[i][1+d_NEQU*k] / (Gamma - 1) + 0.5 * U[i][0+d_NEQU*k] * VMSquare;
        }
    }
}

void LUSGSStrategy::preprocessUpdate()
{

}

void LUSGSStrategy::singleStep(int curStage)
{
    Update();
    // postprocessUpdate();
}

void LUSGSStrategy::singleStepSerial(int curStep)
{
    UpdateSerial();
    postprocessUpdate();
    
}

void LUSGSStrategy::postprocessUpdate()
{
    //Parallel,we only count the send buffer cells and inner region cells.
    int cur_proc;
    MPI_Comm_rank(MPI_COMM_WORLD,&cur_proc);

    //First clear all the cells at recv, check the TopologyHolderStrategy for accessing these cells.
    for(int i = 0;i<d_nmesh;i++)
    {
        for(auto iter = d_hder->recvBegin(i);iter!=d_hder->recvEnd(i);iter++)
        {
            d_deltaW[i][d_NEQU*iter->d_localId] = 0;
        }
        double residual = 0;
        for(int k = 0;k<d_hder->nCells(i);k++)
        {
            residual += d_deltaW[i][d_NEQU*k]*d_deltaW[i][d_NEQU*k];
        }
        std::cout<<"residual at mesh "<<i <<" on proc "<<cur_proc<<" is "<<residual<<'\n';
    }

    
}


bool LUSGSStrategy::converged()
{

    

}

void LUSGSStrategy::preprocessLUSGS()
{

    for(int i = 0;i<d_nmesh;i++)
    {
        auto& curBlk = d_hder->blk2D[i];
        for(int j = 0;j<nForwardSweep_Only[i];j++)
        {
            int cellId = ForwardSweep_Only[i][j];
            auto& curCell = curBlk.d_localCells[cellId];
            for(int k = 0;k<curCell.edge_size();k++)
            {
                int edgeInd = curCell.edgeInd(k);
                auto& curEdge = curBlk.d_localEdges[edgeInd];
                int lC = curEdge.lCInd();
                int rC = curEdge.rCInd();

                if(lC==cellId)
                {
                    curEdge.setLeft(lC + d_hder->nCells(i));
                }
                else if(rC==cellId)
                {
                    curEdge.setRight(rC + d_hder->nCells(i));
                }
                else
                {
                    throw std::runtime_error("Left right not found at preprocess\n");
                }
            }
        }
    }
}

void LUSGSStrategy::Update()
{
    d_hder_strategy->preprocessAdvance(0);
    

    d_hder_strategy->SolveTime(d_stableDt,d_CFL_number);
    //In the middle is the routine for LUSGS

    UpdateConservativeForRecvCells();

    SolveDiagOperator();
    // Forward Sweeop Procedure:

    // 1.For all {SendCell_lv1} + {RecvCell_lv1}, change its edge's left right in preprocessLUSGS.
    preprocessLUSGS();

    // 2.For all {localCell} - {SendCell_lv1}, using its original ranking perform forward Sweep.
    // 3.SendRecv DWS computed in (Procedure 2.) from {SendCell_lv2} to {RecvCell_lv2}.
    // 4.For all {SendCell_lv1} + {RecvCell_lv1}, using its updated ranking perform forward Sweep.
    SolveForwardSweep();

    // 5.For all {SendCell_lv1} + {RecvCell_lv1}, using its updated ranking perform backward Sweep.
    // 6.For all {localCell} - {SendCell_lv1}, using its updated ranking perform backward Sweep.
    SolveBackwardSweep();

    // 7.For all {SendCell_lv1} + {RecvCell_lv1}, change its edge's left right in postprocessLUSGS.
    postprocessLUSGS();
    
    // Update W_scratch and U.
    UpdatePrimitiveVariable();
    
}

void LUSGSStrategy::UpdateSerial()
{
    d_hder_strategy->preprocessAdvanceSerial(0);
    d_hder_strategy->SolveTime(d_stableDt,d_CFL_number);
    // std::cin.get();
    SolveDiagOperator();
    // std::cout<<"Entering ForwardSweep\n";
    // std::cin.get();
    // std::cin.get();
    SolveForwardSweepSerial();
    // std::cout<<"Entering BackwardSweep\n";
    // std::cin.get();
    SolveBackwardSweepSerial();
    // std::cin.get();
    // std::cout<<"Entering UpdatePrimitiveVariable\n";
    UpdatePrimitiveVariable();
}

void LUSGSStrategy::UpdatePrimitiveVariable()
{
    const double Gamma = 1.4;
    double** U = d_hder_strategy->getU();
    for(int i = 0;i<d_nmesh;i++)
    {
        for(int k = 0;k<d_hder->nCells(i);k++)
        {
            for(int j = 0;j<d_NEQU;j++)
            {
                // std::cout<<d_deltaW[i][j+d_NEQU*k]<<'\n';
                // std::cin.get();
                W_scratch[i][j+d_NEQU*k] = W_scratch[i][j+d_NEQU*k] + d_deltaW[i][j+d_NEQU*k];
            }
            U[i][0+d_NEQU*k] = W_scratch[i][0+d_NEQU*k];
            double VMSquare = 0;
            for (int l = 0; l < d_dim; l++)
            {
                U[i][l + 2+d_NEQU*k] = W_scratch[i][l + 2+d_NEQU*k] / U[i][0+d_NEQU*k];
                VMSquare += U[i][l + 2+d_NEQU*k] * U[i][l + 2+d_NEQU*k];            
            }
            U[i][1+d_NEQU*k] = (Gamma - 1) * (W_scratch[i][1+d_NEQU*k] - 0.5 * U[i][0+d_NEQU*k] * VMSquare);
        }
    }
}

void LUSGSStrategy::postprocessLUSGS()
{
    for(int i = 0;i<d_nmesh;i++)
    {
        auto& curBlk = d_hder->blk2D[i];
        for(int j = 0;j<nForwardSweep_Only[i];j++)
        {
            int cellId = ForwardSweep_Only[i][j];
            auto& curCell = curBlk.d_localCells[cellId];
            for(int k = 0;k<curCell.edge_size();k++)
            {
                int edgeInd = curCell.edgeInd(k);
                auto& curEdge = curBlk.d_localEdges[edgeInd];
                int lC = curEdge.lCInd();
                int rC = curEdge.rCInd();

                if(lC==(cellId + d_hder->nCells(i)))
                {
                    curEdge.setLeft(lC - d_hder->nCells(i));
                }
                else if(rC==(cellId+d_hder->nCells(i)))
                {
                    curEdge.setRight(rC - d_hder->nCells(i));
                }
                else
                {
                    throw std::runtime_error("Left right not found at postProcessing\n");
                }
            }
            
        }
    }
}

void LUSGSStrategy::SolveDiagOperator()
{
    SolveDiag(d_diagOperator,d_stableDt);
}

void LUSGSStrategy::SolveForwardSweep()
{

    // 2.For all {localCell} - {SendCell_lv1}, using its original ranking perform forward Sweep.
    SolveLocalForwardSweep(d_diagOperator,W_scratch,d_deltaW1,ForwardSweep_exceptions,nForwardSweep_exceptions);
    // 3.SendRecv DWS computed in (Procedure 2.) from {SendCell_lv2} to {RecvCell_lv2}.

    if(d_hder_strategy->cur_proc==-1)std::cout<<"Finishing Loacal forward Sweep\n";
    d_hder_strategy->RemoteCellCommunication(d_deltaW1);

    d_hder_strategy->NearCellCommunication(d_hder_strategy->getResidual());

    if(d_hder_strategy->cur_proc==-1)std::cout<<"FInishing Remote cell communication\n";
    // 4.For all {SendCell_lv1} + {RecvCell_lv1}, using its updated ranking perform forward Sweep.
    SolveBoundaryForwardSweep(d_diagOperator,W_scratch,d_deltaW1,ForwardSweep_Only,nForwardSweep_Only);
    if(d_hder_strategy->cur_proc==-1)std::cout<<"FInishing boundary forward Sweep\n";
}

void LUSGSStrategy::SolveBackwardSweep()
{
    // 5.For all {SendCell_lv1} + {RecvCell_lv1}, using its updated ranking perform backward Sweep.
    SolveBoundaryBackwardSweep(d_diagOperator,W_scratch,d_deltaW1,d_deltaW,ForwardSweep_Only,nForwardSweep_Only);
    if(d_hder_strategy->cur_proc==-1)std::cout<<"Finishing boundary backward Sweep\n";
    // 6.For all {localCell} - {SendCell_lv1}, using its updated ranking perform backward Sweep.
    SolveLocalBackwardSweep(d_diagOperator,W_scratch,d_deltaW1,d_deltaW,ForwardSweep_exceptions,nForwardSweep_exceptions);
    if(d_hder_strategy->cur_proc==-1)std::cout<<"Finishing local backward Sweep\n";
}

void LUSGSStrategy::SolveForwardSweepSerial()
{
    SolveForwardSweep(d_diagOperator,W_scratch,d_deltaW1);
}

void LUSGSStrategy::SolveBackwardSweepSerial()
{
    SolveBackwardSweep(d_diagOperator,W_scratch,d_deltaW1,d_deltaW);

}

void LUSGSStrategy::UpdateConservativeForRecvCells()
{
    double Gamma = 1.4;
    double** U = d_hder_strategy->getU();
    for(int i = 0;i<d_nmesh;i++)
    {
        for(auto iter = d_hder->recvBegin(i);iter!=d_hder->recvEnd(i);iter++)
        {
            int cellId = iter->d_localId;
            W_scratch[i][0+d_NEQU*cellId] = U[i][0+d_NEQU*cellId];
            W_scratch[i][2+d_NEQU*cellId] = U[i][0+d_NEQU*cellId]*U[i][2+d_NEQU*cellId];
            W_scratch[i][3+d_NEQU*cellId] = U[i][0+d_NEQU*cellId]*U[i][3+d_NEQU*cellId];
            W_scratch[i][1+d_NEQU*cellId] = U[i][1+d_NEQU*cellId]/(Gamma - 1) 
            + 0.5*U[i][0+d_NEQU*cellId]*(powf64(U[i][2+d_NEQU*cellId],2)+powf64(U[i][3+d_NEQU*cellId],2));
        }
    }
}


void LUSGSStrategy::SolveDiag(double** de_diag, double** dt)
{
    double** spectrum_convec = d_hder_strategy->getSpectrumConvective();
    double OMEGAN = 1.5;
    for(int i = 0;i<d_nmesh;i++)
    {
        auto& curBlk = d_hder->blk2D[i];
        for(int k = 0;k<d_hder->nCells(i);k++)
        {
            auto& curCell = curBlk.d_localCells[k];
            double diag = curCell.volume()/dt[i][k] + 0.5*OMEGAN*spectrum_convec[i][k];
            if(!d_hder_strategy->isInvicid())
            {
                double ** spectrum_viscous = d_hder_strategy->getSpectrumViscous();
                // std::cout<<"diag before "<<diag<<'\t';
                diag += spectrum_viscous[i][k];
                // std::cout<<"diag after "<<diag<<'\n';
            }
            for(int j = 0;j<d_NEQU;j++)
            {
                de_diag[i][j+d_NEQU*k] = diag;
            }
        }
    }
}



void LUSGSStrategy::SolveDF(int curMesh,int curCell,int curEdge,
            double** W_scratch ,double** deltaW0,double* DF,
            int ForB,int cellOffset)
{
    auto& edg = d_hder->blk2D[curMesh].d_localEdges[curEdge];
    double Wp[d_NEQU];
    double W[d_NEQU];
    double dW[d_NEQU];
    double W0[d_NEQU];
    double Fcp[d_NEQU];
    double Fc[d_NEQU];

    double rA;
    double Fv[d_NEQU];
    double DFc[d_NEQU];
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
            for(int j = 0;j<d_NEQU;j++)
            {
                W0[j] = W_scratch[curMesh][j+d_NEQU*rC];
                dW[j] = deltaW0[curMesh][j+d_NEQU*rC];
                Wp[j] =W0[j] + dW[j];
                W[j] = W0[j];
            }
            //Get convectiveFlux
            GeomElements::vector3d<2,double> norm = edg.normal_vector();
            ConservativeParam2ConvectiveFlux(W,Fc,norm);
            ConservativeParam2ConvectiveFlux(Wp,Fcp,norm);
            SolveViscousFlux(curMesh,curCell,curEdge,rC,dW,Fv);
            

            //Newly added:compute Viscous FLux
            //i is lc and j is rc
            //Get left cell's Vn
            

            for(int j = 0;j<d_NEQU;j++)
            {
                DFc[j] = Fcp[j] - Fc[j] ;
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
            for(int j = 0;j<d_NEQU;j++)
            {
                W0[j] = W_scratch[curMesh][j+d_NEQU*lC];
                dW[j] = deltaW0[curMesh][j+d_NEQU*lC];
                Wp[j] = W0[j] + dW[j];
                W[j] = W0[j];
            }
            //Get convectiveFlux
            GeomElements::vector3d<2,double> norm = edg.normal_vector();
            ConservativeParam2ConvectiveFlux(W,Fc,norm);
            ConservativeParam2ConvectiveFlux(Wp,Fcp,norm);
            SolveViscousFlux(curMesh,curCell,curEdge,lC,dW,Fv);
            for(int j = 0;j<d_NEQU;j++)
            {
                DFc[j] =  Fc[j] - Fcp[j];
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
        // int cur_proc;
        // MPI_Comm_rank(MPI_COMM_WORLD,&cur_proc);
        // if(cur_proc==0)
        // {
        std::cout<<curCell<<" with offset "<< cellOffset <<" not found in "<<lC<<" and "<<rC<<'\n';
        throw std::runtime_error("Cur cell is not on either side,impossible\n");
    // }
    }
}

void LUSGSStrategy::SolveViscousFlux(int curMesh,int curCell,int curEdge,int anotherCell,double* detlaW,double* Fv)
{
    const double Gamma = 1.4;
    const double OMEGAN = 1.5;
    auto& curBlk = d_hder->blk2D[curMesh];
    double rA;
    double** U = d_hder_strategy->getU();
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
    GeomElements::vector3d<2,double> curVelocity(U[curMesh][2+d_NEQU*anotherCell],U[curMesh][3+d_NEQU*anotherCell]);
    GeomElements::vector3d<2,double> norm = curE.normal_vector();
    double Vn = curVelocity.dot_product(norm);
    double c = sqrtf64(Gamma * U[curMesh][1+d_NEQU*anotherCell] / U[curMesh][0+d_NEQU*anotherCell]);
    
    rA  = OMEGAN*(fabsf64(Vn)+c);
    if(!d_hder_strategy->isInvicid())
    {
        // std::cout<<"NOt invicid\n";
        // std::cout<<"rA before"<<rA<<'\t';
        rA += std::max<double>(4/(3*U[curMesh][0+d_NEQU*anotherCell]),Gamma/U[curMesh][0+d_NEQU*anotherCell])*
        (d_hder_strategy->getMuOverPrCell(curMesh,curCell))/sqrtf64(centerLine.L2Square());
        // std::cout<<"rA after"<<rA<<"\n";
        // std::cout<<std::max<double>(4/(3*U[curMesh][0+d_NEQU*anotherCell]),Gamma/U[curMesh][0+d_NEQU*anotherCell])<<"*"
        // <<(d_hder_strategy->getMuOverPrCell(curMesh,curCell))<<"/"
        // <<centerLine.L2Square()<<"="<<std::max<double>(4/(3*U[curMesh][0+d_NEQU*anotherCell]),Gamma/U[curMesh][0+d_NEQU*anotherCell])*
        // (d_hder_strategy->getMuOverPrCell(curMesh,curCell))/sqrtf64(centerLine.L2Square())<<'\n';
        // std::cin.get();
    }
    for(int j = 0;j<d_NEQU;j++)
    {
        Fv[j] = rA*detlaW[j];
    }
}

//The skip point in here are level 1 sendcells and level 1,2 recvcells.
void LUSGSStrategy::SolveLocalForwardSweep(double** diag,
                            double** W_scratch,double** deltaW1,
                            int** skip_point,int* nskip)
{
    int ForB = 1;
    double** Residual = d_hder_strategy->getResidual();
    
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
                    deltaW1[i][j+d_NEQU*cellId] = (-Residual[i][j+d_NEQU*cellId] - dFi[j])/diag[i][j+d_NEQU*cellId];
                }
            }
        }
    }
}
//The skip point in here are level 1 cells.
void LUSGSStrategy::SolveBoundaryForwardSweep(double** diag,
                            double** W_scratch,double** deltaW1,
                            int** skip_point,int * nskip)
{
    int ForB = 1;
    double** Residual = d_hder_strategy->getResidual();
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
                deltaW1[i][j+d_NEQU*cellId] = (-Residual[i][j+d_NEQU*cellId]-dFi[j])/diag[i][j+d_NEQU*cellId];
;           }   
        }
    }
}
void LUSGSStrategy::SolveBoundaryBackwardSweep(double** diag,
                                double** W_scratch,double** deltaW1,double** deltaW,
                                int ** skip_point,int * nskip)
{
    int ForB = 0;
    double dFi[d_NEQU];
    double dF[d_NEQU];
    for(int i = 0;i< d_nmesh;i++)
    {
        int cellIdOffset = d_hder->nCells(i);
        auto& curBlk = d_hder->blk2D[i];
        int cur_proc;

        for(int k = nskip[i]-1;k>=0;k--)
        // for(int k = 0;k<nskip[i];k++)
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
                deltaW[i][j+d_NEQU*cellId] = deltaW1[i][j+d_NEQU*cellId] - dFi[j]/diag[i][j+d_NEQU*cellId];
            }
        }
    }
}
void LUSGSStrategy::SolveLocalBackwardSweep(double** diag,
                            double** W_scratch,double** deltaW1,double** deltaW,
                            int** skip_point,int* nskip)
{
    int ForB = 0;
    double dFi[d_NEQU];
    double dF[d_NEQU];
    for(int i = 0;i< d_nmesh;i++)
    {
        int cellIdOffset = 0;
        auto& curBlk = d_hder->blk2D[i];
        for(int k = nskip[i]-2;k>=0;k--)
        // for(int k = 0;k<nskip[i]-1;k++)
        {
            for(int l = skip_point[i][k+1]-1;l>skip_point[i][k];l--)
            // for(int l = skip_point[i][k]+1;l<skip_point[i][k+1];l++)
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
                    deltaW[i][j+d_NEQU*cellId] = deltaW1[i][j+d_NEQU*cellId] - dFi[j]/diag[i][j+d_NEQU*cellId];
                }
            }
        }
    }
}

void LUSGSStrategy::SolveDfSimple(double *W,double* DW,
                    const GeomElements::cell3d<2>&anotherCell,
                    const GeomElements::edge3d<2>& curEdge,
                    int LorR, double* dF)
{
    double Wp[5];
    double Fc[5];
    double Fcp[5];
    double Ds = curEdge.area();
    GeomElements::vector3d<2,double> norm= curEdge.normal_vector();
    norm[0] = -norm[0];
    norm[1] = -norm[1];
    for(int j = 0;j<d_NEQU;j++)
    {
        Wp[j] = W[j] + DW[j];
    }
    ConservativeParam2ConvectiveFlux(W,Fc,norm);
    ConservativeParam2ConvectiveFlux(Wp,Fcp,norm);

    for(int j = 0;j<d_NEQU;j++)
    {
        dF[j] = 0.5*(Fcp[j] - Fc[j])*Ds;
    }
}
// void LUSGSStrategy::SolveForwardSweep(double** diag,
//                     double** W_scratch,double** deltaW1)
// {
//     double** Residual = d_hder_strategy->getResidual();
//     double dFi[d_NEQU];
//     double dF[d_NEQU];
//     double anotherW[d_NEQU];
//     double anotherDWs[d_NEQU];
//     for(int i = 0;i<d_nmesh;i++)
//     {
//         auto& curBlk = d_hder->blk2D[i];
//         for(int k = 0;k<d_hder->nCells(i);k++)
//         {
//             auto& curCell = curBlk.d_localCells[k];
//             for(int j = 0;j<d_NEQU;j++)
//             {
//                 dFi[j] = 0;
//             }
//             for(int j = 0;j<curCell.edge_size();j++)
//             {
//                 auto& curEdge  = curBlk.d_localEdges[curCell.edgeInd(j)];
//                 int lC = curEdge.lCInd();
//                 int rC = curEdge.rCInd();
//                 if(k==lC)//Left Side
//                 {
//                     if(rC>=0 and rC<k)
//                     {
//                         auto& anotherCell = curBlk.d_localCells[rC];
//                         for(int l = 0;l<d_NEQU;l++)
//                         {
//                             anotherW[l] = W_scratch[i][l+d_NEQU*rC];
//                             anotherDWs[l] = deltaW1[i][l+d_NEQU*rC];
//                         }
//                         SolveDfSimple(anotherW,anotherDWs,anotherCell,curEdge,1,dF);
//                         for(int l =0;l<d_NEQU;l++)
//                         {
//                             dFi[l] += dF[l];
//                         }
//                     }
//                 }
//                 else if(k==rC)
//                 {
//                     if(lC<k)
//                     {
//                         auto& anotherCell = curBlk.d_localCells[lC];
//                         for(int l = 0;l<d_NEQU;l++)
//                         {
//                             anotherW[l] = W_scratch[i][l+d_NEQU*lC];
//                             anotherDWs[l] = deltaW1[i][l+d_NEQU*lC];
//                         }
//                         SolveDfSimple(anotherW,anotherDWs,anotherCell,curEdge,-1,dF);
//                         for(int l =0;l<d_NEQU;l++)
//                         {
//                             dFi[l] += dF[l];
//                         }                        
//                     }
//                 }
//                 else
//                 {
//                     std::cout<<"Something Wronmg\n";
//                     std::cin.get();
//                 }
//             }
//             for(int j = 0;j<d_NEQU;j++)
//             {
//                 deltaW1[i][j+d_NEQU*k] = (-Residual[i][j+d_NEQU*k]-dFi[j])/diag[i][j+d_NEQU*k];
//             }
//         }
//     }
// }

void LUSGSStrategy::SolveForwardSweep(double** diag,
                    double** W_scratch,double** deltaW1)
{
    int ForB = 1;//1 indicates forward sweep,that is dF=0 when cur < another
    double** Residual = d_hder_strategy->getResidual();
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
                    //  std::cout<<"dFi[j]"<<dFi[l]<<'\n';
                }
            }
            for(int j = 0;j<d_NEQU;j++)
            {
                deltaW1[i][j+d_NEQU*k] = (-Residual[i][j+d_NEQU*k]-dFi[j])/diag[i][j+d_NEQU*k];

                // std::cout<<"deltaW1[i][j+d_NEQU*k]"<<deltaW1[i][j+d_NEQU*k]<<'\n';
                // std::cout<<"Residual[i][j+d_NEQU*k]"<<Residual[i][j+d_NEQU*k]<<'\n';
                // std::cout<<"dFi[j]"<<dFi[j]<<'\n';
            }
            // if(k>10000)
            
            // std::cin.get();
        }
    }
}

// void LUSGSStrategy::SolveBackwardSweep(double** diag,
//                         double** W_scratch,double** deltaW1,double** deltaW)
// {
//     double dFi[d_NEQU];
//     double dF[d_NEQU];
//     double anotherW[d_NEQU];
//     double anotherDW[d_NEQU];
//     for(int i = 0;i<d_nmesh;i++)
//     {
//         auto& curBlk = d_hder->blk2D[i];
//         for(int k = 0;k<d_hder->nCells(i);k++)
//         {
//             auto& curCell = curBlk.d_localCells[k];
//             for(int j = 0;j<d_NEQU;j++)
//             {
//                 dFi[j] = 0;
//             }
//             for(int j = 0;j<curCell.edge_size();j++)
//             {
//                 auto& curEdge  = curBlk.d_localEdges[curCell.edgeInd(j)];
//                 int lC = curEdge.lCInd();
//                 int rC = curEdge.rCInd();
//                 if(k==lC)//Left Side
//                 {
//                     if(rC>=0 and rC>k)
//                     {
//                         auto& anotherCell = curBlk.d_localCells[rC];
//                         for(int l = 0;l<d_NEQU;l++)
//                         {
//                             anotherW[l] = W_scratch[i][l+d_NEQU*rC];
//                             anotherDW[l] = deltaW[i][l+d_NEQU*rC];
//                         }
//                         SolveDfSimple(anotherW,anotherDW,anotherCell,curEdge,1,dF);
//                         for(int l =0;l<d_NEQU;l++)
//                         {
//                             dFi[l] += dF[l];
//                         }
//                     }
//                 }
//                 else if(k==rC)
//                 {
//                     if(lC>k)
//                     {
//                         auto& anotherCell = curBlk.d_localCells[lC];
//                         for(int l = 0;l<d_NEQU;l++)
//                         {
//                             anotherW[l] = W_scratch[i][l+d_NEQU*lC];
//                             anotherDW[l] = deltaW[i][l+d_NEQU*lC];
//                         }
//                         SolveDfSimple(anotherW,anotherDW,anotherCell,curEdge,-1,dF);
//                         for(int l =0;l<d_NEQU;l++)
//                         {
//                             dFi[l] += dF[l];
//                         }                        
//                     }
//                 }
//                 else
//                 {
//                     std::cout<<"Something Wronmg\n";
//                     std::cin.get();
//                 }
//             }
//             for(int j = 0;j<d_NEQU;j++)
//             {
//                 deltaW[i][j+d_NEQU*k] = deltaW1[i][j+d_NEQU*k] - dFi[j]/diag[i][j+d_NEQU*k];
//             }
//         }
//     }
// }

void LUSGSStrategy::SolveBackwardSweep(double** diag,
                        double** W_scratch,double** deltaW1,double** deltaW)
{
    int ForB = 0;
    double dFi[d_NEQU]={0};
    double dF[d_NEQU]={0};
    for(int i = 0;i<d_nmesh;i++)
    {
        auto& curBlk = d_hder->blk2D[i];
        for(int k = d_hder->nCells(i)-1;k>=0;k--)
        // for(int k = 0;k<d_hder->nCells(i);k++)
        {
            // std::cout<<"k is"<<k<<" out of "<<d_hder->nCells(i)<<'\n';
            auto& curCell = curBlk.d_localCells[k];
            for(int j = 0;j<d_NEQU;j++)
            {
                dFi[j] = 0;
            }
            // std::cout<<"curCell.edge_size() is "<<curCell.edge_size()<<'\n';
            for(int j = 0;j<curCell.edge_size();j++)
            {
                // std::cout<<"edgeINd is "<<curCell.edgeInd(j)<<'\n';
                auto& curEdge = curBlk.d_localEdges[curCell.edgeInd(j)];
                // std::cout<<"both sides ind are :"<<curEdge.lCInd()<<" "<<curEdge.rCInd()<<'\n';
                SolveDF(i,k,curCell.edgeInd(j),W_scratch,deltaW,dF,ForB,0);
                for(int l = 0;l<d_NEQU;l++)
                {
                    dFi[l] += dF[l];
                }
            }
            for(int j = 0;j<d_NEQU;j++)
            {
                deltaW[i][j+d_NEQU*k] = deltaW1[i][j+d_NEQU*k] - dFi[j]/diag[i][j+d_NEQU*k];
                // std::cout<<"deltaW[i][j+d_NEQU*k]"<<deltaW[i][j+d_NEQU*k]<<'\n';
                // std::cout<<"dFi[j]"<<dFi[j]<<'\n';
                // std::cout<<"diag[i][j+d_NEQU*k]"<<diag[i][j+d_NEQU*k]<<'\n';
                
                // std::cin.get();
            }
        }
    }
}