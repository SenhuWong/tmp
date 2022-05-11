#include "LU_SGS_Strategy.h"
#include "UnstructIntegrator.h"
#include "TopologyHolderStrategy.h"
LUSGSStrategy::LUSGSStrategy(UnstructTopologyHolder *hder, TopologyHolderStrategy* hder_strategy)
    :TimeStrategy(hder,hder_strategy)
{
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
    //Do nothing cuz U is obtained from W_scratch.
    return;
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

void LUSGSStrategy::singleStep(int curStage)
{
    Update();
    postprocessUpdate();
}

void LUSGSStrategy::singleStepSerial(int curStep)
{
    UpdateSerial();
    postprocessUpdate();
    
}

void LUSGSStrategy::postprocessUpdate()
{
    return;
    d_hder_strategy->cleaning(W_scratch);
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
    
    SolveDiagOperator();

    SolveForwardSweepSerial();

    SolveBackwardSweepSerial();

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
    double** spectrum_cell = d_hder_strategy->getSpectrum();
    double OMEGAN = 5;
    for(int i = 0;i<d_nmesh;i++)
    {
        auto& curBlk = d_hder->blk2D[i];
        for(int k = 0;k<d_hder->nCells(i);k++)
        {
            auto& curCell = curBlk.d_localCells[k];
            double diag = curCell.volume()/dt[i][k] + 0.5*OMEGAN*spectrum_cell[i][k];
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
                W0[j] = W_scratch[curMesh][j+d_NEQU*rC];
                Wp[j] =W0[j] + deltaW0[curMesh][j+d_NEQU*rC];
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
                W0[j] = W_scratch[curMesh][j+d_NEQU*lC];
                Wp[j] = W0[j] + deltaW0[curMesh][j+d_NEQU*lC];
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
        throw std::runtime_error("Cur cell is not on either side,impossible\n");
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
;              }   
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
                    deltaW[i][j+d_NEQU*cellId] = deltaW1[i][j+d_NEQU*cellId] - dFi[j]/diag[i][j+d_NEQU*cellId];
                }
            }
        }
    }
}
void LUSGSStrategy::SolveForwardSweep(double** diag,
                    double** W_scratch,double** deltaW1)
{
    int ForB = 1;
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
                }
            }
            for(int j = 0;j<d_NEQU;j++)
            {
                deltaW1[i][j+d_NEQU*k] = (-Residual[i][j+d_NEQU*k]-dFi[j])/diag[i][j+d_NEQU*k];
            }
        }
    }
}
void LUSGSStrategy::SolveBackwardSweep(double** diag,
                        double** W_scratch,double** deltaW1,double** deltaW)
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
                deltaW[i][j+d_NEQU*k] = deltaW1[i][j+d_NEQU*k] - dFi[j]/diag[i][j+d_NEQU*k];
            }
        }
    }
}