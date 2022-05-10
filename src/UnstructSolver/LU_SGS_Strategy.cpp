#include "LU_SGS_Strategy.h"
#include "UnstructIntegrator.h"
#include "Euler2D.h"
LUSGSStrategy::LUSGSStrategy(UnstructTopologyHolder *hder, Euler2D* hder_strategy)
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

    d_diagOperator = new double**[d_nmesh];
    
    d_deltaW = new double**[d_nmesh];
    d_deltaW1 = new double**[d_nmesh];
    W_scratch = new double**[d_nmesh];
    d_stableDt = new double*[d_nmesh];

    for(int i = 0 ;i<d_nmesh;i++)
    {
        d_diagOperator[i] = new double*[d_NEQU];
        d_deltaW[i] = new double*[d_NEQU];
        d_deltaW1[i] = new double*[d_NEQU];
        W_scratch[i] = new double*[d_NEQU];

        d_stableDt[i] = new double[d_hder->nCells(i)];
        for(int j = 0;j<d_NEQU;j++)
        {
            d_diagOperator[i][j] = new double[d_hder->nCells(i)];
            d_deltaW[i][j] = new double[d_hder->nCells(i)];
            d_deltaW1[i][j] = new double[d_hder->nCells(i)];
            W_scratch[i][j] = new double[d_hder->nCells(i)];
        }
    }
}
void LUSGSStrategy::initialize()
{
    d_hder_strategy->initializeData();
    
    const double Gamma = 1.4;
    double ***U = d_hder_strategy->getU();
    for(int i = 0;i<d_nmesh;i++)
    {
        for(int k = 0; k< d_hder->nCells(i); k++)
        {
            double VMSquare = 0;
            for (int l = 0; l < d_dim; l++)
            {
                VMSquare += U[i][l + 2][k] * U[i][l + 2][k];
                W_scratch[i][l + 2][k] = U[i][0][k] * U[i][l + 2][k];
            }
            W_scratch[i][0][k] = U[i][0][k];
            W_scratch[i][1][k] = U[i][1][k] / (Gamma - 1) + 0.5 * U[i][0][k] * VMSquare;
        }
    }
}

void LUSGSStrategy::preprocessUpdate()
{
    //Do nothing cuz U is obtained from W_scratch.
    return;
    const double Gamma = 1.4;
    double ***U = d_hder_strategy->getU();
    for(int i = 0;i<d_nmesh;i++)
    {
        for(int k = 0; k< d_hder->nCells(i); k++)
        {
            double VMSquare = 0;
            for (int l = 0; l < d_dim; l++)
            {
                VMSquare += U[i][l + 2][k] * U[i][l + 2][k];
                W_scratch[i][l + 2][k] = U[i][0][k] * U[i][l + 2][k];
            }
            W_scratch[i][0][k] = U[i][0][k];
            W_scratch[i][1][k] = U[i][1][k] / (Gamma - 1) + 0.5 * U[i][0][k] * VMSquare;
        }
    }
}

void LUSGSStrategy::singleStep(int curStage)
{
    //preprocessUpdate();
    Update();
    //postprocessUpdate();
    //if(d_hder_strategy->cur_proc==-1)
    //std::cout<<"FInishing singleStep?\n";

}

void LUSGSStrategy::singleStepSerial(int curStep)
{
    //preprocessUpdate();
    UpdateSerial();
    //postprocessUpdate();
    
}

void LUSGSStrategy::postprocessUpdate()
{
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
    double*** U = d_hder_strategy->getU();
    for(int i = 0;i<d_nmesh;i++)
    {
        for(int k = 0;k<d_hder->nCells(i);k++)
        {
            for(int j = 0;j<d_NEQU;j++)
            {
                W_scratch[i][j][k] = W_scratch[i][j][k] + d_deltaW[i][j][k];
            }
            U[i][0][k] = W_scratch[i][0][k];
            double VMSquare = 0;
            for (int l = 0; l < d_dim; l++)
            {
                U[i][l + 2][k] = W_scratch[i][l + 2][k] / U[i][0][k];
                VMSquare += U[i][l + 2][k] * U[i][l + 2][k];            
            }
            U[i][1][k] = (Gamma - 1) * (W_scratch[i][1][k] - 0.5 * U[i][0][k] * VMSquare);
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
    d_hder_strategy->SolveDiag(d_diagOperator,d_stableDt);
}

void LUSGSStrategy::SolveForwardSweep()
{

    // 2.For all {localCell} - {SendCell_lv1}, using its original ranking perform forward Sweep.
    d_hder_strategy->SolveLocalForwardSweep(d_diagOperator,W_scratch,d_deltaW1,ForwardSweep_exceptions,nForwardSweep_exceptions);
    // 3.SendRecv DWS computed in (Procedure 2.) from {SendCell_lv2} to {RecvCell_lv2}.

    if(d_hder_strategy->cur_proc==-1)std::cout<<"Finishing Loacal forward Sweep\n";
    d_hder_strategy->RemoteCellCommunication(d_deltaW1);

    d_hder_strategy->NearCellCommunication(d_hder_strategy->getResidual());

    if(d_hder_strategy->cur_proc==-1)std::cout<<"FInishing Remote cell communication\n";
    // 4.For all {SendCell_lv1} + {RecvCell_lv1}, using its updated ranking perform forward Sweep.
    d_hder_strategy->SolveBoundaryForwardSweep(d_diagOperator,W_scratch,d_deltaW1,ForwardSweep_Only,nForwardSweep_Only);
    if(d_hder_strategy->cur_proc==-1)std::cout<<"FInishing boundary forward Sweep\n";
}

void LUSGSStrategy::SolveBackwardSweep()
{
    // 5.For all {SendCell_lv1} + {RecvCell_lv1}, using its updated ranking perform backward Sweep.
    d_hder_strategy->SolveBoundaryBackwardSweep(d_diagOperator,W_scratch,d_deltaW1,d_deltaW,ForwardSweep_Only,nForwardSweep_Only);
    if(d_hder_strategy->cur_proc==-1)std::cout<<"Finishing boundary backward Sweep\n";
    // 6.For all {localCell} - {SendCell_lv1}, using its updated ranking perform backward Sweep.
    d_hder_strategy->SolveLocalBackwardSweep(d_diagOperator,W_scratch,d_deltaW1,d_deltaW,ForwardSweep_exceptions,nForwardSweep_exceptions);
    if(d_hder_strategy->cur_proc==-1)std::cout<<"Finishing local backward Sweep\n";
}

void LUSGSStrategy::SolveForwardSweepSerial()
{
    d_hder_strategy->SolveForwardSweep(d_diagOperator,W_scratch,d_deltaW1);
}

void LUSGSStrategy::SolveBackwardSweepSerial()
{
    d_hder_strategy->SolveBackwardSweep(d_diagOperator,W_scratch,d_deltaW1,d_deltaW);

}

void LUSGSStrategy::UpdateConservativeForRecvCells()
{
    double Gamma = 1.4;
    double*** U = d_hder_strategy->getU();
    for(int i = 0;i<d_nmesh;i++)
    {
        for(auto iter = d_hder->recvBegin(i);iter!=d_hder->recvEnd(i);iter++)
        {
            int cellId = iter->d_localId;
            W_scratch[i][0][cellId] = U[i][0][cellId];
            W_scratch[i][2][cellId] = U[i][0][cellId]*U[i][2][cellId];
            W_scratch[i][3][cellId] = U[i][0][cellId]*U[i][3][cellId];
            W_scratch[i][1][cellId] = U[i][1][cellId]/(Gamma - 1) 
            + 0.5*U[i][0][cellId]*(powf64(U[i][2][cellId],2)+powf64(U[i][3][cellId],2));
        }
        // for(int k = 0;k<nForwardSweep_Only[i];k++)
        // {
        //     int cellId = ForwardSweep_Only[i][k];
        //     W_scratch[i][0][cellId] = U[i][0][cellId];
        //     W_scratch[i][2][cellId] = U[i][0][cellId]*U[i][2][cellId];
        //     W_scratch[i][3][cellId] = U[i][0][cellId]*U[i][3][cellId];
        //     W_scratch[i][1][cellId] = U[i][1][cellId]/(Gamma - 1) 
        //     + 0.5*U[i][0][cellId]*(powf64(U[i][2][cellId],2)+powf64(U[i][3][cellId],2));
        // }
    }
}