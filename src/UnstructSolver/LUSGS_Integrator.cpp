#include "LUSGS_Integrator.h"
#include "UnstructIntegrator.h"
#include "TopologyHolderStrategy.h"
#include <algorithm>
#include <iostream>
#include <mpi/mpi.h>
#include "toolBox/edge3d_int.h"

LUSGSIntegrator::LUSGSIntegrator()
{

}


void LUSGSIntegrator::initializeSubStrategy()
{
    initializeParallelSweep();
    d_diagOperator = new double*[d_nmesh];
    d_deltaW = new double*[d_nmesh];
    d_deltaW1 = new double*[d_nmesh];
    W_scratch = new double*[d_nmesh];

    for(int i = 0 ;i<d_nmesh;i++)
    {
        d_diagOperator[i] = new double[d_NEQU*d_hder->nCells(i)];
        d_deltaW[i] = new double[d_NEQU*d_hder->nCells(i)];
        d_deltaW1[i] = new double[d_NEQU*d_hder->nCells(i)];
        W_scratch[i] = new double[d_NEQU*d_hder->nCells(i)];
    }

}

void LUSGSIntegrator::initializeParallelSweep()
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
}


LUSGSIntegrator::LUSGSIntegrator(UnstructTopologyHolder *hder, TopologyHolderStrategy* hder_strategy)
    :TimeIntegrator(hder,hder_strategy)
{
    initializeParallelSweep();

    d_diagOperator = new double*[d_nmesh];
    d_deltaW = new double*[d_nmesh];
    d_deltaW1 = new double*[d_nmesh];
    W_scratch = new double*[d_nmesh];

    for(int i = 0 ;i<d_nmesh;i++)
    {
        d_diagOperator[i] = new double[d_NEQU*d_hder->nCells(i)];
        d_deltaW[i] = new double[d_NEQU*d_hder->nCells(i)];
        d_deltaW1[i] = new double[d_NEQU*d_hder->nCells(i)];
        W_scratch[i] = new double[d_NEQU*d_hder->nCells(i)];
    }
}
void LUSGSIntegrator::initializeData()
{
    d_hder_strategy->initializeData();
    d_hder_strategy->initializeConsData(W_scratch);
}

void LUSGSIntegrator::preprocessUpdate()
{

}

void LUSGSIntegrator::singleStep(int curStage)
{
    Update();
    // postprocessUpdate();
}

void LUSGSIntegrator::singleStepSerial(int curStep)
{
    UpdateSerial();
    //postprocessUpdate();
    
}

void LUSGSIntegrator::postprocessUpdate()
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


bool LUSGSIntegrator::converged()
{

    

}

void LUSGSIntegrator::preprocessLUSGS()
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

void LUSGSIntegrator::Update()
{
    d_hder_strategy->preprocessAdvance(0);
    
    d_hder_strategy->SolveTime(d_CFL_number);
    
    d_hder_strategy->UpdateConservativeForRecvCells(W_scratch);

    d_hder_strategy->SolveDiag(d_diagOperator);
    d_hder_strategy->postprocessAdvance(0);
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

void LUSGSIntegrator::UpdateSerial()
{
    d_hder_strategy->preprocessAdvanceSerial(0);
    d_hder_strategy->SolveTime(d_CFL_number);

    d_hder_strategy->SolveDiag(d_diagOperator);

    SolveForwardSweepSerial();

    SolveBackwardSweepSerial();

    UpdatePrimitiveVariable();
}

void LUSGSIntegrator::UpdatePrimitiveVariable()
{
    for(int i = 0;i<d_nmesh;i++)
    {
        for(int k = 0;k<d_hder->nCells(i);k++)
        {
            for(int j = 0;j<d_NEQU;j++)
            {
                W_scratch[i][j+d_NEQU*k] = W_scratch[i][j+d_NEQU*k] + d_deltaW[i][j+d_NEQU*k];
            }
        }
    }
    d_hder_strategy->UpdatePrimitiveVariables(W_scratch);
}

void LUSGSIntegrator::postprocessLUSGS()
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


void LUSGSIntegrator::SolveForwardSweep()
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

void LUSGSIntegrator::SolveBackwardSweep()
{
    // 5.For all {SendCell_lv1} + {RecvCell_lv1}, using its updated ranking perform backward Sweep.
    SolveBoundaryBackwardSweep(d_diagOperator,W_scratch,d_deltaW1,d_deltaW,ForwardSweep_Only,nForwardSweep_Only);
    if(d_hder_strategy->cur_proc==-1)std::cout<<"Finishing boundary backward Sweep\n";
    // 6.For all {localCell} - {SendCell_lv1}, using its updated ranking perform backward Sweep.
    SolveLocalBackwardSweep(d_diagOperator,W_scratch,d_deltaW1,d_deltaW,ForwardSweep_exceptions,nForwardSweep_exceptions);
    if(d_hder_strategy->cur_proc==-1)std::cout<<"Finishing local backward Sweep\n";
}

void LUSGSIntegrator::SolveForwardSweepSerial()
{
    SolveForwardSweep(d_diagOperator,W_scratch,d_deltaW1);
}

void LUSGSIntegrator::SolveBackwardSweepSerial()
{
    SolveBackwardSweep(d_diagOperator,W_scratch,d_deltaW1,d_deltaW);

}




//The skip point in here are level 1 sendcells and level 1,2 recvcells.
void LUSGSIntegrator::SolveLocalForwardSweep(double** diag,
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
                    d_hder_strategy->SolveDF(i,cellId,curCell.edgeInd(j),W_scratch,deltaW1,dF,ForB,cellIdOffset);
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
void LUSGSIntegrator::SolveBoundaryForwardSweep(double** diag,
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
                d_hder_strategy->SolveDF(i,cellId,curCell.edgeInd(j),W_scratch,deltaW1,dF,ForB,cellIdOffset);
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
void LUSGSIntegrator::SolveBoundaryBackwardSweep(double** diag,
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

        for(int k = nskip[i]-1;k>=0;k--)
        {
            
            int cellId = skip_point[i][k];
            auto& curCell = curBlk.d_localCells[cellId];
            for(int j = 0;j<d_NEQU;j++)
            {
                dFi[j] = 0;
            }
            for(int j = 0;j<curCell.edge_size();j++)
            {
                d_hder_strategy->SolveDF(i,cellId,curCell.edgeInd(j),W_scratch,deltaW,dF,ForB,cellIdOffset);
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
void LUSGSIntegrator::SolveLocalBackwardSweep(double** diag,
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
        {
            for(int l = skip_point[i][k+1]-1;l>skip_point[i][k];l--)
            {
                int cellId = l;
                auto& curCell = curBlk.d_localCells[cellId];
                for(int j =0;j<d_NEQU;j++)
                {
                    dFi[j] = 0;
                }
                for(int j = 0;j<curCell.edge_size();j++)
                {
                    d_hder_strategy->SolveDF(i,cellId,curCell.edgeInd(j),W_scratch,deltaW,dF,ForB,cellIdOffset);
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


void LUSGSIntegrator::SolveForwardSweep(double** diag,
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
                d_hder_strategy->SolveDF(i,k,curCell.edgeInd(j),W_scratch,deltaW1,dF,ForB,0);
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


void LUSGSIntegrator::SolveBackwardSweep(double** diag,
                        double** W_scratch,double** deltaW1,double** deltaW)
{
    int ForB = 0;
    double dFi[d_NEQU]={0};
    double dF[d_NEQU]={0};
    for(int i = 0;i<d_nmesh;i++)
    {
        auto& curBlk = d_hder->blk2D[i];
        for(int k = d_hder->nCells(i)-1;k>=0;k--)
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
                d_hder_strategy->SolveDF(i,k,curCell.edgeInd(j),W_scratch,deltaW,dF,ForB,0);
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