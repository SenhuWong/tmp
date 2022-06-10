#pragma once

#include "TimeIntegrator.h"
#include <time.h>
class RungeKuttaStrategy : public TimeIntegrator
{
private:
    int d_nstage = -1;
    double d_CFL_number = 0.0;
    double d_stage_efficient[10];
    double Residual = 0;
    // Members that used to advance the solution

    double **W_scratch = NULL;
    double **W_now = NULL;

public:
    void initializeData()override
    {
        d_hder_strategy->initializeData();
    }
    void writeCell()
    {
        d_hder_strategy->AllCellCommunication(d_hder_strategy->getU());
        //d_hder_strategy->writeCellData("AfterCommFinal",d_hder_strategy->cur_proc, 0, d_hder_strategy->getU());
    }

    void singleStepSerial(int curStep)
    {
        preprocessUpdate();
        for (int i = 0;i<d_nstage;i++)
        {
            UpdateSerial(i);

        }
        postprocessUpdate();
    }

    void UpdateSerial(int iStage)
    {
        double** stableDt = d_hder_strategy->getStableTime();
        const double Gamma = 1.4;
        d_hder_strategy->preprocessAdvanceSerial(iStage);
        if(iStage==0)
        {
            d_hder_strategy->SolveTime(d_CFL_number);
        }
        double ** Residual = d_hder_strategy->getResidual();
        double **U = d_hder_strategy->getU();
        for (int l = 0; l < d_nstage; l++)
        {
            for (int i = 0; i < d_nmesh; i++)
            {
                auto &curBlk = d_hder->blk2D[i];

                for (int k = 0; k < d_hder->nCells(i); k++)
                {
                    auto &curCell = curBlk.d_localCells[k];

                    for (int j = 0; j < d_NEQU; j++)
                    {
                        // W_scratch is in conservative form
                        W_now[i][j+d_NEQU*k] = W_scratch[i][j+d_NEQU*k] + d_stage_efficient[iStage] * stableDt[i][k] * (-Residual[i][j+d_NEQU*k]) / curCell.volume();
                    }
                    // W[i][0][k] is the same
                }
            }
        }
        d_hder_strategy->initializePrimData(W_now);
        //d_hder_strategy->postprocessAdvance(iStage); // Do nothing

    }

    void singleStep(int curStep)
    {
        // clock_t start = clock();
        // These are thing should be done when ir==0
        preprocessUpdate();
        // clock_t end = clock();
        // std::cout<<"Time cost for preprocess is "<<end-start<<'\n';
        for (int i = 0; i < d_nstage; i++)
        {
            Update(i);
        }
        postprocessUpdate();
    }
    bool converged()
    {
        double **U = d_hder_strategy->getU();
        for (int i = 0; i < d_nmesh; i++)
        {
            Residual = 0;
            auto &curBlk = d_hder->blk2D[i];
            for (int j = 0; j < curBlk.nCells(); j++)
            {
                Residual += std::abs<double>(U[i][0+d_NEQU*j] - W_scratch[i][0+d_NEQU*j]);

            }
            if (i == 0)
            {
                std::cout << Residual << '\n';
            }
        }
        return true;
    }
    RungeKuttaStrategy(int nstage)
        : TimeIntegrator()
    {
        d_nstage = nstage;

    }

    void initializeSubStrategy()override
    {
        W_scratch = new double *[d_nmesh];
        W_now = new double*[d_nmesh];
        for (int i = 0; i < d_nmesh; i++)
        {
            W_scratch[i] = new double [d_NEQU*d_hder->nCells(i)];
            W_now[i] = new double[d_NEQU*d_hder->nCells(i)];
        }
        switch (d_nstage)
        {
        case 5:
            d_CFL_number = 0.8;
            d_stage_efficient[0] = 1.0 / 4;
            d_stage_efficient[1] = 1.0 / 6;
            d_stage_efficient[2] = 3.0 / 8;
            d_stage_efficient[3] = 1.0 / 2;
            d_stage_efficient[4] = 1.0;

            break;
        default:
            throw std::runtime_error("RungeKUtta Stage count out of range\n");
            break;
        }

    }
    
    RungeKuttaStrategy(UnstructTopologyHolder *hder, TopologyHolderStrategy *hder_strategy, int nstage)
        : TimeIntegrator(hder, hder_strategy)
    {
        W_scratch = new double *[d_nmesh];
        W_now = new double*[d_nmesh];
        for (int i = 0; i < d_nmesh; i++)
        {
            W_scratch[i] = new double [d_NEQU*hder->nCells(i)];
            W_now[i] = new double[d_NEQU*hder->nCells(i)];
            
        }
        d_nstage = nstage;
        switch (d_nstage)
        {
        case 5:
            d_CFL_number = 0.8;
            d_stage_efficient[0] = 1.0 / 4;
            d_stage_efficient[1] = 1.0 / 6;
            d_stage_efficient[2] = 3.0 / 8;
            d_stage_efficient[3] = 1.0 / 2;
            d_stage_efficient[4] = 1.0;

            break;
        default:
            throw std::runtime_error("RungeKUtta Stage count out of range\n");
            break;
        }
    }
    ~RungeKuttaStrategy()
    {
        for (int i = 0; i < d_nmesh; i++)
        {
            delete[] W_scratch[i];
            delete[] W_now[i];
        }
        delete[] W_scratch;
        delete[] W_now;
    }

public:
    void preprocessUpdate()
    {
        d_hder_strategy->initializeConsData(W_scratch);
    }
    void postprocessUpdate()
    {
        // Do A Cleaning of BufferCells
        // d_hder_strategy->cleaning(W_scratch);
    }
    void Update(int iStage)
    {
        double** stableDt = d_hder_strategy->getStableTime();
        const double Gamma = 1.4;
        // Before advance, call Euler's preprocessAdvance to recompute the modified
        d_hder_strategy->preprocessAdvance(iStage);
        if(iStage==0)
        {
            d_hder_strategy->SolveTime(d_CFL_number);
        }
        double ** Residual = d_hder_strategy->getResidual();
        double ** U = d_hder_strategy->getU();
        for (int l = 0; l < d_nstage; l++)
        {
            for (int i = 0; i < d_nmesh; i++)
            {
                auto &curBlk = d_hder->blk2D[i];

                for (int k = 0; k < d_hder->nCells(i); k++)
                {
                    auto &curCell = curBlk.d_localCells[k];

                    for (int j = 0; j < d_NEQU; j++)
                    {
                        // W_scratch is in conservative form
                        W_now[i][j+d_NEQU*k] = W_scratch[i][j+d_NEQU*k] + d_stage_efficient[iStage] * stableDt[i][k] * (-Residual[i][j+d_NEQU*k]) / curCell.volume();
                    }
                }
            }
        }
        d_hder_strategy->initializePrimData(W_now);
        // d_hder_strategy->postprocessAdvance(iStage); // Do nothing
    }
};