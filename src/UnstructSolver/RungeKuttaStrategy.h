#pragma once

#include "TimeStrategy.h"
struct FlowParamHolder
{
    double density;
    double pressure;
    double soundSpeed;
    double velocity[3];
};
void computeFlowParams(double ***W, int mesh_ind, int cell_ind, int dim, FlowParamHolder *hder);

class RungeKuttaStrategy : public TimeStrategy
{
private:
    int d_nstage = -1;
    double d_CFL_number = 0.0;
    double d_stage_efficient[10];
    double Residual = 0;
    // Members that used to advance the solution

    //double ***VolumeProjection = NULL;
    double **stableDt = NULL;
    double ***W_scratch = NULL;

    // For debug
    int cur_proc = -1;
    int num_proc = -1;

public:
    void setComm(int cur, int num)
    {
        cur_proc = cur;
        num_proc = num;
    }
    void initialize()
    {
        d_hder_strategy->initializeData();
    }
    void writeCell()
    {
        d_hder_strategy->withinBlockCommunication();
        d_hder_strategy->writeCellData("AfterCommFinal", cur_proc, 0, d_hder_strategy->getU());
        //d_hder_strategy->updateEdgeConservatives_butWithOutput();
    }

    void outPut(int iStep)
    {
        //d_hder_strategy->writeTiogaFormat("Out"+std::to_string(iStep),0,0,d_hder_strategy->getU())
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
        const double Gamma = 1.4;
        // Before advance, call Euler's preprocessAdvance to recompute the modified
        // d_hder_strategy->writeCellData("BeforePreprocessAdvance" + std::to_string(iStage) + "stage", d_hder_strategy->cur_proc, 0, d_hder_strategy->getPrimVariables());
        d_hder_strategy->preprocessAdvanceSerial(iStage);
        // d_hder_strategy->writeCellData("AfterPreprocessAdvance" + std::to_string(iStage) + "stage", d_hder_strategy->cur_proc, 0, d_hder_strategy->getPrimVariables());
        if(iStage==0)
        {
            d_hder_strategy->SolveTime(stableDt,d_CFL_number);
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
                        U[i][j+d_NEQU*k] = W_scratch[i][j][k] + d_stage_efficient[iStage] * stableDt[i][k] * (-Residual[i][j+d_NEQU*k]) / curCell.volume();
                    }
                    // W[i][0][k] is the same
                    double VMSquare = 0;
                    for (int l = 0; l < d_dim; l++)
                    {
                        U[i][l + 2+d_NEQU*k] = U[i][l + 2+d_NEQU*k] / U[i][0+d_NEQU*k];
                        VMSquare += U[i][l + 2+d_NEQU*k] * U[i][l + 2+d_NEQU*k];
                    }
                    U[i][1+d_NEQU*k] = (Gamma - 1) * (U[i][1+d_NEQU*k] - 0.5 * U[i][0+d_NEQU*k] * VMSquare);
                }
            }
        }
        // d_hder_strategy->writeCellData("AfterUpdate" + std::to_string(iStage), d_hder_strategy->cur_proc, 0, d_hder_strategy->getPrimVariables());
        d_hder_strategy->postprocessAdvance(); // Do nothing

    }

    void singleStep(int curStep)
    {
        // These are thing should be done when ir==0
        preprocessUpdate();
        for (int i = 0; i < d_nstage; i++)
        {
            Update(i);
        }
        postprocessUpdate();
        //  d_hder_strategy->writeCellData("AfterCommFinal", cur_proc, 0, d_hder_strategy->getPrimVariables());

        // if (curStep % 10 == 0)
        // {
        //     std::cout << curStep << '\n';
        //     converged();
        //     d_hder_strategy->withinBlockCommunication();
        //     d_hder_strategy->writeCellData("AfterCommFinal" + std::to_string(curStep) + "s", cur_proc, 0, d_hder_strategy->getPrimVariables());
        //     // d_hder_strategy->updateEdgeConservatives_butWithOutput();
        // }
    }
    bool converged()
    {
        double **U = d_hder_strategy->getU();
        for (int i = 0; i < d_nmesh; i++)
        {
            Residual = 0;
            auto &curBlk = d_hder->blk2D[i];
            for (int j = 0; j < curBlk.d_nCs; j++)
            {
                Residual += abs(U[i][0+d_NEQU*j] - W_scratch[i][0][j]);
                // std::cout<<W[i][0][j]<<" "<<W_scratch[i][0][j]<<'\n';
            }
            if (i == 0)
            {
                std::cout << Residual << '\n';
            }
        }
        return true;
    }
    RungeKuttaStrategy(UnstructTopologyHolder *hder, Euler2D *hder_strategy, int nstage)
        : TimeStrategy(hder, hder_strategy)
    {
        //VolumeProjection = new double **[d_nmesh];
        stableDt = new double *[d_nmesh];
        W_scratch = new double **[d_nmesh];
        for (int i = 0; i < d_nmesh; i++)
        {

            // VolumeProjection[i] = new double *[d_dim];
            stableDt[i] = new double[hder->nCells(i)];
            W_scratch[i] = new double *[d_NEQU];
            // for (int j = 0; j < d_dim; j++)
            // {
            //     VolumeProjection[i][j] = new double[hder->nCells(i)];
            // }
            for (int j = 0; j < d_NEQU; j++)
            {
                W_scratch[i][j] = new double[hder->nCells(i)];
            }
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
            for (int j = 0; j < d_NEQU; j++)
            {
                delete[] W_scratch[i][j];
            }
            delete[] stableDt[i];
            delete[] W_scratch[i];
        }
        delete[] stableDt;
        delete[] W_scratch;
    }

private:
    void preprocessUpdate()
    {
        const double Gamma = 1.4;
        // Fill W into W_scratch
        double **U = d_hder_strategy->getU();
        for (int i = 0; i < d_nmesh; i++)
        {
            for (int k = 0; k < d_hder->nCells(i); k++)
            {
                double VMSquare = 0;
                for (int l = 0; l < d_dim; l++)
                {
                    VMSquare += U[i][l + 2+d_NEQU*k] * U[i][l + 2+d_NEQU*k];
                    W_scratch[i][l + 2][k] = U[i][0+d_NEQU*k] * U[i][l + 2+d_NEQU*k];
                }
                W_scratch[i][0][k] = U[i][0+d_NEQU*k];
                W_scratch[i][1][k] = U[i][1+d_NEQU*k] / (Gamma - 1) + 0.5 * U[i][0+d_NEQU*k] * VMSquare;
            }
        }
        // Clear Projection(Spectrum) for this step
        //clearProjection();
        // Compute Projection(since this won't change in static grid)
        //computeProjection();
    }
    void postprocessUpdate()
    {
        // Do A Cleaning of BufferCells
        // d_hder_strategy->cleaning(W_scratch);
    }
    void Update(int iStage)
    {
        const double Gamma = 1.4;
        // Before advance, call Euler's preprocessAdvance to recompute the modified
        // d_hder_strategy->writeCellData("BeforePreprocessAdvance" + std::to_string(iStage) + "stage", d_hder_strategy->cur_proc, 0, d_hder_strategy->getPrimVariables());
        d_hder_strategy->preprocessAdvance(iStage);
        // d_hder_strategy->writeCellData("AfterPreprocessAdvance" + std::to_string(iStage) + "stage", d_hder_strategy->cur_proc, 0, d_hder_strategy->getPrimVariables());
        if(iStage==0)
        {
            d_hder_strategy->SolveTime(stableDt,d_CFL_number);
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
                        U[i][j+d_NEQU*k] = W_scratch[i][j][k] + d_stage_efficient[iStage] * stableDt[i][k] * (-Residual[i][j+d_NEQU*k]) / curCell.volume();
                    }
                    // W[i][0][k] is the same
                    double VMSquare = 0;
                    for (int l = 0; l < d_dim; l++)
                    {
                        U[i][l + 2+d_NEQU*k] = U[i][l + 2+d_NEQU*k] / U[i][0+d_NEQU*k];
                        VMSquare += U[i][l + 2+d_NEQU*k] * U[i][l + 2+d_NEQU*k];
                    }
                    U[i][1+d_NEQU*k] = (Gamma - 1) * (U[i][1+d_NEQU*k] - 0.5 * U[i][0+d_NEQU*k] * VMSquare);
                }
            }
        }
        // d_hder_strategy->writeCellData("AfterUpdate" + std::to_string(iStage), d_hder_strategy->cur_proc, 0, d_hder_strategy->getPrimVariables());
        d_hder_strategy->postprocessAdvance(); // Do nothing
    }
};