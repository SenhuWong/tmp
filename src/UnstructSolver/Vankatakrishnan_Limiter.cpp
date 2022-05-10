#include "Vankatakrishnan_Limiter.h"
#include "Euler2D.h"
#include <fstream>

Vankatakrishnan_Limiter::Vankatakrishnan_Limiter(UnstructTopologyHolder *hder, Euler2D *hder_strategy)
    : LimiterStrategy(hder, hder_strategy)
{
    Umax = new double **[d_nmesh];
    Umin = new double **[d_nmesh];
    for (int i = 0; i < d_nmesh; i++)
    {
        Umax[i] = new double *[d_nequ];
        Umin[i] = new double *[d_nequ];
        for (int j = 0; j < d_nequ; j++)
        {
            Umax[i][j] = new double[d_hder->nCells(i)];
            Umin[i][j] = new double[d_hder->nCells(i)];
        }
    }
}

void Vankatakrishnan_Limiter::computeLimiter() // cell_ind starts from 0
{
    double ***UL = d_hder_strategy->getUL();
    double ***UR = d_hder_strategy->getUR();
    // if (cur_proc == 0)
    //     std::cout << "Compute Limiter being called\n";
    // std::ofstream recorder;
    // std::string localFilename = "VankaRecorder" + std::to_string(cur_proc);
    // recorder.open(localFilename);
    //  FIrst get Umax and Umin
    double ***U = d_hder_strategy->getU();
    double ***Ue = d_hder_strategy->getUEdge();
    // double ***W = d_hder_strategy->getPrimVariables();
    // double ***Wi = d_hder_strategy->getWi();
    for (int i = 0; i < d_nmesh; i++)
    {
        auto &curBlk = d_hder->blk2D[i];
        for (int k = 0; k < d_hder->nEdges(i); k++)
        {
            auto &curEdge = curBlk.d_localEdges[k];
            int lc = curEdge.lCInd();
            int rc = curEdge.rCInd();
            if (rc >= 0)
            {
                for (int j = 0; j < d_nequ; j++)
                {
                    Umax[i][j][lc] = std::max(Umax[i][j][lc], U[i][j][rc]);
                    Umin[i][j][lc] = std::min(Umin[i][j][lc], U[i][j][rc]);
                    Umax[i][j][rc] = std::max(Umax[i][j][rc], U[i][j][lc]);
                    Umin[i][j][rc] = std::min(Umin[i][j][rc], U[i][j][lc]);
                }
            }
        }
    }
    // Then get delta1(max,min),which are stored in Umax and Umin
    for (int i = 0; i < d_nmesh; i++)
    {
        auto &curBlk = d_hder->blk2D[i];
        for (int j = 0; j < d_nequ; j++)
        {
            for (int k = 0; k < d_hder->nCells(i); k++)
            {
                Umax[i][j][k] = Umax[i][j][k] - U[i][j][k];
                Umin[i][j][k] = Umin[i][j][k] - U[i][j][k];
            }
        }
    }
    GeomElements::vector3d<2, double> ***grad = d_hder_strategy->getGradient();

    for (int i = 0; i < d_nmesh; i++)
    {
        auto &curBlk = d_hder->blk2D[i];
        for (int k = 0; k < curBlk.d_nEs; k++)
        {
            auto &curEdge = curBlk.d_localEdges[k];
            int lc = curEdge.lCInd();
            int rc = curEdge.rCInd();
            //recorder << "Edge " << k << "'s left is " << lc << " right is " << rc << '\n';
            //recorder << "curEdge's center is " << curEdge.center()[0] << "," << curEdge.center()[1] << '\n';
            if (lc >= 0)
            {
                GeomElements::vector3d<2, double> leftVec = curEdge.center() - curBlk.d_localCells[lc].center();
                //recorder << "Left cell's center is " << curBlk.d_localCells[lc].center()[0] << "," << curBlk.d_localCells[lc].center()[1] << '\n';
                //recorder << "Left vec is " << leftVec[0] << "," << leftVec[1] << '\n';
                double epsilonSquare = powf64(LIMITER_K, d_dim) * d_hder->CellVolume(i, lc);
                double result_limit;
                for (int j = 0; j < d_nequ; j++)
                {

                    double delta2 = grad[i][j][lc].dot_product(leftVec);
                    // recorder << "\t " << j << "th grad and product is " << grad[i][j][lc][0] << "," << grad[i][j][lc][1] << "\t" << delta2 << '\n';

                    UL[i][j][k] = delta2; // The left increment.
                    double delta1_max = Umax[i][j][lc];
                    double delta1_min = Umin[i][j][lc];
                    // recorder << "\tUmax and Umin are " << delta1_max << "," << delta1_min << '\n';
                    if (delta2 > 0)
                    {
                        result_limit = (delta1_max * delta1_max + epsilonSquare + 2 * delta2 * delta1_max) / (delta1_max * delta1_max + 2 * delta2 * delta2 + delta1_max * delta2 + epsilonSquare);
                    }
                    else if (delta2 < 0)
                    {
                        result_limit = (delta1_min * delta1_min + epsilonSquare + 2 * delta2 * delta1_min) / (delta1_min * delta1_min + 2 * delta2 * delta2 + delta1_min * delta2 + epsilonSquare);
                    }
                    else
                    {
                        result_limit = 1;
                    }
                    // recorder << "\tResult limit is " << result_limit << '\n';
                    // mesh_var_cell_limit[i][j][lc] = std::min(mesh_var_cell_limit[i][j][lc], result_limit);
                }
            }
            else
            {
                throw std::runtime_error("Left cell can't be null\n");
            }
            if (rc >= 0)
            {
                GeomElements::vector3d<2, double> rightVec = curEdge.center() - curBlk.d_localCells[rc].center();
                // recorder << "Right cell's center is " << curBlk.d_localCells[rc].center()[0] << "," << curBlk.d_localCells[rc].center()[1] << '\n';
                // recorder << "Right vec is " << rightVec[0] << "," << rightVec[1] << '\n';
                double epsilonSquare = powf64(LIMITER_K, d_dim) * d_hder->CellVolume(i, rc);
                double result_limit;
                for (int j = 0; j < d_nequ; j++)
                {
                    double delta2 = grad[i][j][rc].dot_product(rightVec);
                    // recorder << "\t " << j << "th grad and product is " << grad[i][j][rc][0] << "," << grad[i][j][rc][1] << "\t" << delta2 << '\n';
                    UR[i][j][k] = delta2;

                    double delta1_max = Umax[i][j][rc];
                    double delta1_min = Umin[i][j][rc];
                    // recorder << "\tUmax and Umin are " << delta1_max << "," << delta1_min << '\n';
                    if (delta2 > 0)
                    {
                        result_limit = (delta1_max * delta1_max + epsilonSquare + 2 * delta2 * delta1_max) / (delta1_max * delta1_max + 2 * delta2 * delta2 + delta1_max * delta2 + epsilonSquare);
                    }
                    else if (delta2 < 0)
                    {
                        result_limit = (delta1_min * delta1_min + epsilonSquare + 2 * delta2 * delta1_min) / (delta1_min * delta1_min + 2 * delta2 * delta2 + delta1_min * delta2 + epsilonSquare);
                    }
                    else
                    {
                        result_limit = 1;
                    }
                    // recorder << "\tResult limit is " << result_limit << '\n';
                    mesh_var_cell_limit[i][j][rc] = std::min(mesh_var_cell_limit[i][j][rc], result_limit);
                }
            }
        }
    }
    // recorder.close();
}

double Vankatakrishnan_Limiter::getLimiter(int meshInd, int equInd, int cellInd)
{
    return mesh_var_cell_limit[meshInd][equInd][cellInd];
}

// Set Umax Umin and limiter to initial value.

void Vankatakrishnan_Limiter::init()
{
    double ***U = d_hder_strategy->getU();
    for (int i = 0; i < d_nmesh; i++)
    {
        for (int j = 0; j < d_nequ; j++)
        {

            for (int k = 0; k < d_hder->nCells(i); k++)
            {
                Umax[i][j][k] = U[i][j][k];
                Umin[i][j][k] = U[i][j][k];
                mesh_var_cell_limit[i][j][k] = 1.0;
            }
        }
    }
}


void Vankatakrishnan_Limiter::write_Lmitera()
{
    for (int i = 0; i < d_nmesh; i++)
    {
        std::ofstream fout;
        fout.open("Limiter" + std::to_string(i) + std::to_string(cur_proc));
        auto &curBlk = d_hder->blk2D[i];
        fout << "TITLE =\"Tioga output\"\n";
        fout << "VARIABLES=\"X\",\"Y\",\"Limiter0\",\"Limiter1\",\"Limiter2\",\"Limiter3\"\n";
        fout << "ZONE T=\"VOL_MIXED\",N=" << curBlk.d_nPs << " E=" << curBlk.d_nCs << " ET = QUADRILATERAL, F = FEBLOCK\n";
        fout << "VARLOCATION =  (1=NODAL,2=NODAL,3=CELLCENTERED,4=CELLCENTERED,5=CELLCENTERED,6=CELLCENTERED)\n";
        for (int j = 0; j < d_dim; j++)
        {
            for (int k = 0; k < curBlk.d_nPs; k++)
            {
                fout << curBlk.d_localPoints[k][j] << '\n';
            }
        }
        for (int j = 0; j < d_nequ; j++)
        {
            for (int k = 0; k < curBlk.d_nCs; k++)
            {
                fout << (mesh_var_cell_limit[i][j][k] > 1) << '\n';
            }
        }
        for (int k = 0; k < curBlk.d_nCs; k++)
        {
            auto &curCell = curBlk.d_localCells[k];
            for (int j = 0; j < curCell.size(); j++)
            {
                fout << curCell.pointInd(j) + 1 << " ";
            }
            for (int j = curCell.size(); j < 4; j++)
            {
                fout << curCell.pointInd(curCell.size() - 1) + 1 << " ";
            }
            fout << '\n';
        }
        fout.close();
    }
}