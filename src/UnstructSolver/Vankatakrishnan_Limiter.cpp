#include "Vankatakrishnan_Limiter.h"
#include "UnstructIntegrator.h"
#include "toolBox/vector3d.h"
#include <fstream>

Vankatakrishnan_Limiter::Vankatakrishnan_Limiter(UnstructTopologyHolder *hder, TopologyHolderStrategy *hder_strategy)
    : LimiterStrategy(hder, hder_strategy)
{
    mesh_var_cell_limit = new double*[d_nmesh];
    Umax = new double *[d_nmesh];
    Umin = new double *[d_nmesh];
    for (int i = 0; i < d_nmesh; i++)
    {
        mesh_var_cell_limit[i] = new double [d_NEQU*d_hder->nCells(i)];
        Umax[i] = new double [d_NEQU*d_hder->nCells(i)];
        Umin[i] = new double [d_NEQU*d_hder->nCells(i)];
    }
}

Vankatakrishnan_Limiter::~Vankatakrishnan_Limiter()
{
    for(int i = 0;i < d_nmesh;i++)
    {
        delete[] mesh_var_cell_limit[i];
        delete[] Umax[i];
        delete[] Umin[i];
    }
    delete[] mesh_var_cell_limit;
    delete[] Umax;
    delete[] Umin;
}

void Vankatakrishnan_Limiter::computeLimiter() // cell_ind starts from 0
{
    double **UL = d_hder_strategy->getUL();
    double **UR = d_hder_strategy->getUR();
    //  FIrst get Umax and Umin
    double **U = d_hder_strategy->getU();

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
                for (int j = 0; j < d_NEQU; j++)
                {
                    Umax[i][j+d_NEQU*lc] = std::max(Umax[i][j+d_NEQU*lc], U[i][j+d_NEQU*rc]);
                    Umin[i][j+d_NEQU*lc] = std::min(Umin[i][j+d_NEQU*lc], U[i][j+d_NEQU*rc]);
                    Umax[i][j+d_NEQU*rc] = std::max(Umax[i][j+d_NEQU*rc], U[i][j+d_NEQU*lc]);
                    Umin[i][j+d_NEQU*rc] = std::min(Umin[i][j+d_NEQU*rc], U[i][j+d_NEQU*lc]);
                }
            }
        }
    }
    // Then get delta1(max,min),which are stored in Umax and Umin
    for (int i = 0; i < d_nmesh; i++)
    {
        auto &curBlk = d_hder->blk2D[i];
        for (int j = 0; j < d_NEQU; j++)
        {
            for (int k = 0; k < d_hder->nCells(i); k++)
            {
                Umax[i][j+d_NEQU*k] = Umax[i][j+d_NEQU*k] - U[i][j+d_NEQU*k];
                Umin[i][j+d_NEQU*k] = Umin[i][j+d_NEQU*k] - U[i][j+d_NEQU*k];
            }
        }
    }
    GeomElements::vector3d<2, double> **grad =(GeomElements::vector3d<2, double> **)(d_hder_strategy->getGradient());
    for (int i = 0; i < d_nmesh; i++)
    {
        auto &curBlk = d_hder->blk2D[i];
        for (int k = 0; k < curBlk.nEdges(); k++)
        {
            auto &curEdge = curBlk.d_localEdges[k];
            int lc = curEdge.lCInd();
            int rc = curEdge.rCInd();
            if (lc >= 0)
            {
                GeomElements::vector3d<2, double> leftVec = curEdge.center() - curBlk.d_localCells[lc].center();
                double epsilonSquare = powf64(LIMITER_K, d_dim) * d_hder->CellVolume(i, lc);
                double result_limit;
                for (int j = 0; j < d_NEQU; j++)
                {
                    double delta2 = grad[i][j+d_NEQU*lc].dot_product(leftVec);
                    UL[i][j+d_NEQU*k] = delta2; // The left increment.
                    double delta1_max = Umax[i][j+d_NEQU*lc];
                    double delta1_min = Umin[i][j+d_NEQU*lc];
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
                    mesh_var_cell_limit[i][j+d_NEQU*lc] = std::min(mesh_var_cell_limit[i][j+d_NEQU*lc], result_limit);
                }
            }
            else
            {
                throw std::runtime_error("Left cell can't be null\n");
            }
            if (rc >= 0)
            {
                GeomElements::vector3d<2, double> rightVec = curEdge.center() - curBlk.d_localCells[rc].center();
                double epsilonSquare = powf64(LIMITER_K, d_dim) * d_hder->CellVolume(i, rc);
                double result_limit;
                for (int j = 0; j < d_NEQU; j++)
                {
                    double delta2 = grad[i][j+d_NEQU*rc].dot_product(rightVec);
                    UR[i][j+d_NEQU*k] = delta2;
                    double delta1_max = Umax[i][j+d_NEQU*rc];
                    double delta1_min = Umin[i][j+d_NEQU*rc];
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
                    mesh_var_cell_limit[i][j+d_NEQU*rc] = std::min(mesh_var_cell_limit[i][j+d_NEQU*rc], result_limit);
                }
            }
        }
    }
    // recorder.close();
}

double Vankatakrishnan_Limiter::getLimiter(int meshInd, int equInd, int cellInd)
{
    return mesh_var_cell_limit[meshInd][equInd+d_NEQU*cellInd];
}

// Set Umax Umin and limiter to initial value.

void Vankatakrishnan_Limiter::init()
{
    double **U = d_hder_strategy->getU();
    for (int i = 0; i < d_nmesh; i++)
    {
        for (int j = 0; j < d_NEQU; j++)
        {
            for (int k = 0; k < d_hder->nCells(i); k++)
            {
                Umax[i][j+d_NEQU*k] = U[i][j+d_NEQU*k];
                Umin[i][j+d_NEQU*k] = U[i][j+d_NEQU*k];
                mesh_var_cell_limit[i][j+d_NEQU*k] = 1.0;
            }
        }
    }
}


// void Vankatakrishnan_Limiter::write_Lmitera()
// {
//     for (int i = 0; i < d_nmesh; i++)
//     {
//         std::ofstream fout;
//         fout.open("Limiter" + std::to_string(i) + std::to_string(cur_proc));
//         auto &curBlk = d_hder->blk2D[i];
//         fout << "TITLE =\"Tioga output\"\n";
//         fout << "VARIABLES=\"X\",\"Y\",\"Limiter0\",\"Limiter1\",\"Limiter2\",\"Limiter3\"\n";
//         fout << "ZONE T=\"VOL_MIXED\",N=" << curBlk.d_nPs << " E=" << curBlk.d_nCs << " ET = QUADRILATERAL, F = FEBLOCK\n";
//         fout << "VARLOCATION =  (1=NODAL,2=NODAL,3=CELLCENTERED,4=CELLCENTERED,5=CELLCENTERED,6=CELLCENTERED)\n";
//         for (int j = 0; j < d_dim; j++)
//         {
//             for (int k = 0; k < curBlk.d_nPs; k++)
//             {
//                 fout << curBlk.d_localPoints[k][j] << '\n';
//             }
//         }
//         for (int j = 0; j < d_NEQU; j++)
//         {
//             for (int k = 0; k < curBlk.d_nCs; k++)
//             {
//                 fout << (mesh_var_cell_limit[i][j][k] > 1) << '\n';
//             }
//         }
//         for (int k = 0; k < curBlk.d_nCs; k++)
//         {
//             auto &curCell = curBlk.d_localCells[k];
//             for (int j = 0; j < curCell.size(); j++)
//             {
//                 fout << curCell.pointInd(j) + 1 << " ";
//             }
//             for (int j = curCell.size(); j < 4; j++)
//             {
//                 fout << curCell.pointInd(curCell.size() - 1) + 1 << " ";
//             }
//             fout << '\n';
//         }
//         fout.close();
//     }
// }