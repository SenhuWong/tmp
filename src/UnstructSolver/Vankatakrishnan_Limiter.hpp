
#pragma once
#define LIMITER_K 5
// This should be a subclass for limiter .

// Looks like I have to write limiter one for unstruct and another for struct.
// Or if all limiter want to know max and min, we could implement the filling
// of max and min in Strategy.
// And after that is filled ,Wmax and Wmin could be directly used to store
// delta1,max and delta1,min.
// For every edge that we want to compute limiter, we need to pass in the delta2
#include "LimiterStrategy.h"
#include "UnstructIntegrator.h"
#include "toolBox/vector3d.h"
template<int ndim>
class Vankatakrishnan_Limiter : public LimiterStrategy
{
private:
    double ** Umax=NULL;
    double ** Umin=NULL;
    double ** mesh_var_cell_limit = NULL;
    
public:
    Vankatakrishnan_Limiter(UnstructTopologyHolder *hder, TopologyHolderStrategy *hder_strategy);

    ~Vankatakrishnan_Limiter();

    void computeLimiter(double** UL, double** UR); // cell_ind starts from 0
    
    double getLimiter(int meshInd, int equInd, int cellInd);

    // void write_Lmitera();

    void computeLeftRight();

    void init();

};


template<int ndim>
Vankatakrishnan_Limiter<ndim>::Vankatakrishnan_Limiter(UnstructTopologyHolder *hder, TopologyHolderStrategy *hder_strategy)
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
template<int ndim>
Vankatakrishnan_Limiter<ndim>::~Vankatakrishnan_Limiter()
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
template<int ndim>
void Vankatakrishnan_Limiter<ndim>::computeLimiter(double **UL,double **UR) // cell_ind starts from 0
{
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
    GeomElements::vector3d<ndim, double> **grad =(GeomElements::vector3d<ndim, double> **)(d_hder_strategy->getGradientPrimitive());
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
                GeomElements::vector3d<ndim, double> leftVec = curEdge.center() - curBlk.d_localCells[lc].center();
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
                GeomElements::vector3d<ndim, double> rightVec = curEdge.center() - curBlk.d_localCells[rc].center();
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

template<int ndim>
double Vankatakrishnan_Limiter<ndim>::getLimiter(int meshInd, int equInd, int cellInd)
{
    return mesh_var_cell_limit[meshInd][equInd+d_NEQU*cellInd];
}

// Set Umax Umin and limiter to initial value.
template<int ndim>
void Vankatakrishnan_Limiter<ndim>::init()
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
