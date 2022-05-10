#include "LimiterStrategy.h"
#include "Euler2D.h"

LimiterStrategy::LimiterStrategy(UnstructTopologyHolder *hder, Euler2D *hder_strategy)
    : d_nmesh(hder->d_nmesh), d_hder(hder), d_hder_strategy(hder_strategy)
{
    d_dim = d_hder->d_dim;
    d_nequ = d_hder_strategy->d_NEQU;
    mesh_var_cell_limit = new double **[d_nmesh];
    for (int i = 0; i < d_nmesh; i++)
    {
        mesh_var_cell_limit[i] = new double *[d_nequ];
        for (int j = 0; j < d_nequ; j++)
        {
            mesh_var_cell_limit[i][j] = new double[d_hder->nCells(i)];
        }
    }
}
LimiterStrategy::~LimiterStrategy()
{
    if (mesh_var_cell_limit)
    {
        for (int i = 0; i < d_nmesh; i++)
        {
            for (int j = 0; j < d_nequ; j++)
            {
                delete[] mesh_var_cell_limit[i][j];
            }
            delete[] mesh_var_cell_limit[i];
        }
        delete[] mesh_var_cell_limit;
    }
}