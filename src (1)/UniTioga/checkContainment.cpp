#include"MeshBlock.h"
#include<iostream>
void computeNodalWeights(int dim, double xv[8][3], double* xp, double* frac, int nvert);
void MeshBlock::checkContainment(int* cellIndex, int adtElement, double* x2search)
{
    int i, id;
    int n;
    int nvert;
    int icell = elementList[adtElement];
    //allocate enough space for 3D, but still can be used for 2D; 
    double xv[8][3];
    double frac[8];
    if (ihigh == 0)
    {
        //locate cell
        int isum = 0;
        for (n = 0; n < ntypes; n++)
        {
            isum += nc[n];
            if (icell < isum)
            {
                i = icell - (isum - nc[n]);
                break;
            }
        }
        nvert = nv[n];
        for (int m = 0; m < nvert; m++)
        {
            id = d_dim * (vconn[n][nvert * i + m] - BASE);
            for (int j = 0; j < d_dim; j++)
            {
                xv[m][j] = x[id + j];
            }

        }
        
        computeNodalWeights(d_dim, xv, x2search, frac, nvert);
        cellIndex[0] = icell;
        cellIndex[1] = 0;
        for (int m = 0; m < nvert; m++)
        {
            if ((frac[m] + TOL) * (frac[m] - 1.0 - TOL) > 0)
            {
                //std::cout << "Exclusiong make\n";
                cellIndex[0] = -1;
                return;
            }
            if (fabs(frac[m]) < TOL && cellRes[icell] == BIGVALUE)
            {
                cellIndex[1] = 1;
                //TODO:When is this necessary?
                //Valid inclusion yet a small frac and unacceptable donor is the element
            }
        }
        //std::cout << "Inclusiong madeeeeeeee\n";
        //iblank_cell[icell] = -1;
        return;
    }
    else
    {
        std::cout << "Error at MeshBlock::checkContainment : ihigh is not implemented yet";
        return;
    }
}
