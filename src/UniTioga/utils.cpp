#include<iostream>
#include"codetypes.h"
#include"../mathTools/vector3d.h"
#define TESTING false
class SquareMaTrix2D
{
public:
    int rank = 2;
    double elements[2][2];
    SquareMaTrix2D(){}
    SquareMaTrix2D(double a00,double a01,double a10,double a11) 
    {
        elements[0][0] = a00;
        elements[0][1] = a01;
        elements[1][0] = a10;
        elements[1][1] = a11;
    }
    SquareMaTrix2D(const SquareMaTrix2D& another)
    {
        elements[0][0] = another.elements[0][0];
        elements[0][1] = another.elements[0][1];
        elements[1][0] = another.elements[1][0];
        elements[1][1] = another.elements[1][1];
    }
    SquareMaTrix2D dotProduct(const SquareMaTrix2D& another)
    {
        SquareMaTrix2D result(0, 0, 0, 0);
        result.elements[0][0] = this->elements[0][0] * another.elements[0][0] + this->elements[0][1] * another.elements[1][0];
        result.elements[0][1] = this->elements[0][0] * another.elements[0][1] + this->elements[0][1] * another.elements[1][1];
        result.elements[1][0] = this->elements[1][0] * another.elements[0][0] + this->elements[1][1] * another.elements[1][0];
        result.elements[1][1] = this->elements[1][0] * another.elements[0][1] + this->elements[1][1] * another.elements[1][1];
        return result;
    }
    void show()
    {
        std::cout << "Matrix result is :\n";
        std::cout << elements[0][0] << '\t' << elements[1][0] << '\n';
        std::cout << elements[1][0] << '\t' << elements[1][1] << '\n';
    }
};
void kaiser_wrap(double** matrixA, const int& nrows, const int& n, double* eigenv, double* trace, double* sume, int* ier);

void findOBB(int dim,double* x, double xc[3], double dxc[3], double vec[3][3], int nnodes)
{
    if (dim > 3) 
    {
        std::cout << "4D is not my business\n";
        return;
    }
    double** aa;
    double* eigenv;
    double trace, sume;
    int ier;
    int id;
    double xd[3];
    double xmin[3], xmax[3];

    for (int i = 0; i < dim; i++)
    {
        xc[i] = 0;
    }
    // find centroid coordinates (not true centroid)
    //
    for (int i = 0; i < nnodes; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            xc[j] += x[dim * i + j];
        }
    }
    //
    for (int i = 0; i < dim; i++)
    {
        xc[i] = xc[i] / nnodes;
    }

    if (nnodes <= dim)
    {
        for (int i = 0; i < dim; i++)
        {
            for (int j = 0; j < dim; j++)
            {
                vec[i][j] = 0;
            }
        }
        for (int i = 0; i < dim; i++)
        {
            vec[i][i] = 1;
        }
        if (nnodes == 1)
        {
            for (int i = 0; i < dim; i++)
            {
                dxc[i] = 1e-3;
            }
            return;
        }
        else if (nnodes == 2)
        {
            for (int i = 0; i < dim; i++)
            {
                dxc[i] = TIOGA_MAX(1e-3, fabs(x[dim + i] - x[i])) * 0.5;

            }
            return;
        }
        else//nnodes is 3 which suggest this is 3D case
        {
            for (int i = 0; i < nnodes; i++)
            {
                id = dim * i;
                for (int j = 0; j < dim; j++)
                    dxc[j] = TIOGA_MAX(1e-3, fabs(x[id + j] - x[0]));
            }
            return;
        }
    }
    //
    // find co-variance matrix
    // aa = [I11 I12 I13;I21 I22 I23;I31 I32 I33]
    //
    aa = new double*[dim];
    for (int i = 0; i < dim; i++)
    {
        aa[i] = new double[dim];
    }
    eigenv = new double[dim];
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            aa[i][j] = 0;
        }
    }
    
    for (int i = 0; i < nnodes; i++)
    {
       /* aa[0][0] += (x[2 * i] - xc[0]) * (x[2 * i] - xc[0]);
        aa[0][1] += (x[2 * i] - xc[0]) * (x[2 * i + 1] - xc[1]);
        aa[1][0] += (x[2 * i + 1] - xc[1]) * (x[2 * i] - xc[0]);
        aa[1][1] += (x[2 * i + 1] - xc[1]) * (x[2 * i + 1] - xc[1]);*/
        int id = dim * i;
        for (int j = 0; j < dim; j++)
        {
            for (int k = j; k < dim; k++)
            {
                aa[j][k] += (x[id + j] - xc[j]) * (x[id + k] - xc[k]);
            }
        }
    }
    
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < i; j++)
        {
            aa[i][j] = aa[j][i];
        }
    }
    double original_aa[3][3];
    if (TESTING)
    {

        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                original_aa[i][j] = aa[i][j];
            }
        }
    }
    // use kaisers method to estimate
    // eigen values and vectors of the covariance matrix
    //

    kaiser_wrap(aa, dim, dim, eigenv, &trace, &sume, &ier);
    std::cout << "Error message is :" << ier<<'\n';
    //
    // copy the eigen vector basis on to vec
    //
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            vec[i][j] = aa[i][j];//original fortran code is column major ,so vec[i][:] is a vector
        }
    }
    //Here is a recalculate
    if (TESTING)
    {
        if (false)
        {
            //Do the right multiplication AA*T:

            double a00 = original_aa[0][0] * vec[0][0] + original_aa[0][1] * vec[1][0];
            double a01 = original_aa[0][0] * vec[0][1] + original_aa[0][1] * vec[1][1];
            double a10 = original_aa[1][0] * vec[0][0] + original_aa[1][1] * vec[1][0];
            double a11 = original_aa[1][0] * vec[0][1] + original_aa[1][1] * vec[0][1];
            double fin00 = a00 * vec[0][0] + a10 * vec[1][0];
            double fin01 = a00 * vec[0][1] + a10 * vec[1][1];
            double fin10 = a01 * vec[0][0] + a11 * vec[1][0];
            double fin11 = a01 * vec[0][1] + a11 * vec[1][1];
            if (abs(fin01) > 1.0e-6 or abs(fin10) > 1.0e-6)
            {
                std::cout << "Definately something wrong\n";
                std::cout << fin01 << "\t" << fin10 << '\n';
                std::cout << fin00 << '\t' << fin11 << '\n';
            }
        }
        if (false)
        {
            double a00 = original_aa[0][0] * vec[0][0] + original_aa[0][1] * vec[0][1];
            double a01 = original_aa[0][0] * vec[1][0] + original_aa[0][1] * vec[1][1];
            double a10 = original_aa[1][0] * vec[0][0] + original_aa[1][1] * vec[0][1];
            double a11 = original_aa[1][0] * vec[1][0] + original_aa[1][1] * vec[1][0];
            double fin00 = a00 * vec[0][0] + a10 * vec[0][1];
            double fin01 = a00 * vec[1][0] + a10 * vec[1][1];
            double fin10 = a01 * vec[0][0] + a11 * vec[0][1];
            double fin11 = a01 * vec[1][0] + a11 * vec[1][1];
            if (abs(fin01) > 1.0e-6 or abs(fin10) > 1.0e-6)
            {
                std::cout << "Definately something wrong\n";
                std::cout << fin01 << "\t" << fin10 << '\n';
                std::cout << fin00 << '\t' << fin11 << '\n';
            }

        }
        if (false)
        {
            SquareMaTrix2D A(original_aa[0][0], original_aa[0][1], original_aa[1][0], original_aa[1][1]);
            SquareMaTrix2D T(vec[0][0], vec[0][1], vec[1][0], vec[1][1]);
            SquareMaTrix2D TT(vec[0][0], vec[1][0], vec[0][1], vec[1][1]);
            SquareMaTrix2D AT = A.dotProduct(T);
            SquareMaTrix2D TTAT = TT.dotProduct(AT);
            TTAT.show();
        }
    }
    // find min and max bounds in the bounding box
    // vector basis
    //
    for (int i = 0; i < dim; i++)
    {
        xmax[i] = -BIGVALUE;
        xmin[i] = BIGVALUE;
    }
    for (int i = 0; i < nnodes; i++)
    {
        id = i * dim;
        for (int j = 0; j < dim; j++) xd[j] = 0;
        //
        for (int j = 0; j < dim; j++)
        {
            for (int k = 0; k < dim; k++)
            {
                xd[j] += (x[id + k] - xc[k]) * vec[j][k];
            }
        }
        //
        for (int j = 0; j < dim; j++)
        {
            xmax[j] = TIOGA_MAX(xmax[j], xd[j]);
            xmin[j] = TIOGA_MIN(xmin[j], xd[j]);
        }
    }
    //
    // find the extents of the box
    // and coordinates of the center w.r.t. xc
    // increase extents by 1% for tolerance
    //
    for (int i = 0; i < dim; i++)
    {
        dxc[i] = (xmax[i] - xmin[i]) * 0.5 * 1.01;
        xd[i] = (xmax[i] + xmin[i]) * 0.5;
    }
    //
    // find the center of the box in 
    // actual cartesian coordinates
    //
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            xc[i] += (xd[j] * vec[j][i]);
        }
    }
    //
    for (int i = 0; i < dim; i++)
    {
        delete[] aa[i];
    }
    delete[] aa;
    delete[] eigenv;
}

int obbIntersectCheck2D(double vA[3][3], double xA[3], double dxA[3],
    double vB[3][3], double xB[3], double dxB[3])
{
    int iflag;
    int i, j, k;
    int i1, i2, j1, j2;
    double r, r0, r1;
    double d1, d2;
    double eps = 1e-12;
    double D[3];
    double c[3][3];
    //
    // D=distance between centers
    // C=scalar product of axes
    //
    for (i = 0; i < 2; i++) D[i] = xB[i] - xA[i];
   
    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
        {
            c[i][j] = 0;
            for (k = 0; k < 2; k++)
                c[i][j] = c[i][j] + vA[i][k] * vB[j][k];
        }
    //
    // separating axes based on the faces of box A
    //
    for (i = 0; i < 2; i++)
    {
        r0 = dxA[i];
        r1 = 0;
        r = 0;
        for (j = 0; j < 2; j++)
        {
            r1 += dxB[j] * fabs(c[i][j]);
            r += fabs(vA[i][j]) * D[j];
        }
        if (r > (r0 + r1 + eps)) return 0;
    }
    //
    // separating axes based on the faces of box B
    //
    for (i = 0; i < 2; i++)
    {
        r1 = dxB[i];
        r0 = 0;
        r = 0;
        for (j = 0; j < 2; j++)
        {
            r0 += dxA[j] * fabs(c[j][i]);
            r += fabs(vB[i][j]) * D[j];
        }
        if (r > (r0 + r1 + eps)) return 0;
    }

    // return zero if no separation can be found
    //
    return 1;
}

int obbIntersectCheck3D(double vA[3][3], double xA[3], double dxA[3],
    double vB[3][3], double xB[3], double dxB[3])
{
    int iflag;
    int i, j, k;
    int i1, i2, j1, j2;
    double r, r0, r1;
    double d1, d2;
    double eps = 1e-12;
    double D[3];
    double c[3][3];
    //
    // D=distance between centers
    // C=scalar product of axes
    //
    for (i = 0; i < 3; i++) D[i] = xB[i] - xA[i];
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
        {
            c[i][j] = 0;
            for (k = 0; k < 3; k++)
                c[i][j] = c[i][j] + vA[i][k] * vB[j][k];
        }
    //
    // separating axes based on the faces of box A
    //
    for (i = 0; i < 3; i++)
    {
        r0 = dxA[i];
        r1 = 0;
        r = 0;
        for (j = 0; j < 3; j++)
        {
            r1 += dxB[j] * fabs(c[i][j]);
            r += fabs(vA[i][j]) * D[j];
        }
        if (r > (r0 + r1 + eps)) return 0;
    }
    //
    // separating axes based on the faces of box B
    //
    for (i = 0; i < 3; i++)
    {
        r1 = dxB[i];
        r0 = 0;
        r = 0;
        for (j = 0; j < 3; j++)
        {
            r0 += dxA[j] * fabs(c[j][i]);
            r += fabs(vB[i][j]) * D[j];
        }
        if (r > (r0 + r1 + eps)) return 0;
    }
    //
    // cross products
    //
    for (i = 0; i < 3; i++)
    {
        i1 = (i + 1) % 3;
        i2 = (i + 2) % 3;
        for (j = 0; j < 3; j++)
        {
            j1 = (j + 1) % 3;
            j2 = (j + 2) % 3;

            r0 = dxA[i1] * fabs(c[i2][j]) + dxA[i2] * fabs(c[i1][j]);
            r1 = dxB[j1] * fabs(c[i][j2]) + dxB[j2] * fabs(c[i][j1]);

            d2 = 0;
            d1 = 0;
            for (k = 0; k < 3; k++)
            {
                d2 += vA[i2][k] * D[k];
                d1 += vA[i1][k] * D[k];
            }

            r = fabs(c[i1][j] * d2 - c[i2][j] * d1);

            if (r > (r0 + r1 + eps)) {
                return 0;
            }
        }
    }
    //
    // return zero if no separation can be found
    //
    return 1;
}
int obbIntersectCheck(int dim, double vA[3][3], double xA[3], double dxA[3],
    double vB[3][3], double xB[3], double dxB[3])
{
    if (dim == 2)
    {
        return obbIntersectCheck2D(vA, xA, dxA, vB, xB, dxB);
    }
    else
    {
        return obbIntersectCheck3D(vA, xA, dxA, vB, xB, dxB);
    }
}
void getobbcoords(int dim, double* xc, double* dxc, double vec[3][3], double xv[8][3])
{
    int i, j, k, ik;
    
    for (i = 0; i < pow(2,dim); i++)
    {
        for (k = 0; k < dim; k++)
            xv[i][k] = xc[k];
        for (k = 0; k < dim; k++)
        {
            ik = (2 * ((i & (1 << k)) >> k) - 1);
            for (j = 0; j < dim; j++)
                xv[i][j] += ik * dxc[k] * vec[k][j];
        }
    }

}
void transform2OBB(int dim, double* xv, double* xc, double vec[3][3], double* xd)
{
    int j, k;
    for (j = 0; j < dim; j++)
    {
        xd[j] = 0;
        for (k = 0; k < dim; k++)
            xd[j] += (xv[k] - xc[k]) * vec[j][k];
    }
}


//modify ADT builder to remove common nodes
void uniqNodesTree(double* coord,
    int* itag, double* rtag, int* meshtag,
    int* elementsAvailable,
    int ndim, int nav)
{
    int ibox;
    int p1, p2;
    int* tmpint;
    int npts[8], iv[3], cft[9];
    double xmax[3], xmin[3], xmid[3], dx[3], xp[3];
    int icheck = 1;
    //
    // if there are more than 10 elements divide the tree
    //
    if (nav > 20) 
    {
        //
        // find the bound of the boxes
        //
        icheck = 0;
        xmin[0] = xmin[1] = xmin[2] = BIGVALUE;
        xmax[0] = xmax[1] = xmax[2] = -BIGVALUE;
        for (int i = 0; i < nav; i++)
        {
            for (int j = 0; j < ndim; j++)
            {
                xmin[j] = TIOGA_MIN(xmin[j], coord[ndim * elementsAvailable[i] + j]);
                xmax[j] = TIOGA_MAX(xmax[j], coord[ndim * elementsAvailable[i] + j]);
            }
        }
        for (int j = 0; j < ndim; j++) 
        {
            xmid[j] = (xmax[j] + xmin[j]) * 0.5;
            dx[j] = (xmax[j] - xmin[j]) * 0.5 + TOL;
        }
        for (int j = 0; j < 8; j++) npts[j] = 0;
        for (int i = 0; i < nav; i++)
        {
            for (int j = 0; j < ndim; j++) 
            {
                xp[j] = coord[ndim * elementsAvailable[i] + j] - xmid[j];
                iv[j] = floor(xp[j] / dx[j]) + 1;
            }
            if (ndim == 3)
            {
                ibox = 4 * iv[0] + 2 * iv[1] + iv[2];
            }
            else if (ndim == 2)
            {
                ibox = 2 * iv[0] + iv[1];
            }
            else
            {
                throw std::runtime_error("Dimension error\n");
            }
            npts[ibox]++;
        }
        for (int j = 0; j < 8; j++) if (npts[j] == nav) icheck = 1;
        if (!icheck) 
        {
            cft[0] = 0;
            for (int j = 0; j < pow(2,ndim); j++)
                cft[j + 1] = cft[j] + npts[j];
            tmpint = new int[nav];
            for (int i = 0; i < nav; i++)
            {
                for (int j = 0; j < ndim; j++) 
                {
                    xp[j] = coord[ndim * elementsAvailable[i] + j] - xmid[j];
                    iv[j] = floor(xp[j] / dx[j]) + 1;
                }
                if (ndim == 2)
                {
                    ibox = 2 * iv[0] + iv[1];
                }
                else if(ndim==3)
                {
                    ibox = 4 * iv[0] + 2 * iv[1] + iv[2];
                }
                else
                {
                    throw std::runtime_error("Dimension error\n");
                }
                
                npts[ibox] = npts[ibox] - 1;
                tmpint[npts[ibox] + cft[ibox]] = elementsAvailable[i];
            }
            for (int i = 0; i < nav; i++)
                elementsAvailable[i] = tmpint[i];
            TIOGA_FREE(tmpint);
            for (int j = 0; j < pow(2,ndim); j++)
                if (cft[j + 1] > cft[j])
                    uniqNodesTree(coord, itag, rtag, meshtag, &(elementsAvailable[cft[j]]),
                        ndim, cft[j + 1] - cft[j]);
        }
    }
    if (icheck) 
    {
        for (int i = 0; i < nav; i++)
        {
            p1 = elementsAvailable[i];
            for (int j = i + 1; j < nav; j++)
            {
                p2 = elementsAvailable[j];
                if (ndim == 2)
                {
                    if (fabs(coord[2 * p1] - coord[2 * p2]) + fabs(coord[2 * p1 + 1] - coord[2 * p2 + 1]) < TOL && meshtag[p1] == meshtag[p2])
                    {
                        if (p1 > p2) 
                        {
                            rtag[p2] = TIOGA_MAX(rtag[p1], rtag[p2]);
                            itag[p1] = itag[p2];
                        }
                        else 
                        {
                            rtag[p1] = TIOGA_MAX(rtag[p1], rtag[p2]);
                            itag[p2] = itag[p1];
                        }
                    }
                }
                if (ndim == 3)
                {
                    if (fabs(coord[3 * p1] - coord[3 * p2]) +
                        fabs(coord[3 * p1 + 1] - coord[3 * p2 + 1]) +
                        fabs(coord[3 * p1 + 2] - coord[3 * p2 + 2]) < TOL &&
                        meshtag[p1] == meshtag[p2])
                    {
                        if (p1 > p2) {
                            rtag[p2] = TIOGA_MAX(rtag[p1], rtag[p2]);
                            itag[p1] = itag[p2];
                        }
                        else {
                            rtag[p1] = TIOGA_MAX(rtag[p1], rtag[p2]);
                            itag[p2] = itag[p1];
                        }
                    }
                }
            }
        }
    }
}
/*
 * Create a unique hash for list of coordinates with duplicates in
 * them. Find the rtag as max of all duplicate samples. itag contains
 * the hash to the real point
 */
void uniquenodes_octree(int dim, double* x, int* meshtag, double* rtag, int* itag,
    int* nn)
{
    int nelem = *nn;
    int* elementsAvailable = new int[nelem];
    for (int i = 0; i < nelem; i++)
    {
        elementsAvailable[i] = i;
        itag[i] = i;
    }
    uniqNodesTree(x, itag, rtag, meshtag, elementsAvailable, dim, nelem);
    TIOGA_FREE(elementsAvailable);
}


int checkHoleMap(int dim, double* x, int* nx, int* sam, double* extents)
{
    int i;
    int mm;
    double dx[3];
    int ix[3];

    for (i = 0; i < dim; i++)
    {
        dx[i] = (extents[dim + i] - extents[i]) / nx[i];
    }
    for (i = 0; i < dim; i++)
    {
        ix[i] = (x[i] - extents[i]) / dx[i];
        if (ix[i] < 0 or ix[i] > nx[i] - 1) { return 0; }
    }
    if (dim == 2)
    {
        mm = ix[1] * nx[0] + ix[0];
    }
    else if (dim == 3)
    {
        mm = ix[2] * nx[1] * nx[0] + ix[1] * nx[0] + ix[0];
    }
    return sam[mm];
}

/**
 fill a given hole map using iterative
 flood fill from outside the marked boundary.
 boundary is marked by "2"
*/
void fillHoleMap3D(int* holeMap, int ix[3], int isym)
{
    int m;
    int ii, jj, kk, mm;
    int ns2;
    int i, j, k;
    int ipaint, npaint, nneig;
    int mk[6];
    //
    // now start from outside and paint the
    // the exterior
    //
    ns2 = ix[0] * ix[1];
    //
    for (kk = 0; kk < ix[2]; kk += (ix[2] - 1))//Along the 2nd direction ,all holemap value is set to 1.
        for (jj = 0; jj < ix[1]; jj++)
            for (ii = 0; ii < ix[0]; ii++)
            {
                mm = kk * ns2 + jj * ix[0] + ii;
                holeMap[mm] = 1;
            }
    for (kk = 0; kk < ix[2]; kk++)
        for (jj = 0; jj < ix[1]; jj += (ix[1] - 1))
            for (ii = 0; ii < ix[0]; ii++)
            {
                mm = kk * ns2 + jj * ix[0] + ii;
                holeMap[mm] = 1;
            }
    for (kk = 0; kk < ix[2]; kk++)
        for (jj = 0; jj < ix[1]; jj++)
            for (ii = 0; ii < ix[0]; ii += (ix[0] - 1))
            {
                mm = kk * ns2 + jj * ix[0] + ii;
                holeMap[mm] = 1;
            }
    npaint = ns2 * ix[2];
    while (npaint > 0)
    {
        npaint = 0;
        for (k = 1; k < ix[2] - 1; k++)
            for (j = 1; j < ix[1] - 1; j++)
                for (i = 1; i < ix[0] - 1; i++)
                {
                    m = k * ns2 + j * ix[0] + i;
                    if (holeMap[m] == 0)
                    {
                        ipaint = 0;
                        if (isym == 1)
                        {
                            mk[0] = m - ns2;
                            mk[1] = m + ns2;
                            mk[2] = m - ix[0];
                            mk[3] = m + ix[0];
                            mk[4] = m - 1;
                            mk[5] = m + 1;
                            nneig = 4;
                        }
                        else if (isym == 2)
                        {
                            mk[0] = m - ns2;
                            mk[1] = m + ns2;
                            mk[4] = m - ix[0];
                            mk[5] = m + ix[0];
                            mk[2] = m - 1;
                            mk[3] = m + 1;
                            nneig = 4;
                        }
                        else if (isym == 3)
                        {
                            mk[4] = m - ns2;
                            mk[5] = m + ns2;
                            mk[0] = m - ix[0];
                            mk[1] = m + ix[0];
                            mk[2] = m - 1;
                            mk[3] = m + 1;
                            nneig = 4;
                        }
                        else
                        {
                            mk[0] = m - ns2;
                            mk[1] = m + ns2;
                            mk[2] = m - ix[0];
                            mk[3] = m + ix[0];
                            mk[4] = m - 1;
                            mk[5] = m + 1;
                            nneig = 6;
                        }
                        for (kk = 0; kk < nneig && ipaint == 0; kk++)
                        {
                            ipaint = (ipaint || holeMap[mk[kk]] == 1);
                        }
                        if (ipaint > 0)
                        {
                            holeMap[m] = 1;
                            npaint++;
                        }
                    }
                }
    }
    for (i = 0; i < ix[2] * ix[1] * ix[0]; i++)
    {
        holeMap[i] = (holeMap[i] == 0 || holeMap[i] == 2);
    }
}

void fillHoleMap2D(int* holeMap, int ix[2], int isym)
{
    int m;
    int ii, jj, kk, mm;
    int i, j;
    int ipaint, npaint, nneig;
    int mk[4];

    for (jj = 0; jj < ix[1]; jj += (ix[1] - 1))
    {
        for (ii = 0; ii < ix[0]; ii++)
        {
            mm = jj * ix[0] + ii;
            holeMap[mm] = 1;
        }
    }
    for (jj = 0; jj < ix[1]; jj++)
    {
        for (ii = 0; ii < ix[0]; ii += (ix[0] - 1))
        {
            mm = jj * ix[0] + ii;
            holeMap[mm] = 1;
        }
    }
    npaint = ix[1] * ix[0];
    while (npaint > 0)
    {
        npaint = 0;
        for (j = 1; j < ix[1] - 1; j++)
        {
            for (i = 1; i < ix[0] - 1; i++)
            {
                m = j * ix[0] + i;
                if (holeMap[m] == 0)
                {
                    ipaint = 0;
                    if (isym == 1)
                    {
                        mk[0] = m - ix[0];
                        mk[1] = m + ix[0];
                        mk[2] = m - 1;
                        mk[3] = m + 1;
                        nneig = 2;
                    }
                    else if (isym == 2)
                    {
                        mk[0] = m - 1;
                        mk[1] = m + 1;
                        mk[2] = m - ix[0];
                        mk[3] = m + ix[0];
                        nneig = 2;
                    }
                    else
                    {
                        mk[0] = m - ix[0];
                        mk[1] = m + ix[0];
                        mk[2] = m - 1;
                        mk[3] = m + 1;
                        nneig = 4;
                    }
                    for (kk = 0; kk < nneig && ipaint == 0; kk++)
                    {
                        ipaint = (ipaint || holeMap[mk[kk]] == 1);
                    }
                    if (ipaint > 0)
                    {
                        holeMap[m] = 1;
                        npaint++;
                    }
                }
            }
        }
    }
    for (i = 0; i < ix[1] * ix[0]; i++)
    {
        holeMap[i] = (holeMap[i] == 0 || holeMap[i] == 2);
    }
}

void fillHoleMap(int dim,int* holeMap, int ix[3], int isym)
{
    if (dim == 2)
    {
        fillHoleMap2D(holeMap, ix, isym);
    }
    else if (dim == 3)
    {
        fillHoleMap3D(holeMap, ix, isym);
    }
}