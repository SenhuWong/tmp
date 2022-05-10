#include<fstream>
#include<iostream>
void get_amr_index_xyz_2d(int nq, int i, int j,
                          int pBasis,
                          int nX, int nY,
                          int nf,
                          double *xlo,
                          double *dx,
                          double qnodes[],
                          int *index,
                          double *xyz)
{
    int two_nf = 2 * nf;
    int qp_1d = pBasis + 1;

    int x_stride = nq * (nX + two_nf) * (nY + two_nf);//The stride of q and that of a grid
    int y_stride = x_stride * (qp_1d);
    //This is the cube inside the patch given by (i,j,k)
    double cell_xlo = xlo[0] + i * dx[0];
    double cell_ylo = xlo[1] + j * dx[1];

    double half_dx = 0.5 * dx[0];
    double half_dy = 0.5 * dx[1];

    //This is the starting entry given by (i,j,k)
    int starting_entry = (i + nf) +
                         (nX + two_nf) * (j + nf);
    double xco = 0;
    double yco = 0;

    int temp_count = 0;
    //Looping through xyz to find amr coordinates and index
    //qnodes[] therefore should be a number from -1 to 1
    for (int y = 0; y < qp_1d; y++)
    {
        yco = cell_ylo + half_dy * (1.0 + qnodes[y]);
        for (int x = 0; x < qp_1d; x++)
        {
            xco = cell_xlo + half_dx * (1.0 + qnodes[x]);
            index[temp_count] = starting_entry + y_stride * (y) + x_stride * (x);
            xyz[2 * temp_count + 0] = xco;
            xyz[2 * temp_count + 1] = yco;
            temp_count++;
        }
    }
}

void amr_index_to_ijklmn_2d(int pBasis, int nX, int nY, int nf, int nq, int index, int *ijklmn)
{
    int remainder;//127
    int two_nf = 2 * nf;//2
    int qp_1d = pBasis + 1;//2
    int i_size = (nX + two_nf);//6
    int j_size = (nY + two_nf);//4

    int i_stride = i_size;//6
    int ij_stride = i_stride * j_size;//24

    int grid_stride = nq * ij_stride;//120

    int x_stride = grid_stride;//120
    int y_stride = x_stride * (qp_1d);//240

    remainder = index;

    int m = remainder / y_stride;//0
    remainder = remainder % y_stride;//127

    int l = remainder / x_stride;//1
    remainder = remainder % x_stride;//7

    //Skip fields,we do not locate for q*grid_stride
    remainder = remainder % ij_stride;

    int j = remainder / i_stride;
    remainder = remainder % i_stride;
    int i = remainder;

    ijklmn[0] = i - nf;
    ijklmn[1] = j - nf;
    ijklmn[2] = l;
    ijklmn[3] = m;
}
void get_amr_index_xyz_3d(int nq, int i, int j, int k,
                          int pBasis,
                          int nX, int nY, int nZ,
                          int nf,
                          double *xlo,
                          double *dx,
                          double qnodes[],
                          int *index,
                          double *xyz)
{
    int two_nf = 2 * nf;
    int qp_1d = pBasis + 1;

    int x_stride = nq * (nX + two_nf) * (nY + two_nf) * (nZ + two_nf);
    int y_stride = x_stride * (qp_1d);
    int z_stride = y_stride * (qp_1d);
    //This is the cube inside the patch given by (i,j,k)
    double cell_xlo = xlo[0] + i * dx[0];
    double cell_ylo = xlo[1] + j * dx[1];
    double cell_zlo = xlo[2] + k * dx[2];
    double half_dx = 0.5 * dx[0];
    double half_dy = 0.5 * dx[1];
    double half_dz = 0.5 * dx[2];
    //This is the starting entry given by (i,j,k)
    int starting_entry = (i + nf) +
                         (nX + two_nf) * (j + nf) +
                         (nY + two_nf) * (nX + two_nf) * (k + nf);
    double xco = 0;
    double yco = 0;
    double zco = 0;
    int temp_count = 0;
    //Looping through xyz to find amr coordinates and index
    for (int z = 0; z < qp_1d; z++)
    {
        zco = cell_zlo + half_dz * (1.0 + qnodes[z]); //qnodes[] therefore should be a number from -1 to 1
        for (int y = 0; y < qp_1d; y++)
        {
            yco = cell_ylo + half_dy * (1.0 + qnodes[y]);
            for (int x = 0; x < qp_1d; x++)
            {
                xco = cell_xlo + half_dx * (1.0 + qnodes[x]);
                index[temp_count] = starting_entry + z_stride * (z) + y_stride * (y) + x_stride * (x);
                //The INDEX from lower(moves faster) to upper(moves slower) is
                //1.i,j,k
                //2.q
                //3.x,y,z
                //
                //This has something to do with the data storage at MeshBlock I guess
                //Q[p+1,p+1,p+1,nq,nZ+2*nf,nY+2*nf,nX+2*nf]-->C++ Storage
                xyz[3 * temp_count + 0] = xco;
                xyz[3 * temp_count + 1] = yco;
                xyz[3 * temp_count + 2] = zco;
                temp_count++;
            }
        }
    }
}
void amr_index_to_ijklmn_3d(int pBasis, int nX, int nY, int nZ, int nf, int nq, int index, int *ijklmn)
{
    int remainder;
    int two_nf = 2 * nf;
    int qp_1d = pBasis + 1;
    int i_size = (nX + two_nf);
    int j_size = (nY + two_nf);
    int k_size = (nZ + two_nf);

    int i_stride = i_size;
    int ij_stride = i_stride * j_size;
    int ijk_stride = ij_stride * k_size;

    int grid_stride = nq * ijk_stride;

    int x_stride = grid_stride;
    int y_stride = x_stride * (qp_1d);
    int z_stride = y_stride * (qp_1d);

    remainder = index;
    int n = remainder / z_stride;
    remainder = remainder % z_stride;

    int m = remainder / y_stride;
    remainder = remainder % y_stride;

    int l = remainder / x_stride;
    remainder = remainder % x_stride;

    //Skip fields,we do not locate for q*grid_stride
    remainder = remainder % ijk_stride;

    int k = remainder / ij_stride;
    remainder = remainder % ij_stride;

    int j = remainder / i_stride;
    remainder = remainder % i_stride;
    int i = remainder;

    ijklmn[0] = i - nf;
    ijklmn[1] = j - nf;
    ijklmn[2] = k - nf;
    ijklmn[3] = l;
    ijklmn[4] = m;
    ijklmn[5] = n;
}
void get_amr_index_xyz(int dim, int nq, int i, int j, int k,
                       int pBasis,
                       int nX, int nY, int nZ,
                       int nf,
                       double *xlo,
                       double *dx,
                       double qnodes[],
                       int *index,
                       double *xyz)
{
    if (dim == 2)
    {
        get_amr_index_xyz_2d(nq, i, j, pBasis, nX, nY, nf, xlo, dx, qnodes, index, xyz);
    }
    else if (dim == 3)
    {
        get_amr_index_xyz_3d(nq, i, j, k, pBasis, nX, nY, nZ, nf, xlo, dx, qnodes, index, xyz);
    }
}

void amr_index_to_ijklmn(int dim, int pBasis, int nX, int nY, int nZ, int nf, int nq, int index, int *ijklmn)
{
    if (dim == 2)
    {
        amr_index_to_ijklmn_2d(pBasis, nX, nY, nf, nq, index, ijklmn);
    }
    else if (dim == 3)
    {
        amr_index_to_ijklmn_3d(pBasis, nX, nY, nZ, nf, nq, index, ijklmn);
    }
}
#include<fstream>
void test_amr_index(int dim)
{
    std::ofstream fout;
    fout.open("test_amr_index");
    int pdegree = 1; //2 point inside
    int qstride = 5;
    int nf = 1;
    double xlo[3];
    double dx[3];
    int index1[8];
    double qnodes[2];
    double xyz[8 * 3];
    for (int i = 0; i < dim; i++)
    {
        xlo[i] = 0.0;
        dx[i] = 1.0;
    }
    qnodes[0] = -0.6;
    qnodes[1] = 0.4;
    int ijklmn[6];
    int nX = 4;
    int nY = 2;
    int nZ = 2;
    for (int i = 0; i < nX; i++)
    {
        for (int j = 0; j < nY; j++)
        {
            for (int k = 0; k < nZ; k++)
            {
                get_amr_index_xyz(dim, qstride, i, j, k, pdegree, nX, nY, nZ, nf, xlo, dx, qnodes, index1, xyz);
                if (dim == 2)
                {
                    fout<<"("<<i<<","<<j<<"):\n";
                    for(int q = 0;q<4;q++)
                    {
                        amr_index_to_ijklmn(dim,pdegree,nX,nY,nZ,nf,qstride,index1[q],ijklmn);
                        fout<<'\t'<<index1[q]<<","<<xyz[2*q]<<","<<xyz[2*q+1]<<"\n";
                        fout<<'\t'<<ijklmn[0]<<","<<ijklmn[1]<<","<<ijklmn[2]<<","<<ijklmn[3]<<"\n";
                    }
                }
                else if(dim==3)
                {
                    fout<<"("<<i<<","<<j<<","<<k<<"):\n";
                    for(int q = 0;q<8;q++)
                    {
                        amr_index_to_ijklmn(dim,pdegree,nX,nY,nZ,nf,qstride,index1[q],ijklmn);
                        fout<<'\t'<<index1[q]<<","<<xyz[3*q]<<","<<xyz[3*q+1]<<","<<xyz[3*q+2]<<"\n";
                        fout<<'\t'<<ijklmn[0]<<","<<ijklmn[1]<<","<<ijklmn[2]<<","<<ijklmn[3]<<","<<ijklmn[4]<<","<<ijklmn[5]<<"\n";
                    }

                }
            }
        }
    }
    fout.close();
}