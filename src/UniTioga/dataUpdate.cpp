#include "MeshBlock.h"
#include <fstream>

#define ROW 0
#define COLUMN 1
void MeshBlock::getInterpolatedSolution(int* nints, int* nreals, int** intData, double** realData, double* q,
		int nvar, int interptype)
{
    double weight;
    double *qq = NULL;
    int icount,dcount;

    (*nints) = (*nreals) = 0;
    for(int i = 0;i<ninterp;i++)
    {
        if(!interpList[i].cancel)
        {
            (*nints)++;
            (*nreals)+=nvar;
        }
    }
    if((*nints)==0) {return;}
    qq = new double[nvar];
    *intData = new int[3*(*nints)];
    *realData = new double[*nreals];
    icount=dcount = 0;
    
    int k,m,inode;
    std::string filename = "DonorCellInds"+std::to_string(meshtag)+"_"+std::to_string(myid);
    std::ofstream fout;
    fout.open(filename);

    if (interptype==ROW)
    {
        for(int i = 0;i<ninterp;i++)
        {
            if(!interpList[i].cancel)
            {
                for(k = 0;k<nvar;k++)
                {
                    qq[k] = 0;
                }
                for(m = 0;m<interpList[i].nweights;m++)
                {
                    inode = interpList[i].inode[m];
                    for(int j = 0;j<d_dim;j++)
                    {
                        fout << x[inode*d_dim+j]<<'\t';
                    }
                    fout<<'\n';
                    weight = interpList[i].weights[m];
                    if(weight < -TOL or weight > 1.0+TOL)
                    {
                        throw std::runtime_error("weighrs are not convex 1\n");
                    }
                    int offset = inode*nvar;
                    for(k = 0;k<nvar;k++)
                    {
                        qq[k] += q[offset+k]*weight;
                    }
                }
                (*intData)[icount++]=interpList[i].receptorInfo[0];//proc
                (*intData)[icount++]=interpList[i].receptorInfo[1];//point
                (*intData)[icount++]=interpList[i].receptorInfo[2];//block
                for(k = 0;k<nvar;k++)
                {
                    (*realData)[dcount++] = qq[k];
                }


            }
        }
    }
    else if(interptype==COLUMN)
    {
        for(int i = 0;i<ninterp;i++)
        {
            if(!interpList[i].cancel)
            {
                for(k = 0;k<nvar;k++)
                {
                    qq[k] =0;
                }
                for(m = 0;m < interpList[i].nweights;m++)
                {
                    inode = interpList[i].inode[m];
                    weight = interpList[i].weights[m];
                    if(weight < -TOL or weight > 1.0+TOL)
                    {
                        throw std::runtime_error("weights are not convex 1\n");
                    }
                    for(k = 0;k<nvar;k++)
                    {
                        qq[k] += q[k*nnodes+inode]*weight;
                    }
                }
                (*intData)[icount++]=interpList[i].receptorInfo[0];
                (*intData)[icount++]=interpList[i].receptorInfo[1];
                (*intData)[icount++]=interpList[i].receptorInfo[2];
                for(k = 0;k<nvar;k++)
                {
                    (*realData)[dcount++] = qq[k];
                }
            }

        }
    }
    fout.close();
    if(qq)
    {
        TIOGA_FREE(qq);
    }
}

void MeshBlock::updateSolnData(int inode, double* qvar, double* q, int nvar, int interptype)
{
    if(interptype==ROW)
    {
        for(int k = 0;k<nvar;k++)
        {
            q[inode*nvar+k] = qvar[k];
        }
    }
    else if(interptype==COLUMN)
    {
        for(int k = 0;k<nvar;k++)
        {
            q[nnodes*k+inode] = qvar[k];

        }
    }

}