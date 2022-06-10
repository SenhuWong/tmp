#include "CartGrid.h"
#include "codetypes.h"
#include<iostream>
void CartGrid::registerData(int dim,int nfin, int qstridein, double *qnodein, int *idata, double *rdata, int ngridsin, int qnodesize)
{

	std::cout<<ngridsin<<"is ngridin at CartGrid::registerData \n";
	ngrids = ngridsin;
	std::cout<<ngrids<<"is ngrid at CartGrid::registerData \n";
	d_dim = dim;
	global_id = new int[ngrids];
	level_num = new int[ngrids];
	proc_id = new int[ngrids];
	ilo = new int[d_dim * ngrids];
	ihi = new int[d_dim * ngrids];
	xlo = new double[d_dim * ngrids];
	dx = new double[d_dim * ngrids];
	porder = new int[ngrids];
	local_id = new int[ngrids];
	qnode = new double[qnodesize]; //Should this size be kept;
	ncells = new int[d_dim * ngrids];
	for (int i = 0; i < qnodesize; i++)
	{
		qnode[i] = qnodein[i];
	}
	nf = nfin;
	qstride = qstridein;
	for (int i = 0; i < ngrids; i++)
	{
		int id = d_dim * i;
		int id2 = 2 * id;
		int iloc = (5 + 2 * d_dim) * i;
		//3D:5+3+3==11
		//2D:5+2+2==9
		//[gid,lvnum,]
		global_id[i] = idata[iloc];
		level_num[i] = idata[iloc + 1];
		proc_id[i] = idata[iloc + 2];
		porder[i] = idata[iloc + 3];
		local_id[i] = idata[iloc + 4];
		for (int n = 0; n < d_dim; n++)
		{
			ilo[id + n] = idata[iloc + 5 + n];
			ihi[id + n] = idata[iloc + 5 + d_dim + n];
			ncells[id + n] = ihi[id + n] - ilo[id + n] + 1; //Acoording to SAMRAI's Convention that [1,2] means mesh on 1 and 2.
			xlo[id + n] = rdata[id2 + n];
			dx[id + n] = rdata[id2 + d_dim + n];
		}
	}
}
#include<iostream>
void CartGrid::preprocess()
{
	xlosup[0] = xlosup[1] = xlosup[2] = BIGVALUE;
	maxlevel = -1;
	std::cout<<"ngrids at preprocess is "<<ngrids<<'\n';

	for (int i = 0; i < ngrids; i++)
	{
		for (int n = 0; n < d_dim; n++)
		{
			xlosup[n] = TIOGA_MIN(xlosup[n], xlo[d_dim * i + n]);
		}
		maxlevel = TIOGA_MAX(maxlevel, level_num[i]);
	}
	

	maxlevel++;
	lcount = new int[maxlevel];
	dxlvl = new double[d_dim * maxlevel];
	for (int i = 0; i < maxlevel; i++)
	{
		lcount[i] = 0;
	}
	for (int i = 0; i < ngrids; i++)
	{
		lcount[level_num[i]]++;
		for (int n = 0; n < d_dim; n++)
		{
			dxlvl[d_dim * level_num[i] + n] = dx[d_dim * i + n];
		}
	}
}
void CartGrid::search(double *x, int *donorid, int npts)
{
	bool flag;
	std::cout<<"curproc:"<<myid<<": point_total="<<npts<<'\n';
	for (int i = 0; i < npts; i++)
	{
		int id = d_dim * i;
		flag = 0;
		donorid[i] = -1;
		//First set donorId to -1, if there is a donor set to the donorid
		for (int l = maxlevel - 1; l >= 0 and flag == 0; l--)
		{
			for (int j = 0; j < ngrids and flag == 0; j++)
			{
				if (level_num[j] == l)
				{
					flag = 1;
					for (int n = 0; n < d_dim; n++)
					{
						flag = flag and (x[d_dim * i + n] - xlo[d_dim * j + n] > -TOL);
						flag = flag and (x[d_dim * i + n] - (xlo[d_dim * j + n] + dx[d_dim * j + n] * ncells[d_dim * j + n]) < TOL);
					}
					if (flag)
					{
						donorid[i] = j;
					}
				}
			}
		}
	}
}