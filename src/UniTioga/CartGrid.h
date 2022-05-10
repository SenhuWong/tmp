#pragma once
#ifndef CARTGRID_H
#define CARTGRID_H

#include<cstdlib>

class CartGrid
{
private:
	int d_dim = 2;
	double xlosup[3] = { 0.0,0 };
	double* dxlvl = NULL;
	int* lcount = NULL;
	int maxlevel = 0;
public:
	int* global_id = NULL;
	int* local_id = NULL;

	int* level_num = NULL;
	int* proc_id = NULL;
	
	int* porder = NULL;
	int* ilo = NULL;
	int* ihi = NULL;
	int* ncells = NULL;
	int myid = -1;
	int nf = 0;
	int qstride = 0;
	double* xlo = NULL;
	double* dx = NULL;
	double* qnode = NULL;
	int ngrids = 0;
	void(*donor_frac)(int*, double*, int*, double*) = NULL;
	CartGrid()
	{
	}
	CartGrid(int in_dim)
		:d_dim(in_dim)
	{

	}
	~CartGrid()
	{
		if (dxlvl) delete[] dxlvl;
		if (lcount) delete[] lcount;
		if (global_id) delete[] global_id;
		if (local_id) delete[] local_id;
		if (level_num) delete[] level_num;
		if (proc_id) delete[] proc_id;
		if (porder) delete[] porder;
		if (ilo) delete[] ilo;
		if (ihi) delete[] ihi;
		if (ncells) delete[] ncells;
		if (xlo) delete[] xlo;
		if (dx) delete[] dx;
		if (qnode) delete[] qnode;
	}
	void registerData(int dim,int nfin, int qstridein, double* qnodein, int* idata, double* rdata, int ngridsin, int qnodesize);
	void preprocess();
	void setDim(int dim)
	{
		d_dim = dim;
	}
	void search(double* x, int* donorid, int nsearch);
	void setcallback(void(*f1)(int *,double *, int*,double*))
	{
		donor_frac = f1;
	}
};
#endif // !CARTGRID_H
