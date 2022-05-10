#pragma once
#ifndef CARTBLOCK_H
#define CARTBLOCK_H
#include <cstdlib>
#include <fstream>
#include<list>
struct INTERPLIST2;
struct DONORLIST;
struct HOLEMAP;
class CartGrid;
class CartBlock
{
private:
	int global_id = -1;
	int local_id = -1;
	int d_dim = 2;
	int ncells[3] = {0, 0, 0};
	int nf = 0;
	int qstride = 0;
	int ndof = 0;
	int pdegree = 0;
	int pd = 1;
	int d1 = 0;
	int d2 = 0;
	int dd = 0;
	int ddnf = 0;
	int myid = 0;
	int *ibl = NULL;
	double *q = NULL;
	double *qnode = NULL;
	double xlo[3] = {0, 0, 0};
	double dx[3] = {0, 0, 0};
	int ndonors = 0;
	int interpListSize = 0;
	std::list<INTERPLIST2*> interpList;
	
	//INTERPLIST2 *listptr = NULL;
	DONORLIST **donorList = NULL;
	void (*donor_frac)(int *, double *, int *, double *) = NULL;

public:
	CartBlock();
	CartBlock(int dim_in);
	~CartBlock();
	void registerData(int dim,int local_id_in, int global_id_in, int *ibl_in, double *qin)
	{
		d_dim = dim;
		local_id = local_id_in;
		global_id = global_id_in;
		ibl = ibl_in;
		q = qin;
	}
	void setDim(int dim)
	{
		d_dim = dim;
	}
	void preprocess(CartGrid *cg);
	void getInterpolatedData(int *nints, int *nreals, int **intData, double **realData, int nvar);
	void update(double *qval, int index, int nq);
	void getCancellationData(int *cancelledData, int *ncancel);
	void processDonors(HOLEMAP *holemap, int nmesh);
	void insertInDonorList(int senderid, int index, int meshtagdonor, int remoteid, int remoteblockid, double cellRes);
	void insertInInterpList(int procid, int remoteid, int remoteblockid, double *xtmp);
	void writeCellFile(const std::string& filename);
	void clearLists();
	void initializeLists();

private:
	void processDonors2D(HOLEMAP *holemap, int nmesh);
	void processDonors3D(HOLEMAP *holemap, int nmesh);
};
#endif // !CARTBLOCK_H
