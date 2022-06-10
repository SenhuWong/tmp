#pragma once
#include<cstdlib>
#include<memory>
#include"codetypes.h"
class MeshBlock;

class ADT
{
private:
	int ndim = 4;
	int nelem = 0;
	int* adtIntegers = NULL;//The structure of the tree [node_element,leftChild,rightChild,parent];
	double* adtReals = NULL;//The boudning box of the subregion of tree node[min,min,max,max]
	double* adtExtents = NULL;//The bounding box of the region of root node[min,max,min,max]
	double* coord = NULL;//The bounding box of cell at the tree node
public:

	ADT() 
	{

	}
	~ADT()
	{
		if (adtIntegers) TIOGA_FREE(adtIntegers);
		if (adtReals) TIOGA_FREE(adtReals);
		if (adtExtents) TIOGA_FREE(adtExtents);
	}
	void setDim(int dim)
	{
		if (dim == 2)
		{
			ndim = 4;
		}
		if (dim == 3)
		{
			ndim = 6;
		}
	}
	void clearData()
	{
		if (adtIntegers) delete[] adtIntegers;
		if (adtReals) delete[] adtReals;
		if (adtExtents) delete[] adtExtents;
		adtIntegers = NULL;
		adtReals = NULL;
		adtExtents = NULL;
	}

	void buildADT(int d, int nelements, double* elementBbox);
	void searchADT(MeshBlock* mb, int* cellindx, double* xsearch);

	void searchNearest(double* minD,double* x_search,int cur_proc);
};