#pragma once
#include<string>
#include"CellBrick.h"
class BlockInfo
{
public:
	int d_dim = 2;
	int ntypes = 0;//
	int* eachNodeCount = NULL;//Unchanged
	int* eachCellCount = NULL;//Unchanged
	int nnodes = 0;//
	int nedges = 0;
	int ncells = 0;//
	double* rxyz = NULL;//
	int* ibl = NULL;//Will be modifies by tioga
	int* wbc = NULL;//Unchanged
	int* obc = NULL;//Unchanged
	int nwbc = 0;
	int nobc = 0;
	int** vconn = NULL;//Unchanged

	int nvars = 0;
	double* q = NULL;
	int ncellDistributed = 0;
	int* cellDistributed = NULL;//Unchanged

	int gNwbc = 0;
	int gNobc = 0;
	int* gWbcNode = NULL;//Unchanged
	int* gObcNode = NULL;//Unchanged
	double* gWxyz = NULL;//Unchanged
	double* gOxyz = NULL;//Unchanged

	BlockInfo()
	{}
	~BlockInfo()//This might not be that necessary cuz tioga will manage it
	{}
	void setDim(int idim)
	{
		d_dim = idim;
	}
	void make_iblank();
	//All these are not necessary
	void move_by(double dx, double dy, double dz)
	{
		if (d_dim == 2)
		{
			for (int i = 0; i < nnodes; i++)
			{
				rxyz[2 * i] += dx;
				rxyz[2 * i + 1] += dy;
			}
			for (int i = 0; i < gNwbc; i++)
			{
				gWxyz[2 * i] += dx;
				gWxyz[2 * i + 1] += dy;
			}
			for (int i = 0; i < gNobc; i++)
			{
				gOxyz[2 * i] += dx;
				gOxyz[2 * i + 1] += dy;
			}
		}
		else if (d_dim == 3)
		{
			for (int i = 0; i < nnodes; i++)
			{
				rxyz[3 * i] += dx;
				rxyz[3 * i + 1] += dy;
				rxyz[3 * i + 2] += dz;
			}
			for (int i = 0; i < gNwbc; i++)
			{
				gWxyz[3 * i] += dx;
				gWxyz[3 * i + 1] += dy;
				gWxyz[3 * i + 2] += dz;
			}
			for (int i = 0; i < gNobc; i++)
			{
				gOxyz[3 * i] += dx;
				gOxyz[3 * i + 1] += dy;
				gOxyz[3 * i + 2] += dz;
			}
		}
	}
	void write_grd(const std::string& filename, int meshtag, int cur_proc);
	void write_scatter_wbc(const std::string& filename, int meshtag, int cur_proc);
	void write_scatter_obc(const std::string& filename, int meshtag, int cur_proc);
	void write_scatter_gwbc(const std::string& filename, int meshtag, int cur_proc);
	void write_scatter_gobc(const std::string& filename, int meshtag, int cur_proc);
	void write_info(const std::string& filename, int meshtag, int cur_proc);
};