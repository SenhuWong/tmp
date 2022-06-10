#pragma once
//I need a reader and a Blockinformation organizer.
#include <vector>
#include <string>
#include "Reader.h"
#include "../samrai/SamraiFeeder.h"
#include "../UniTioga/TiogaFeeder.h"
#include "../UnstructSolver/UnstructFeeder.h"
void findClockWise(int *mapping, double xv[8][3]);

//For now ,lets read one at a time
class CobaltReader : public TiogaFeeder,
					 public SamraiFeeder,
					 public UnstructFeeder
{
public:
	const int cellBase = 1;
	const int nodeBase = 1;
	int signOver = -1; //in BoundaryUtility
	int signWall = -2;
	int signSymm = -3;
	double coord_symm = 0;
	int dim_symm = 2;
	//This is for SAMRAI's adaptive refinement

	CobaltReader() :Reader(),TiogaFeeder(),SamraiFeeder()
	{
	}
	~CobaltReader()
	{
	}
	void setWallSign(int sign)
	{
		signWall = sign;
	}
	void setOverSign(int sign)
	{
		signOver = sign;
	}
	void setSymmSign(int sign)
	{
		signSymm = sign;
	}
	void setSymmMirror(int sdim, double scoord)
	{
		dim_symm = sdim;
		coord_symm = scoord;
	}
	void setDim(int idim)
	{
		d_dim = idim;
		for (auto iter = bis.begin(); iter != bis.end(); iter++)
		{
			iter->setDim(idim);
		}
	}
	void addPath(const std::string &pathname);
	void setComm(int icur_proc, int inum_proc)
	{
		cur_proc = icur_proc;
		num_proc = inum_proc;
	}
	void initializeBlockInfo();
	//void readFile(const MPI_Comm& comm, const std::string& filename);

	//void outPutFile(const std::string& filename);

	

	//Make mirrors for those who are not at the symm boundaries.
	void readSymmFiles();
	void selectReadSymmFile();


	void writeFiles();
	void performCommunication();
	void outPutBlockInfo()
	{
		int mtag = 0;
		for (auto iter = bis.begin(); iter != bis.end(); iter++)
		{
			iter->setDim(d_dim);
			iter->write_info(d_filenames[mtag], mtag, cur_proc);
			mtag++;
		}
	}
	void make_movement(int index, double dx, double dy, double dz)
	{
		auto &cur_bi = bis[index];
		cur_bi.move_by(dx, dy, dz);
	}


	void readAll()
	{
		readFiles();
		readForTioga();
		MPI_Barrier(MPI_COMM_WORLD);
		
		readForSamrai();
		std::cout<<"Entering read for unstruct\n";
		readForUnstruct();
		for(int i = 0;i<d_nmesh;i++)
		{
			UBs2[i].gNwbc = bis[i].gNwbc;
			UBs2[i].gNwbcCoord=bis[i].gWxyz;
		}
		
	}

private:
	void selectReadForUnstruct();
	void readForTioga();
	void readForSamrai();
	void readForUnstruct();
	void readFiles();
	void selectReadFile();
	void readOversetBoundary();
};
