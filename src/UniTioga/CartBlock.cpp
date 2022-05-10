#include "codetypes.h"
#include "CartBlock.h"
#include "CartGrid.h"
#include <assert.h>
void deallocateLinkList(DONORLIST *temp);
void deallocateLinkList2(INTEGERLIST *temp);
void deallocateLinkList3(INTEGERLIST2 *temp);
void deallocateLinkList4(INTERPLIST2 *temp);
void insertInList(DONORLIST **donorList, DONORLIST *temp1);
void get_amr_index_xyz(int dim, int nq, int i, int j, int k,
					   int pBasis,
					   int nX, int nY, int nZ,
					   int nf,
					   double xlo[3], double dx[3],
					   double qnodes[],
					   int *index, double *xyz);
void amr_index_to_ijklmn(int dim, int pBasis, int nX, int nY, int nZ, int nf, int nq,
						 int index, int *ijklmn);
int checkHoleMap(int dim, double *x, int *nx, int *sam, double *extents);
CartBlock::CartBlock()
{
}

CartBlock::CartBlock(int dim_in)
	: d_dim(dim_in)
{
}
CartBlock::~CartBlock()
{
	clearLists();
}

void CartBlock::preprocess(CartGrid *cg)
{
	
	for (int n = 0; n < d_dim; n++)
	{
		xlo[n] = cg->xlo[d_dim * global_id + n];
		dx[n] = cg->dx[d_dim * global_id + n];
		ncells[n] = cg->ncells[d_dim * global_id + n];
	}
	
	pdegree = cg->porder[global_id];
	pd = 1;
	for (int n = 0; n < d_dim; n++)
	{
		pd = pd * (pdegree + 1);
	}
	nf = cg->nf;
	myid = cg->myid;
	qstride = cg->qstride;
	donor_frac = cg->donor_frac;
	qnode = cg->qnode;
	d1 = ncells[0];
	d2 = ncells[0] * ncells[1];
	dd = d2;
	if (d_dim == 3)
	{
		dd = dd * ncells[2];
	}
	ddnf = 1; //Is the dd with nf
	for (int i = 0; i < d_dim; i++)
	{
		ddnf = ddnf * (ncells[i] + 2 * nf);
	}
	for(int i=0;i<ddnf;i++)
	{
		ibl[i] = 1;
	}
	
	ndof = dd * pd; //Product of dd(stride of ijk) and pd(stride of lmn)
}

void CartBlock::getInterpolatedData(int *nints, int *nreals, int **intData, double **realData, int nvar)
{
}

void CartBlock::update(double *qval, int index, int nq)
{
}

void CartBlock::getCancellationData(int *cancelledData, int *ncancel)
{
	int m = 0;
	int idof = -1;
	*ncancel = 0;
	for (int k = 0; k < ncells[2]; k++)
		for (int j = 0; j < ncells[1]; j++)
			for (int i = 0; i < ncells[0]; i++)
				for (int p = 0; p < pd; p++)
				{
					idof++;
					if (donorList[idof] != NULL)
					{
						DONORLIST *temp = donorList[idof];
						while (temp != NULL)
						{
							if (temp->cancel == 1)
							{
								(*ncancel)++;
								cancelledData[m++] = temp->donorData[0];
								cancelledData[m++] = 1;
								cancelledData[m++] = temp->donorData[2];
								cancelledData[m++] = temp->donorData[3];
							}
							temp = temp->next;
						}
					}
				}
}

void CartBlock::processDonors(HOLEMAP *holemap, int nmesh)
{
	if (d_dim == 2)
	{
		processDonors2D(holemap, nmesh);
	}
	else if (d_dim == 3)
	{
		processDonors3D(holemap, nmesh);
	}
}

void CartBlock::processDonors2D(HOLEMAP *holemap, int nmesh)
{
	int holeFlag;
	int *iflag = new int[nmesh];
	int *index = new int[pd];
	double *xtmp = new double[d_dim * pd];
	int ibcount = -1;
	int idof = -1;
	int ploc = pdegree * (pdegree + 1) / 2;
	for (int j = 0; j < ncells[1]; j++)
		for (int i = 0; i < ncells[0]; i++)
		{
			ibcount++;
			get_amr_index_xyz(d_dim, qstride, i, j, 0, pdegree, ncells[0], ncells[1], 0, nf, xlo, dx, &qnode[ploc], index, xtmp);
			holeFlag = 1;
			idof = ibcount * pd - 1;
			for (int p = 0; p < pd and holeFlag; p++)
			{
				idof++;
				if (donorList[idof] == NULL) //There was no valid donor
				{
					for (int h = 0; h < nmesh; h++)
					{
						if (holemap[h].existWall)
						{
							if (checkHoleMap(d_dim, &xtmp[d_dim * p], holemap[h].nx, holemap[h].sam, holemap[h].extents))
							{
								int ibindex = (j + nf) * (ncells[0] + 2 * nf) + (i + nf);
								ibl[ibindex] = 0;
								holeFlag = 0;
								break; //Inside any one wall then break
							}
						}
					}
				}
				else //There is valid donor
				{
					DONORLIST *temp = donorList[idof];
					for (int h = 0; h < nmesh; h++)
					{
						iflag[h] = 0;
					}
					while (temp != NULL)
					{
						int meshtag_donor = temp->donorData[1];
						iflag[meshtag_donor - BASE] = 1;
						temp = temp->next;
					}
					for (int h = 0; h < nmesh; h++)
					{
						if (holemap[h].existWall)
						{
							if (!iflag[h]) //Exist wall and have no donor from it
							{
								if (checkHoleMap(d_dim, &(xtmp[d_dim * p]), holemap[h].nx, holemap[h].sam, holemap[h].extents))
								{
									int ibindex = (j + nf) * (ncells[0] + 2 * nf) + (i + nf);
									ibl[ibindex] = 0;
									holeFlag = 0;
									break;
								}
							}
						}
					}
				}
			}
		}
	//Now that all cube with amr points inside wall is tagged
	ibcount = -1;
	idof = -1;
	for (int j = 0; j < ncells[1]; j++)
		for (int i = 0; i < ncells[0]; i++)
		{
			ibcount++;
			int ibindex = (nf + j) * (ncells[0] + 2 * nf) + (nf + i);
			get_amr_index_xyz(d_dim, qstride, i, j, 0, pdegree, ncells[0], ncells[1], 0, nf, xlo, dx, &qnode[ploc], index, xtmp);
			if (ibl[ibindex] == 0)
			{
				idof = ibcount * pd - 1;
				for (int p = 0; p < pd; p++)
				{
					idof++;
					if (donorList[idof] != NULL)
					{
						DONORLIST *temp = donorList[idof];
						while (temp != NULL)
						{
							temp->cancel = 1;
							temp = temp->next;
						}
					}
				}
			}
			else
			{
				int temp_count = 0;
				idof = ibcount * pd - 1;
				for (int p = 0; p < pd; p++)
				{
					idof++;
					if (donorList[idof] != NULL)
					{
						DONORLIST *temp = donorList[idof];
						while (temp != NULL)
						{
							if (temp->donorRes < BIGVALUE)
							{
								temp_count++;
								break;
							}
							temp = temp->next;
						}
					}
				}
				if (temp_count == pd)
				{
					ibl[ibindex] = -1;
				}
				else
				{
					//std::cout<<temp_count<<'\n';
					idof = ibcount * pd - 1;
					for (int p = 0; p < pd; p++)
					{
						idof++;
						if (donorList[idof] != NULL)
						{
							DONORLIST *temp = donorList[idof];
							while (temp != NULL)
							{
								temp->cancel = 1;
								temp = temp->next;
							}
						}
					}
				}
			}
		}
	if (iflag)
		TIOGA_FREE(iflag);
	if (xtmp)
		TIOGA_FREE(xtmp);
	if (index)
		TIOGA_FREE(index);
}

void CartBlock::processDonors3D(HOLEMAP *holemap, int nmesh)
{
	int holeFlag;
	int *iflag = new int[nmesh];
	int *index = new int[pd];
	double *xtmp = new double[d_dim * pd];
	int ibcount = -1;
	int idof = -1;
	int ploc = pdegree * (pdegree + 1) / 2;
	for (int k = 0; k < ncells[2]; k++)
		for (int j = 0; j < ncells[1]; j++)
			for (int i = 0; i < ncells[0]; i++)
			{
				ibcount++;
				get_amr_index_xyz(d_dim, qstride, i, j, k, pdegree, ncells[0], ncells[1], ncells[2], nf, xlo, dx, &qnode[ploc], index, xtmp);
				holeFlag = 1;
				idof = ibcount * pd - 1; //Beginning of the amr index at (i,j,k)
				for (int p = 0; p < pd and holeFlag; p++)
				{
					idof++;
					if (donorList[idof] == NULL)
					{
						for (int h = 0; h < nmesh; h++)
						{
							if (holemap[h].existWall)
							{
								if (checkHoleMap(d_dim, &xtmp[d_dim * p], holemap[h].nx, holemap[h].sam, holemap[h].extents))
								{
									int ibindex = (k + nf) * (ncells[1] + 2 * nf) * (ncells[0] + 2 * nf) + (j + nf) * (ncells[0] + 2 * nf) + (i + nf);
									ibl[ibindex] = 0;
									holeFlag = 0;
									break;
								}
							}
						}
					}
					else
					{
						DONORLIST *temp = donorList[idof];
						while (temp != NULL)
						{
							int meshtag_donor = temp->donorData[1];
							iflag[meshtag_donor - BASE] = 1;
							temp = temp->next;
						}
						for (int h = 0; h < nmesh; h++)
						{
							if (holemap[h].existWall)
							{
								if (!iflag[h])
								{
									if (checkHoleMap(d_dim, &xtmp[d_dim * p], holemap[h].nx, holemap[h].sam, holemap[h].extents))
									{
										int ibindex = (k + nf) * (ncells[1] + 2 * nf) * (ncells[0] + 2 * nf) + (j + nf) * (ncells[0] + 2 * nf) + (i + nf);
										ibl[ibindex] = 0;
										holeFlag = 0;
										break;
									}
								}
							}
						}
					}
				}
			}
	ibcount = -1;
	idof = -1;
	for (int k = 0; k < ncells[2]; k++)
		for (int j = 0; j < ncells[1]; j++)
			for (int i = 0; i < ncells[0]; i++)
			{
				ibcount++;
				int ibindex = (nf + k) * (ncells[1] + 2 * nf) * (ncells[0] + 2 * nf) + (nf + j) * (ncells[0] + 2 * nf) + (nf + i);
				get_amr_index_xyz(d_dim, qstride, i, j, k, pdegree, ncells[0], ncells[1], ncells[2], nf, xlo, dx, &qnode[ploc], index, xtmp);
				if (ibl[ibindex] == 0)
				{
					idof = ibcount * pd - 1;
					for (int p = 0; p < pd; p++)
					{
						idof++;
						if (donorList[idof] != NULL)
						{
							DONORLIST *temp = donorList[idof];
							while (temp != NULL)
							{
								temp->cancel = 1;
								temp = temp->next;
							}
						}
					}
				}
				else
				{
					int temp_count = 0;
					idof = ibcount * pd - 1;
					for (int p = 0; p < pd; p++)
					{
						idof++;
						if (donorList[idof] != NULL)
						{
							DONORLIST *temp = donorList[idof];
							while (temp != NULL)
							{
								if (temp->donorRes < BIGVALUE)
								{
									temp_count++;
									break;
								}
								temp = temp->next;
							}
						}
					}
					if (temp_count == pd)
					{
						ibl[ibindex] = -1;
					}
					else
					{
						idof = ibcount * pd - 1;
						for (int p = 0; p < pd; p++)
						{
							idof++;
							if (donorList[idof] != NULL)
							{
								DONORLIST *temp = donorList[idof];
								while (temp != NULL)
								{
									temp->cancel = 1;
									temp = temp->next;
								}
							}
						}
					}
				}
			}

	if (iflag)
		TIOGA_FREE(iflag);
	if (xtmp)
		TIOGA_FREE(xtmp);
	if (index)
		TIOGA_FREE(index);
}
#include <iostream>
void CartBlock::insertInDonorList(int senderid, int index, int meshtagdonor, int remoteid, int remoteblockid, double cellRes)
{
	//int maximum_pointid = ncells[0]
	DONORLIST *temp1 = new DONORLIST[1];
	int ijklmn[6];
	int pointid;
	amr_index_to_ijklmn(d_dim, pdegree, ncells[0], ncells[1], ncells[2], nf, qstride, index, ijklmn);
	//The mapping relationship from ijklmn to pointid inside (ndof) is specified here
	if (d_dim == 2)
	{
		pointid = (ijklmn[1] * d1 + ijklmn[0]) * pd + ijklmn[3] * (pdegree + 1) + ijklmn[2];
	}
	else if (d_dim == 3)
	{
		pointid = (ijklmn[2] * d2 + ijklmn[1] * d1 + ijklmn[0]) * pd + ijklmn[5] * (pdegree + 1) * (pdegree + 1) + ijklmn[4] * (pdegree + 1) + ijklmn[3];
	}
	temp1->donorData[0] = senderid;
	temp1->donorData[1] = meshtagdonor;
	temp1->donorData[2] = remoteid;
	temp1->donorData[3] = remoteblockid;
	temp1->donorRes = cellRes;
	temp1->cancel = 0;
	insertInList(&(donorList[pointid]), temp1);
}

void CartBlock::insertInInterpList(int procid, int remoteid, int remoteblockid, double *xtmp)
{

	int ix[3];
	double rst[3];
	interpList.push_back(new INTERPLIST2[1]);
	// if (myid == 1)
	// {
	// 	std::cout << interpList.size() << '\n';
	// }
	interpList.back()->inode = NULL;
	interpList.back()->weights = NULL;
	interpList.back()->receptorInfo[0] = procid;
	interpList.back()->receptorInfo[1] = remoteid;

	for (int n = 0; n < d_dim; n++)
	{
		ix[n] = (xtmp[n] - xlo[n]) / dx[n];
		rst[n] = (xtmp[n] - xlo[n] - ix[n] * dx[n]) / dx[n];
		if (ix[n] == ncells[n])
		{
			if (fabs(rst[n]) < TOL)
			{
				ix[n]--;
				rst[n] = (xtmp[n] - xlo[n] - ix[n] * dx[n]) / dx[n];
			}
		}
	}
	if (d_dim == 2)
	{
		interpList.back()->nweights = (pdegree + 1) * (pdegree + 1);
	}
	else if (d_dim == 3)
	{
		interpList.back()->nweights = (pdegree + 1) * (pdegree + 1) * (pdegree + 1);
	}
	interpList.back()->weights = new double[interpList.back()->nweights];
	interpList.back()->inode = new int[d_dim];
	for (int i = 0; i < d_dim; i++)
	{
		interpList.back()->inode[i] = ix[i];
	}
	//TODO::Find a donor_frac here
	return;
	//donor_frac(&pdegree, rst, &(listptr->nweights), listptr->weights);
}

void CartBlock::writeCellFile(const std::string &filename)
{
	int w_nnodes = (ncells[0] + 1) * (ncells[1] + 1);
	int w_ncells = ncells[0] * ncells[1];
	if (d_dim == 3)
	{
		w_nnodes = w_nnodes * (ncells[2] + 1);
		w_ncells = w_ncells * ncells[2];
	}
	std::string local_filename = filename + std::to_string(myid)+"_"+std::to_string(local_id);
	std::ofstream fout;
	fout.open(local_filename);
	if (!fout.is_open())
	{
		fout << "Open file failure\n";
		return;
	}
	fout << "TITLE = \"CART\"\n";
	fout << "VARIABLES=\"X\",\"Y\",";
	if (d_dim == 3)
	{
		fout << "\"Z\",";
	}
	fout << "\"IBLANK_CELL\"\n";
	fout << "ZONE T=\"VOL_MIXED\",N=" << w_nnodes << " E=" << w_ncells << " ";
	if (d_dim == 2)
	{
		fout << "ET=QUADRILATERAL,"
			 << "F=FEBLOCK\n";
		fout << "VARLOCATION = (1=NODAL,2=NODAL,3=CELLCENTERED)\n";
	}
	else if (d_dim == 3)
	{
		fout << "ET=BRICK,"
			 << "F=FEBLOCK\n";
		fout << "VARLOCATION = (1=NODAL,2=NODAL,3=NODAL,4=CELLCENTERED)\n";
	}
	if (d_dim == 2)
	{
		for (int j = 0; j < ncells[1] + 1; j++)
		{
			for (int i = 0; i < ncells[0] + 1; i++)
			{
				fout << xlo[0] + dx[0] * i << "\n";
			}
		}
		for (int j = 0; j < ncells[1] + 1; j++)
		{
			for (int i = 0; i < ncells[0] + 1; i++)
			{
				fout << xlo[1] + dx[1] * j << "\n";
			}
		}
		for (int j = 0; j < ncells[1]; j++)
		{
			for (int i = 0; i < ncells[0]; i++)
			{
				int ibindex = (j + nf) * (ncells[0] + 2 * nf) + (i + nf);
				fout << ibl[ibindex] << "\n";
			}
		}
		int id = 0;
		int dd1 = ncells[0] + 1;

		for (int j = 0; j < ncells[1]; j++)
		{
			for (int i = 0; i < ncells[0]; i++)
			{
				id = j * dd1 + (i + 1);
				fout << id << " "
					 << id + 1 << " "
					 << id + 1 + dd1 << " "
					 << id + dd1 << "\n";
			}
		}
	}
	else if (d_dim == 3)
	{
		for (int k = 0; k < ncells[2] + 1; k++)
		{
			for (int j = 0; j < ncells[1] + 1; j++)
			{
				for (int i = 0; i < ncells[0] + 1; i++)
				{
					fout << xlo[0] + dx[0] * i << "\n";
				}
			}
		}
		for (int k = 0; k < ncells[2] + 1; k++)
		{
			for (int j = 0; j < ncells[1] + 1; j++)
			{
				for (int i = 0; i < ncells[0] + 1; i++)
				{
					fout << xlo[1] + dx[1] * j << "\n";
				}
			}
		}
		for (int k = 0; k < ncells[2] + 1; k++)
		{
			for (int j = 0; j < ncells[1] + 1; j++)
			{
				for (int i = 0; i < ncells[0] + 1; i++)
				{
					fout << xlo[2] + dx[2] * k << "\n";
				}
			}
		}
		for (int k = 0; k < ncells[2]; k++)
		{
			for (int j = 0; j < ncells[1]; j++)
			{
				for (int i = 0; i < ncells[0]; i++)
				{
					int ibindex = (k + nf) * (ncells[1] + 2 * nf) * (ncells[0] + 2 * nf) + (j + nf) * (ncells[0] + 2 * nf) + (i + nf);
					fout << ibl[ibindex] << "\n";
				}
			}
		}
		int id = 0;
		int dd1 = ncells[0] + 1;
		int dd2 = dd1 * (ncells[1] + 1);
		for (int k = 0; k < ncells[2]; k++)
		{
			for (int j = 0; j < ncells[1]; j++)
			{
				for (int i = 0; i < ncells[0]; i++)
				{
					id = k * dd2 + j * dd1 + (i + 1);
					fout << id << " "
						 << id + 1 << " "
						 << id + 1 + dd1 << " "
						 << id + dd1 << " "
						 << id + dd2 << " "
						 << id + 1 + dd2 << " "
						 << id + 1 + dd1 + dd2 << " "
						 << id + dd1 + dd2 << "\n";
				}
			}
		}
	}
	fout.close();
}

void CartBlock::clearLists()
{
	if (donorList)
	{
		for (int i = 0; i < ndof; i++)
		{
			deallocateLinkList(donorList[i]);
			donorList[i] = NULL;
		}
		TIOGA_FREE(donorList);
	}
	interpList.clear();
}

void CartBlock::initializeLists()
{
	donorList = new DONORLIST *[ndof];
	for (int i = 0; i < ndof; i++)
	{
		donorList[i] = NULL;
	}
	std::cout << ndof << '\n';
}