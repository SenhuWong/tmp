#include "CobaltReader.h"
#include "metis.h"
#include "../mathTools/LocalGlobalMap.h"
#include "../UniTioga/tioga.h"
#include "../toolBox/edge3d_int.h"
#include "../toolBox/cell3d_int.h"
#include <vector>
#include <regex>
//#include<io.h>
#include <sys/io.h>
#include <set>
#include <iostream>
#include <fstream>
#include <numeric>
#include <algorithm>
#include "hdf5.h"
#include "H5FDmpi.h"
#define UNSTRUCT_DEBUG
void CobaltReader::initializeBlockInfo()
{
	d_nmesh = d_filenames.size();

	TiogaFeeder::initialize();
	SamraiFeeder::initialize();
	UnstructFeeder::initialize();
	std::cout << "End of initialize BlockInfo\n";
}

void CobaltReader::readFiles()
{
	initializeBlockInfo();
	if (cur_proc == 0)
	{
		std::cout << "Entering ReadFile()\n";
	}
	// Assume that the number of input files is less than that of processors.

	for (int m = 0; m < bis.size(); m++)
	{
		auto &cur_bi = bis[m];
		int to_throw_away;

		if (cur_proc == m)
		{
			std::string filename = d_filenames[cur_proc];
			std::ifstream fin;
			fin.open(filename);
			if (!fin.is_open())
			{
				throw std::runtime_error("Open file Failure\n");
			}
			int temp_nnodes;
			int temp_nedges;
			int temp_ncells;
			fin >> to_throw_away >> to_throw_away >> to_throw_away;
			fin >> temp_nnodes >> temp_nedges >> temp_ncells >> to_throw_away >> to_throw_away;

			double *xyzs = new double[d_dim * temp_nnodes];

			for (int i = 0; i < temp_nnodes; i++)
			{
				for (int j = 0; j < d_dim; j++)
				{
					fin >> xyzs[d_dim * i + j];
				}
			}
			int *edge2Cell = new int[2 * temp_nedges];
			int temp_counts;

			int *temp_nodes3D = NULL;
			int temp_nodes[4]; 
			Brick4Cobalt3D *temp_cells3D = NULL;
			Brick4Cobalt2D *temp_cells2D = NULL;

			if (d_dim == 2)
			{
				temp_cells2D = new Brick4Cobalt2D[temp_ncells];
				for (int i = 0; i < temp_nedges; i++)
				{
					fin >> temp_counts;
					for (int j = 0; j < temp_counts; j++)
					{
						fin >> temp_nodes[j];
					}
					fin >> edge2Cell[2 * i] >> edge2Cell[2 * i + 1];

					temp_cells2D[edge2Cell[2 * i] - cellBase].insert(temp_nodes, (double*)&edge2Cell[2*i], 2);
					
					if (edge2Cell[2 * i + 1] > 0)
					{
						temp_cells2D[edge2Cell[2 * i + 1] - cellBase].insert(temp_nodes,(double*)&edge2Cell[2*i+1], 2);
					}
				}
			}
			else if (d_dim == 3)
			{
				temp_nodes3D = new int[4 * temp_nedges];//4 times the num of edges
				temp_cells3D = new Brick4Cobalt3D[temp_ncells];
				for (int i = 0; i < temp_nedges; i++)
				{
					int indx = 4*i;

					fin >> temp_counts;
					for (int j = 0; j < temp_counts; j++)
					{
						fin >> temp_nodes3D[indx + j];
					}
					fin >> edge2Cell[2 * i] >> edge2Cell[2 * i + 1];
					temp_cells3D[edge2Cell[2 * i] - cellBase].insert(&temp_nodes3D[indx],1, xyzs, temp_counts);
					if (edge2Cell[2 * i + 1] > 0)
					{
						temp_cells3D[edge2Cell[2 * i + 1] - cellBase].insert(&temp_nodes3D[indx],-1, xyzs, temp_counts);
					}
				}
				for(int i = 0;i<temp_ncells;i++)
				{
					temp_cells3D[i].reOrdering();
					temp_cells3D[i].checkRepeat();
				}
				delete[] temp_nodes3D;
			}
			fin.close();
			//std::cin.get();
			// Splitting begin
			idx_t num_elements = temp_ncells;
			idx_t num_nodes = temp_nnodes;
			idx_t index_offset = 0;
			std::vector<idx_t> elemPtr;
			std::vector<idx_t> elemInd;
			int *int_iter;
			for (int i = 0; i < temp_ncells; i++)
			{
				elemPtr.push_back(index_offset);
				if (d_dim == 2)
				{
					index_offset += static_cast<idx_t>(temp_cells2D[i].sizeIs());
					int_iter = temp_cells2D[i].begin();
					for (int j = 0; j < temp_cells2D[i].sizeIs(); j++)
					{
						temp_counts = *int_iter;
						// temp_counts should start from 1;
						elemInd.push_back(static_cast<idx_t>(temp_counts - 1)); // Splitting index start from 0
						int_iter++;
					}
				}
				else if (d_dim == 3)
				{
					index_offset += static_cast<idx_t>(temp_cells3D[i].sizeIs());
					int_iter = temp_cells3D[i].begin();
					for (int j = 0; j < temp_cells3D[i].sizeIs(); j++)
					{
						temp_counts = *int_iter;
						// temp_counts should start from 1;
						elemInd.push_back(static_cast<idx_t>(temp_counts - 1)); // Splitting index start from 0
						int_iter++;
					}
				}
			}
			elemPtr.push_back(static_cast<idx_t>(index_offset));
			// if(cur_proc==0)
			// {
			// 	int counter = 0;
			// 	for(auto iter = elemPtr.begin();iter!=elemPtr.end();iter++)
			// 	{
			// 		std::cout<<*iter<<'\t';
			// 		counter++;
			// 		if(counter%6==5)
			// 		{
			// 			std::cout<<"\n";
			// 		}
			// 	}
			// 	counter = 0;
			// 	for(auto iter = elemInd.begin();iter!=elemInd.end();iter++)
			// 	{
			// 		std::cout<<*iter<<'\t';
			// 		counter++;
			// 		if(counter%6==5)
			// 		{
			// 			std::cout<<"\n";
			// 		}
			// 	}
			// 	std::cin.get();
			// }
			std::vector<idx_t> cellPartition(num_elements, -1);
			idx_t NCOMM = 2; // 2D case
			idx_t NPART = static_cast<idx_t>(num_proc);
			idx_t objVal = 0;
			std::vector<idx_t> nodePartition(num_nodes, 0);
			idx_t options[METIS_NOPTIONS];
			METIS_SetDefaultOptions(options);
			options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;
			options[METIS_OPTION_CONTIG] = 1;
			options[METIS_OPTION_NUMBERING] = 0; // Start from 0;
			std::cout << "Calling Mestis part\n";
			if (num_proc > 1)
			{
				METIS_PartMeshDual(&num_elements, &num_nodes, elemPtr.data(), elemInd.data(), NULL, NULL, &NCOMM, &NPART, NULL, options, &objVal, cellPartition.data(), nodePartition.data());
			}
			else
			{
				for (auto iter = cellPartition.begin(); iter != cellPartition.end(); iter++)
				{
					*iter = cur_proc;
				}
			}
			std::cout << "Called Mestis part\n";

			int *eachCount = new int[num_proc];
			for (int i = 0; i < num_proc; i++)
			{
				eachCount[i] = 0;
			}
			
			for (auto itera = cellPartition.begin(); itera != cellPartition.end(); itera++)
			{
				
				// std::cout << *itera << '\n';
				++eachCount[*(itera)];
			}
			
			MPI_Request *send_request = new MPI_Request[num_proc - 1];
			MPI_Status *send_status = new MPI_Status[num_proc - 1];
			
			std::cout << "Sending eachCount\n";
			for (int i = 0; i < num_proc; i++)
			{
				if (i != cur_proc)
				{
					std::cout << "eachCount " << eachCount[i] << '\n';
					MPI_Ssend(&(eachCount[i]), 1, MPI_INT, i, 0, MPI_COMM_WORLD);
				}
			}

			cur_bi.ncellDistributed = eachCount[cur_proc];
			cur_bi.cellDistributed = new int[eachCount[cur_proc]];
			
			delete[] xyzs;
			delete[] edge2Cell;

			delete[] eachCount;
			
			
			// outout a cell file so that we don't need to reinsert to get all the cellsS
			if (true)
			{
				std::string ofilename = "_cell" + filename + ".h5";

				hid_t file_id = H5Fcreate(ofilename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
				
				int eachCellSize = 1 + 1 + (d_dim == 2 ? 4 : 8);
				hsize_t dim[2] = {temp_ncells, eachCellSize};
				hid_t dataspace_id = H5Screate_simple(2, dim, NULL);
				std::string dataset_name = "cell_distribution";
				int *buffer = new int[eachCellSize];
				hid_t dset_id = H5Dcreate(file_id, dataset_name.c_str(), H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

				hsize_t count[2] = {1, eachCellSize};
				hsize_t offset[2] = {0, 0};
				hsize_t stride[2] = {1, 1};
				hsize_t block[2] = {1, 1};
				
				hid_t memspace_id = H5Screate_simple(2, count, NULL);

				if (d_dim == 2)
				{
					for (int i = 0; i < temp_ncells; i++)
					{
						offset[0] = i;
						hid_t status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, stride, count, block);

						buffer[0] = temp_cells2D[i].sizeIs();
						buffer[1] = cellPartition[i];
						int_iter = temp_cells2D[i].begin();
						for (int j = 0; j < temp_cells2D[i].sizeIs(); j++)
						{
							buffer[2 + j] = (*int_iter);
							//std::cout<<buffer[2+j]<<" ";
							++int_iter;
						}
						//std::cout<<'\n';
						//std::cin.get();
						for (int j = temp_cells2D[i].sizeIs(); j < 4; j++)
						{
							buffer[2 + j] = -1;
						}
						status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT, buffer);
						//std::cout<<"Status: "<<status<<'\n';
					}
				}
				else if (d_dim == 3)
				{
					for (int i = 0; i < temp_ncells; i++)
					{
						offset[0] = i;
						hid_t status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, stride, count, block);
						buffer[0] = temp_cells3D[i].sizeIs();
						buffer[1] = cellPartition[i];
						int_iter = temp_cells3D[i].begin();
						for (int j = 0; j < temp_cells3D[i].sizeIs(); j++)
						{
							buffer[2 + j] = (*int_iter);
							++int_iter;
						}
						for (int j = temp_cells3D[i].sizeIs(); j < 8; j++)
						{
							buffer[2 + j] = -1;
						}
						status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT, buffer);
					}
				}
				H5Sclose(memspace_id);
				H5Dclose(dset_id);
				H5Sclose(dataspace_id);
				H5Fclose(file_id);
				delete[] buffer;
			}
			if (temp_cells2D)
				delete[] temp_cells2D;
			if (temp_cells3D)
				delete[] temp_cells3D;
			delete[] send_request;
			delete[] send_status;
		}
		else
		{
			MPI_Request recv_request;
			MPI_Status recv_status;
			MPI_Recv(&(cur_bi.ncellDistributed), 1, MPI_INT, m, 0, MPI_COMM_WORLD, &recv_status);

			int *cellDistributed = new int[cur_bi.ncellDistributed];
			cur_bi.cellDistributed = cellDistributed;
		}
	}
}

void CobaltReader::selectReadForUnstruct()
{
	for (int m = 0; m < bis.size(); m++)
	{
		std::string localFilename = "reader" + std::to_string(cur_proc) + std::to_string(m) + ".ddt";
		std::ofstream recorder;
		recorder.open(localFilename);
		auto &cur_bi = bis[m];
		UnstructBlock2D<2> *cur_ub2D;
		UnstructBlock2D<3> *cur_ub3D;
		if (d_dim == 2)
		{
			cur_ub2D = &(UBs2[m]);
		}
		else if (d_dim == 3)
		{
			cur_ub3D = &(UBs3[m]);
		}
		std::string filename = d_filenames[m];
		std::ifstream finGrid;
		//std::ifstream finCell;

		finGrid.open(filename, std::ios::in);
		if (!finGrid.is_open())
		{
			throw std::runtime_error("Open grid file failure\n");
		}

		int temp_nnodes, temp_nedges, temp_ncells;
		int to_throw_away;
		finGrid >> to_throw_away >> to_throw_away >> to_throw_away;
		finGrid >> temp_nnodes >> temp_nedges >> temp_ncells >> to_throw_away >> to_throw_away;

		double *temp_xyz = new double[temp_nnodes * d_dim];
		int temp_count;
		for (temp_count = 0; temp_count < temp_nnodes; temp_count++)
		{
			for (int j = 0; j < d_dim; j++)
			{
				finGrid >> temp_xyz[d_dim * temp_count + j];
			}
		}

		int *edge2Cell = new int[2 * temp_nedges];
		std::vector<int> edgeNodePtr(temp_nedges + 1, -1);
		std::vector<int> edgeNodeInd;
		edgeNodeInd.reserve(2 * (d_dim - 1) * temp_nedges);
		int temp_nodes[4];
		edgeNodePtr[0] = 0;
		for (int i = 0; i < temp_nedges; i++)
		{
			finGrid >> temp_count;
			edgeNodePtr[i + 1] = edgeNodePtr[i] + temp_count;
			for (int j = 0; j < temp_count; j++)
			{
				finGrid >> temp_nodes[j];
				edgeNodeInd.push_back(temp_nodes[j]);
			}
			finGrid >> edge2Cell[2 * i] >> edge2Cell[2 * i + 1];
		}
		finGrid.close();
		// Extend cells distributed on cur_proc by nfringe.
		int *buffered_flag = new int[temp_ncells];
		int *buffered_flagtmp = new int[temp_ncells];
		for (int i = 0; i < temp_ncells; i++)
		{
			buffered_flag[i] = buffered_flagtmp[i] = -1;
		}
		for (int i = 0; i < cur_bi.ncellDistributed; i++)
		{
			// cellDIstributed starts from 1
			buffered_flag[cur_bi.cellDistributed[i] - 1] = buffered_flagtmp[cur_bi.cellDistributed[i] - 1] = 0;
		}
		
		int nfringe = 2;
		for (int i = 0; i < nfringe; i++)
		{
			for (int j = 0; j < temp_nedges; j++)
			{
				if (edge2Cell[2 * j] < 0 or edge2Cell[2 * j + 1] < 0)
				{
					continue;
				}
				bool leftInside = buffered_flag[edge2Cell[2 * j] - cellBase] >= 0;
				bool rightInside = buffered_flag[edge2Cell[2 * j + 1] - cellBase] >= 0;
				if (leftInside and !rightInside)
				{
					buffered_flagtmp[edge2Cell[2 * j + 1] - cellBase] = i+1;
				}
				if (!leftInside and rightInside)
				{
					buffered_flagtmp[edge2Cell[2 * j] - cellBase] = i+1;
				}
			}
			for (int j = 0; j < temp_ncells; j++)
			{
				buffered_flag[j] = buffered_flagtmp[j];
			}
		}
		// Find all Edges and Nodes related to the extended localCells;
		std::set<int, std::less<int>> *hereNodes = new std::set<int, std::less<int>>[1];
		std::set<int, std::less<int>> *hereCells = new std::set<int, std::less<int>>[1];

		//Wall nodes needed by SST Model
		std::set<int, std::less<int>> *globalWbc = new std::set<int, std::less<int>>[1];
		int nEdges;
		int nCells;

		int size;
		int inode;
		int *cellRealDistribution = new int[temp_ncells];
		
		std::string curIntFilename = "_cell" + filename + ".h5";
		/*
			Here comes the parallel HDF5
		*/
		
		MPI_Comm comm = MPI_COMM_WORLD;
		MPI_Info info = MPI_INFO_NULL;
		hid_t plist_id=H5Pcreate(H5P_FILE_ACCESS);
		H5Pset_fapl_mpio(plist_id,comm,info);
		hid_t file_id = H5Fopen(curIntFilename.c_str(),H5P_DEFAULT,H5FD_MPIO_INDEPENDENT);
		hid_t status = H5Pclose(plist_id);
		std::string curDatasetname = "cell_distribution";
		//H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_INDEPENDENT);
		hid_t dataset_id = H5Dopen(file_id,curDatasetname.c_str(),H5P_DEFAULT);


		hid_t file_dataspace = H5Dget_space(dataset_id);
		int eachCellSize = 1 + 1 + (d_dim == 2 ? 4 : 8);
		int *buffer = new int[eachCellSize];
		hsize_t count[2] = {1, eachCellSize};
		hsize_t offset[2] = {0, 0};
		hsize_t stride[2] = {1, 1};
		hsize_t block[2] = {1, 1};
		hid_t memspace_id = H5Screate_simple(2, count, NULL);
		
		
		for (int i = 0; i < temp_ncells; i++)
		{
			
			offset[0] = i;
			hid_t status = H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, offset, stride, count, block);
			H5Dread(dataset_id, H5T_NATIVE_INT, memspace_id, file_dataspace, H5P_DEFAULT, buffer);
			cellRealDistribution[i] = buffer[1];
			if(buffered_flag[i]>=0)
			{
				for (int j = 0; j < buffer[0]; j++)
				{
					hereNodes->insert(buffer[2+j]);
				}
				hereCells->insert(i + cellBase);	
			}

		}
		

		LocalGlobalMap hereNodeMapping(hereNodes->size(), hereNodes->begin(), hereNodes->end());
		LocalGlobalMap hereCellMapping(hereCells->size(), hereCells->begin(), hereCells->end());
		// Allocate Point3d.
		GeomElements::point3d<2> *pt2d;
		GeomElements::point3d<3> *pt3d;
		if (d_dim == 2)
		{
			pt2d = new GeomElements::point3d<2>[hereNodes->size()];
			int tmp_count = 0;
			for (auto iter = hereNodes->begin(); iter != hereNodes->end(); iter++)
			{
				pt2d[tmp_count++] = GeomElements::vector3d<2, double>(temp_xyz[2 * ((*iter) - nodeBase)], temp_xyz[2 * ((*iter) - nodeBase) + 1]);
			}
		}
		else if (d_dim == 3)
		{
		}
		// Allocate Cell3d
		GeomElements::cell3d<2> *cel2d;
		GeomElements::cell3d<3> *cel3d;
		// Sort out the distribution of communication cells on each Proc
		std::vector<int> nEachProcBufferCell(num_proc, 0);
		int nlocalCells = 0;
		// How many communication cells there is on each proc
		for (int i = 0; i < temp_ncells; i++)
		{
			if (buffered_flag[i] == 0)
			{
				++nlocalCells;
			}
			else if (buffered_flag[i] > 0)
			{
				++nEachProcBufferCell[cellRealDistribution[i]];
			}
		}
		int nbufferCells = std::accumulate(nEachProcBufferCell.begin(), nEachProcBufferCell.end(), 0);
		// How many related proc there is on this proc.

		temp_count = 0;

		std::map<int, int> *index2Proc = new std::map<int, int>[1];
		std::map<int, int> *proc2Index = new std::map<int, int>[1];

		for (int i = 0; i < num_proc; i++)
		{
			if (nEachProcBufferCell[i] != 0)
			{

				(*index2Proc)[temp_count] = i;
				(*proc2Index)[i] = temp_count;
				temp_count++;
			}
		}
		int nRelatedProcs = index2Proc->size();

		// Map is formed here, we should give it to BLock later.
		//%TODO::This should be modified later,we don't want to use BlockInfo.
		if (nlocalCells != cur_bi.ncells)
		{
			std::cout << "SOmething wrong\n";
			std::cin.get();
		}
		nlocalCells += nbufferCells;

		if (d_dim == 2)
		{
			cel2d = new GeomElements::cell3d<2>[nlocalCells];
			
			int tmp_count = 0;
			int tmp_count_scratch = 0;
			for(int i = 0;i<temp_ncells;i++)
			{
				if(buffered_flag[i]>=0)
				{
					offset[0]=i;
					hid_t status = H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, offset, stride, count, block);
					H5Dread(dataset_id, H5T_NATIVE_INT, memspace_id, file_dataspace, H5P_DEFAULT, buffer);
					for (int j = 0; j < buffer[0]; j++)
					{
						//cel2d[tmp_count].push_back(-1);
						cel2d[tmp_count].push_back(hereNodeMapping.global2local(buffer[2+j]) - 1);
						//std::cout<<buffer[2+j]<<"->"<<hereNodeMapping.global2local(buffer[2+j])<<'\n';
						//std::cin.get();
					}
					
					if (buffered_flag[i] > 0) // Is a buffered cell,we could compute its local id and globalid(stored in remote id for now) and proc for sending
					{
						int localid = hereCellMapping.global2local(i + 1) - 1; // Starts from 0;
						cur_ub2D->recvCommCells().push_back(CommunicationCell(localid, proc2Index->find(cellRealDistribution[i])->second, i + 1,buffered_flag[i]));
						// recv_cells[tmp_count_scratch].setCommCell(localid,proc2Index->find(cellRealDistribution[i])->second,i);
						tmp_count_scratch++;
					}
					tmp_count++;
				}
			}
		}
		status = H5Sclose(memspace_id);
		status = H5Sclose(file_dataspace);
		status = H5Dclose(dataset_id);
		status = H5Fclose(file_id);
		delete[] buffer;

		recorder << "nRelatedProc of this block is " << nRelatedProcs << '\n';

		// Now for Cell communication,we need
		// 1.Proc associated with current proc and its mapping  (done)
		// 2.How many CommunicationCell there is on each proc
		// 3.Send the receive request and make those send request registered on remote side.

		// 1.make relevant proc mapping

		// 3.Send the receive request by 2 steps.
		// First let each other know how many cells there is.
		std::vector<int> nEachProcSendCells(nRelatedProcs, 0);
		std::vector<int> nEachProcRecvCells(nRelatedProcs, 0);
		for (auto iter = index2Proc->begin(); iter != index2Proc->end(); iter++)
		{
			nEachProcRecvCells[iter->first] = nEachProcBufferCell[iter->second];
			recorder << iter->first << " " << iter->second << " " << nEachProcRecvCells[iter->first] << '\n';
		}

		std::vector<MPI_Request> EachPrcoRequestsSend(nRelatedProcs);
		std::vector<MPI_Request> EachPrcoRequestsRecv(nRelatedProcs);
		std::vector<MPI_Status> EachProcStatusRecv(nRelatedProcs);

		for (int i = 0; i < nRelatedProcs; i++)
		{
			MPI_Irecv(&nEachProcSendCells[i], 1, MPI_INT, index2Proc->find(i)->second, 0, MPI_COMM_WORLD, &(EachPrcoRequestsRecv[i]));
			MPI_Isend(&nEachProcRecvCells[i], 1, MPI_INT, index2Proc->find(i)->second, 0, MPI_COMM_WORLD, &(EachPrcoRequestsSend[i]));
		}

		MPI_Waitall(nRelatedProcs, EachPrcoRequestsRecv.data(), EachProcStatusRecv.data());

		for (int i = 0; i < nRelatedProcs; i++)
		{
			recorder << nEachProcSendCells[i] << '\n';
		}
		// Then Allocate space for sending and recving the communicationCells
		// First Alocate satially contiguous memory for send Recv at EachProc;
		int nlocal_recvCell_buffer = std::accumulate(nEachProcRecvCells.begin(), nEachProcRecvCells.end(), 0);
		int nlocal_sendCell_buffer = std::accumulate(nEachProcSendCells.begin(), nEachProcSendCells.end(), 0);
		int *local_recvCell_buffer = new int[3 * nlocal_recvCell_buffer];
		int *local_sendCell_buffer = new int[3 * nlocal_sendCell_buffer];
		// Calculate the offset for each Proc
		std::vector<int> ProcOffset_s(nRelatedProcs + 1, 0);
		std::vector<int> ProcOffset_r(nRelatedProcs + 1, 0);

		for (int i = 0; i < nRelatedProcs; i++)
		{
			ProcOffset_s[i + 1] = ProcOffset_s[i] + (nEachProcSendCells[i]);
			ProcOffset_r[i + 1] = ProcOffset_r[i] + (nEachProcRecvCells[i]);
		}

		if (d_dim == 2)
		{
			temp_count = 0;
			std::vector<CommunicationCell> &cur_recvCell = cur_ub2D->getRecvCells();
			std::sort(cur_recvCell.begin(), cur_recvCell.end(), rankWithSource); // Rank with the local id of the receive cells
			for (auto iter = cur_recvCell.begin(); iter != cur_recvCell.end(); iter++)
			{
				local_recvCell_buffer[temp_count++] = iter->d_localId;
				local_recvCell_buffer[temp_count++] = iter->d_remoteId;
				local_recvCell_buffer[temp_count++] = iter->d_layer_level;
			}
		}
		for (int i = 0; i < nlocal_recvCell_buffer; i++)
		{
			recorder << local_recvCell_buffer[3 * i] << " " << local_recvCell_buffer[3 * i + 1]<< " " << local_recvCell_buffer[3 * i + 2] << '\n';
		}
		// Sending receive cells to the sender side
		for (int i = 0; i < nRelatedProcs; i++)
		{
			MPI_Irecv(&local_sendCell_buffer[3 * ProcOffset_s[i]], 3 * nEachProcSendCells[i], MPI_INT, index2Proc->find(i)->second, 0, MPI_COMM_WORLD, &(EachPrcoRequestsRecv[i]));
			MPI_Isend(&local_recvCell_buffer[3 * ProcOffset_r[i]], 3 * nEachProcRecvCells[i], MPI_INT, index2Proc->find(i)->second, 0, MPI_COMM_WORLD, &(EachPrcoRequestsSend[i]));
		}

		MPI_Waitall(nRelatedProcs, EachPrcoRequestsRecv.data(), EachProcStatusRecv.data());
		recorder << "---------------------------------------------\n";
		for (int i = 0; i < nlocal_sendCell_buffer; i++)
		{
			recorder << local_sendCell_buffer[3 * i] << " " << local_sendCell_buffer[3 * i + 1] << " " << local_sendCell_buffer[3 * i + 2] << '\n';
		}
		temp_count = 0;
		for (int i = 0; i < nRelatedProcs; i++)
		{
			int proc_id = i;
			for (int j = ProcOffset_s[i]; j < ProcOffset_s[i + 1]; j++)
			{
				int remote_id = local_sendCell_buffer[3*temp_count];
				int global_id = local_sendCell_buffer[3*temp_count+1];
				int local_id = hereCellMapping.global2local(global_id) - 1;
				int layer_level = local_sendCell_buffer[3*temp_count+2];
				if (d_dim == 2)
				{
					cur_ub2D->pushSendCell(CommunicationCell(local_id, proc_id, remote_id,layer_level));
				}
				local_sendCell_buffer[temp_count] = local_id;
				++temp_count;
			}
		}
		for (int i = 0; i < nRelatedProcs; i++)
		{
			MPI_Irecv(&local_recvCell_buffer[ProcOffset_r[i]], nEachProcRecvCells[i], MPI_INT, index2Proc->find(i)->second, 1, MPI_COMM_WORLD, &(EachPrcoRequestsRecv[i]));
			MPI_Isend(&local_sendCell_buffer[ProcOffset_s[i]], nEachProcSendCells[i], MPI_INT, index2Proc->find(i)->second, 1, MPI_COMM_WORLD, &(EachPrcoRequestsSend[i]));
		}
		MPI_Waitall(nRelatedProcs, EachPrcoRequestsRecv.data(), EachProcStatusRecv.data());

		if (d_dim == 2)
		{
			temp_count = 0;
			std::vector<CommunicationCell> &cur_recvCell = cur_ub2D->getRecvCells();
			for (auto iter = cur_recvCell.begin(); iter != cur_recvCell.end(); iter++)
			{
				iter->d_remoteId = local_recvCell_buffer[temp_count++];
			}
			std::vector<CommunicationCell> &cur_sendCell = cur_ub2D->getSendCells();
			std::sort(cur_sendCell.begin(), cur_sendCell.end(), rankWithDestin); // Rank with the local id of the recv cells on the other side.
		}

		// Allocate Edge3d
		GeomElements::edge3d<2> *edg2d;
		GeomElements::edge3d<3> *edg3d;
		std::vector<int> local_edge_inds;
		local_edge_inds.reserve(temp_nedges / num_proc);

		std::vector<int> unFluxedEdgeInd;
		std::vector<int> fluxedEdgeInd;

		std::vector<int> wallBoundEdgeInd;
		std::vector<int> overBoundEdgeInd;
		// std::vector<int> commBoundEdgeInd;
		for (int i = 0; i < temp_nedges; i++)
		{
			// natural Boundary Type just find a minimal between left and right.
			// So if the minimum is a negative, then it is a natural boundary.
			int naturalBoundaryTypes = std::min(edge2Cell[2 * i], edge2Cell[2 * i + 1]);
			int naturalBoundaryOtherSide = std::max(edge2Cell[2 * i], edge2Cell[2 * i + 1]);

			bool naturalBoundaryExists = naturalBoundaryTypes < 0;
			bool withInLocalBufferedArea = false;
			if (naturalBoundaryExists) // Then check if the other side is in buffered area.
			{
				
				// If buffered_flag tag is 0, then it is inside the split of Metis
				// If is greater or equal to 1, then it is inside the buffer region.
				withInLocalBufferedArea = buffered_flag[naturalBoundaryOtherSide - cellBase] >= 0;
				
				if (withInLocalBufferedArea) // The other side is in buffered area, then
				{
					// std::cout << "WithInLocalBufferedArea";
					if (naturalBoundaryTypes == signWall)
					{
						// std::cout << "SignWall found\n";
						wallBoundEdgeInd.push_back(local_edge_inds.size());
					}
					else if (naturalBoundaryTypes == signOver)
					{
						// std::cout << "SignOver found\n";
						overBoundEdgeInd.push_back(local_edge_inds.size());
					}
					else
					{
						throw std::runtime_error("Undefined Boundary Type\n");
					}
					local_edge_inds.push_back(i);
				}
				else // The other side is not in buffered area, leave it alone
				{
				}
			}
			else // Natural Boundary doesn't exists, then all we need to do is check if bothside is in buffered area.
			{

				withInLocalBufferedArea = buffered_flag[naturalBoundaryTypes - cellBase] >= 0 and buffered_flag[naturalBoundaryOtherSide - cellBase] >= 0;
				if (withInLocalBufferedArea) // Still need to check if it is a fluxed edge or unfluxed
				{
					bool isFluxedEdge = buffered_flag[naturalBoundaryTypes - cellBase] == 0 or buffered_flag[naturalBoundaryOtherSide - cellBase] == 0;

					if (isFluxedEdge) // std::cout<<"withInBUfferedLocal found\n";
					{
						local_edge_inds.push_back(i);
					}
					else
					{
						unFluxedEdgeInd.push_back(i);
					}
				}
				else
				{
				}
			}
		}
		int nFluxedEdge = local_edge_inds.size();
		for (int i = 0; i < unFluxedEdgeInd.size(); i++)
		{
			local_edge_inds.push_back(unFluxedEdgeInd[i]);
		}
		

		if (d_dim == 2)
		{
			
			int tmp_count = 0;

			edg2d = new GeomElements::edge3d<2>[local_edge_inds.size()];


			for (auto iter = local_edge_inds.begin(); iter != local_edge_inds.end(); iter++)
			{
				std::cout<<temp_count<<'\n';
				for (int j = edgeNodePtr[*iter]; j < edgeNodePtr[*iter + 1]; j++)
				{
					edg2d[tmp_count].push_back(hereNodeMapping.global2local(edgeNodeInd[j]) - 1);

				}
				std::cout<<temp_count<<'\n';
				int lC_ind = edge2Cell[2 * (*iter)];
				int rC_ind = edge2Cell[2 * (*iter) + 1]; // hereCellMapping.global2local(edge2Cell[2 * (*iter)]); // LeftSide of
				if (lC_ind > 0)
				{
					int lC = hereCellMapping.global2local(lC_ind);
					if (lC > 0)
					{
						edg2d[tmp_count].setLeft(lC - cellBase);
						cel2d[lC - cellBase].push_edge(tmp_count);
						
					}
					else
					{
						throw std::runtime_error("Can't be the case where left is NULL\n");
					}
				}
				if (rC_ind > 0)
				{
					int rC = hereCellMapping.global2local(rC_ind);
					if (rC > 0)
					{
						edg2d[tmp_count].setRight(rC - cellBase);
						cel2d[rC-cellBase].push_edge(tmp_count);
					}
				}
				tmp_count++;	
			}

			for (auto iter = wallBoundEdgeInd.begin(); iter != wallBoundEdgeInd.end(); iter++)
			{
				if (edg2d[*iter].rCInd() > 0)
				{
					throw std::runtime_error("Right side of the wallBoundary is already set\n");
				}
				edg2d[*iter].setRight(GeomElements::edge3d<2>::BoundaryType::WALL);
			}
			for (auto iter = overBoundEdgeInd.begin(); iter != overBoundEdgeInd.end(); iter++)
			{
				if (edg2d[*iter].rCInd() > 0)
				{
					throw std::runtime_error("Right side of the overBoundary is already set\n");
				}
				// std::cout<<"Over found\n";
				edg2d[*iter].setRight(GeomElements::edge3d<2>::BoundaryType::FARFIELD);
			}
			//Here we output the edge 2 cell map of 
		}
		else if (d_dim == 3)
		{
		}

		if (d_dim == 2)
		{
			cur_ub2D->setLocalStructure(pt2d, cel2d, edg2d, hereNodes->size(), hereCells->size(), local_edge_inds.size(), nFluxedEdge);
			cur_ub2D->setLocalCommunication(proc2Index, index2Proc);
		}
		else if (d_dim == 3)
		{
		}
		recorder.close();
	}
}
void CobaltReader::selectReadFile()
{

	if (cur_proc == 0)
	{
		std::cout << "Entering SelectReadFile() " << cur_proc << "\n";
	}

	for (int m = 0; m < bis.size(); m++)
	{
		std::string localFilename = "reader" + std::to_string(cur_proc) + std::to_string(m) + ".ddt";
		std::ofstream recorder;
		recorder.open(localFilename);
		auto &cur_bi = bis[m];

		std::string filename = d_filenames[m];

		Brick4Cobalt2D *hereCells2D = NULL;
		Brick4Cobalt3D *hereCells3D = NULL;
		if (d_dim == 2)
		{
			hereCells2D = new Brick4Cobalt2D[cur_bi.ncellDistributed];
		}
		else if (d_dim == 3)
		{
			hereCells3D = new Brick4Cobalt3D[cur_bi.ncellDistributed];
		}
		std::set<int, std::less<int>> *overNodes = new std::set<int, std::less<int>>[1];
		std::set<int, std::less<int>> *wallNodes = new std::set<int, std::less<int>>[1];
		std::set<int, std::less<int>> *hereNodes = new std::set<int, std::less<int>>[1];
		std::set<int, std::less<int>> *hereWbcNodes = new std::set<int, std::less<int>>[1];
		std::set<int, std::less<int>> *hereObcNodes = new std::set<int, std::less<int>>[1];
		std::ifstream finGrid;
		std::ifstream finCell;

		// First read cellfile to fill the nodes

		// We need to first read how many cells there are.

		finGrid.open(filename, std::ios::in);
		if (!finGrid.is_open())
		{
			throw std::runtime_error("Open grid file failure\n");
		}
		int temp_nnodes, temp_nedges, temp_ncells;
		int to_throw_away;
		finGrid >> to_throw_away >> to_throw_away >> to_throw_away;
		finGrid >> temp_nnodes >> temp_nedges >> temp_ncells >> to_throw_away >> to_throw_away;

		std::string curIntFilename = "_cell" + filename + ".h5";
		/*
			Here comes the parallel HDF5
		*/
		
		MPI_Comm comm = MPI_COMM_WORLD;
		MPI_Info info = MPI_INFO_NULL;

		hid_t acc_tpl1 = H5Pcreate(H5P_FILE_ACCESS);
		hid_t status = H5Pset_fapl_mpio(acc_tpl1,comm,info);

		hid_t file_id = H5Fopen(curIntFilename.c_str(),H5F_ACC_RDONLY,acc_tpl1);
		status = H5Pclose(acc_tpl1);
		
		std::string curDatasetname = "cell_distribution";
		//H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_INDEPENDENT);
		hid_t dataset_id = H5Dopen(file_id,curDatasetname.c_str(),H5P_DEFAULT);
		
		std::cout<<dataset_id<<" at -1\n";

		// hid_t xfer_plist = H5Pcreate(H5P_DATASET_XFER);
		// xfer_plist = H5Pcreate(H5P_DATASET_XFER);
		// status = H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
		int eachCellSize = 1 + 1 + (d_dim == 2 ? 4 : 8);

		hid_t file_dataspace = H5Dget_space(dataset_id);
		hsize_t count[2] = {1, eachCellSize};
		
		
		hid_t memspace_id = H5Screate_simple(2, count, NULL);
		std::cout<<file_dataspace<<" "<<memspace_id<<" at 0\n";

		// hid_t astatus = H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET,NULL,NULL,NULL,NULL);
		std::cout<<file_dataspace<<" "<<memspace_id<<" at 1\n";
		
		int *buffer = new int[eachCellSize];
		hsize_t offset[2] = {0, 0};
		hsize_t stride[2] = {1, 1};
		hsize_t block[2] = {1, 1};
		
		int temp_count = 0;
		int temp_count_scratch = 0;
		if (d_dim == 2)
		{
			for (int i = 0; i < temp_ncells; i++)
			{
				offset[0] = i;
				hid_t status = H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, offset, stride, count, block);
				
				H5Dread(dataset_id, H5T_NATIVE_INT, memspace_id, file_dataspace, H5P_DEFAULT, buffer);
				for (int j = 0; j < eachCellSize; j++)
				{
					recorder << buffer[j] << "\t";
				}
				recorder << '\n';

				if (buffer[1] == cur_proc)
				{
					hereCells2D[temp_count].setSize(buffer[0]);
					for (int j = 0; j < buffer[0]; j++)
					{
						hereNodes->insert(buffer[2+j]);
						hereCells2D[temp_count].setNode(buffer[2 + j], j);
						//std::cout<<buffer[2 + j]<<"\t";
					}
					//std::cout<<'\n';
					//std::cin.get();
					++temp_count;
					cur_bi.cellDistributed[temp_count_scratch++] = i + 1;
				}
				
			}
		}
		else if (d_dim == 3)
		{
			for (int i = 0; i < temp_ncells; i++)
			{
				offset[0] = i;
				hid_t status = H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, offset, stride, count, block);
				// std::cout<<file_dataspace<<" "<<memspace_id<<" at 0\n";

				H5Dread(dataset_id, H5T_NATIVE_INT, memspace_id, file_dataspace,H5P_DEFAULT, buffer);

				if (buffer[1] == cur_proc)
				{
					hereCells3D[temp_count].setSize(buffer[0]);
					for (int j = 0; j < buffer[0]; j++)
					{
						hereNodes->insert(buffer[2+j]);
						hereCells3D[temp_count].setNode(buffer[2 + j], j);
					}
					++temp_count;
					// Make it inside herecells
					cur_bi.cellDistributed[temp_count_scratch++] = i + 1;
				}
				
			}
		}
		// continue;
		status = H5Sclose(memspace_id);
		status = H5Sclose(file_dataspace);
		status = H5Dclose(dataset_id);
		status = H5Fclose(file_id);
		delete[] buffer;
		if (temp_count != cur_bi.ncellDistributed)
		{
			throw std::runtime_error("number of cells not compatible ");
		}

		// Making the mapping from global to local
		LocalGlobalMap hereMapping(hereNodes->size(), hereNodes->begin(), hereNodes->end());
		int nc[10] = {0};
		int size2Type[10] = {-1};
		int ntypes = 0;
		if (d_dim == 2)
		{
			for (int i = 0; i < cur_bi.ncellDistributed; i++)
			{
				nc[hereCells2D[i].sizeIs() - 1]++;
			}
		}
		else if (d_dim == 3)
		{
			for (int i = 0; i < cur_bi.ncellDistributed; i++)
			{
				nc[hereCells3D[i].sizeIs() - 1]++;
			}
		}

		// for (int i = 0; i < 10; i++)
		// {
		// 	std::cout << nc[i] << '\n';
		// }

		for (int i = 0; i < 10; i++)
		{
			if (nc[i] > 0)
			{
				ntypes++;
			}
		}

		cur_bi.ntypes = ntypes;
		cur_bi.vconn = new int *[ntypes];
		cur_bi.ncells = cur_bi.ncellDistributed;
		cur_bi.eachNodeCount = new int[ntypes];
		cur_bi.eachCellCount = new int[ntypes];

		temp_count = 0;
		for (int i = 0; i < 10; i++)
		{
			if (nc[i] > 0)
			{
				size2Type[i] = temp_count;
				cur_bi.eachNodeCount[temp_count] = i + 1;
				cur_bi.eachCellCount[temp_count] = nc[i];
				temp_count++;
			}
		}

		// std::cout << "ntypes are " << ntypes << '\n';

		for (int i = 0; i < ntypes; i++)
		{
			// std::cout<<cur_bi.eachCellCount[i]<<" "<<cur_bi.eachNodeCount[i]<<'\n';
			// std::cout<<cur_bi.vconn[i]<<'\n';
			int size = cur_bi.eachCellCount[i] * cur_bi.eachNodeCount[i];
			std::cout << "size is " << size << '\n';
			cur_bi.vconn[i] = new int[size];
		}

		std::vector<int> eachTypeCount(ntypes, 0);
		if (d_dim == 2)
		{
			for (int i = 0; i < cur_bi.ncellDistributed; i++)
			{
				// std::cout<<hereCells2D[i].sizeIs();
				int curtype = size2Type[hereCells2D[i].sizeIs() - 1];
				// std::cout<<curtype<<" "<<cur_bi.eachNodeCount[0]<<'\n';
				for (int j = 0; j < cur_bi.eachNodeCount[curtype]; j++)
				{
					cur_bi.vconn[curtype][eachTypeCount[curtype] * cur_bi.eachNodeCount[curtype] + j] = hereMapping.global2local(hereCells2D[i].nodeIs(j));
				}
				eachTypeCount[curtype] += 1;
			}
		}

		else if (d_dim == 3)
		{

			for (int i = 0; i < cur_bi.ncellDistributed; i++)
			{
				int curtype = size2Type[hereCells3D[i].sizeIs() - 1];
				for (int j = 0; j < cur_bi.eachNodeCount[curtype]; j++)
				{
					cur_bi.vconn[curtype][eachTypeCount[curtype] * cur_bi.eachNodeCount[curtype] + j] = hereMapping.global2local(hereCells3D[i].nodeIs(j));
				}
				eachTypeCount[curtype] += 1;
			}
		}

		if (hereCells2D)
			delete[] hereCells2D;
		if (hereCells3D)
			delete[] hereCells3D;
		// Now we have got vconn and nc and nv and ntypes all set, it's time to get node informations.
		double *temp_xyz = new double[temp_nnodes * d_dim];
		for (temp_count = 0; temp_count < temp_nnodes; temp_count++)
		{
			for (int j = 0; j < d_dim; j++)
			{
				finGrid >> temp_xyz[d_dim * temp_count + j];
			}
		}

		int edgeNodes[10];
		int edgeCells[2];
		int temp_int = 0; // Used to record ptr

		for (int i = 0; i < temp_nedges; i++)
		{

			finGrid >> temp_count;

			// recorder << temp_count;
			for (int j = 0; j < temp_count; j++)
			{
				finGrid >> edgeNodes[j];
				// recorder << edgeNodes[j] << '\t';
			}
			finGrid >> edgeCells[0] >> edgeCells[1];
			recorder << edgeCells[0] << '\t' << edgeCells[1] << '\n';

			if (edgeCells[1] == signWall)
			{
				for (int j = 0; j < temp_count; j++)
				{
					wallNodes->insert(edgeNodes[j]);
				}
			}
			else if (edgeCells[1] == signOver)
			{

				for (int j = 0; j < temp_count; j++)
				{
					overNodes->insert(edgeNodes[j]);
				}
			}
			else if (edgeCells[1] > 0)
			{
			}
			else
			{
				std::cout << "Pre is " << edgeNodes[0] << " " << edgeNodes[1] << '\n';
				std::cout << "proc and meshtag is " << cur_proc << " " << m + 1 << " " << i << " out of " << temp_nedges << " : " << edgeCells[0] << " " << edgeCells[1] << '\n';
				std::cout << "Something wrong,stop\n";
			}
		}
		finGrid.close();

		int local_nnodes = hereNodes->size();
		std::cout<<"Local nnode is "<<local_nnodes<<'\n';

		double *local_xyz = new double[local_nnodes * d_dim];
		cur_bi.gNwbc = wallNodes->size();
		cur_bi.gNobc = overNodes->size();
		cur_bi.gWbcNode = new int[cur_bi.gNwbc];
		cur_bi.gObcNode = new int[cur_bi.gNobc];
		cur_bi.gWxyz = new double[d_dim * cur_bi.gNwbc];
		cur_bi.gOxyz = new double[d_dim * cur_bi.gNobc];

		std::set<int, std::less<int>>::iterator itera = hereNodes->begin();

		// Now find all the wbc obc and local nodes coordinates
		temp_count = 0;

		for (auto iter = hereNodes->begin(); iter != hereNodes->end(); iter++)
		{
			for (int j = 0; j < d_dim; j++)
			{
				local_xyz[temp_count++] = temp_xyz[d_dim * ((*iter) - 1) + j];
			}
		}

		int temp_count_int = 0;
		temp_count = 0;
		for (auto iter = wallNodes->begin(); iter != wallNodes->end(); iter++)
		{
			cur_bi.gWbcNode[temp_count_int++] = *iter;
			for (int j = 0; j < d_dim; j++)
			{
				cur_bi.gWxyz[temp_count++] = temp_xyz[d_dim * ((*iter) - 1) + j];
			}
		}
		temp_count_int = 0;
		temp_count = 0;
		for (auto iter = overNodes->begin(); iter != overNodes->end(); iter++)
		{
			cur_bi.gObcNode[temp_count_int++] = *iter;
			for (int j = 0; j < d_dim; j++)
			{
				cur_bi.gOxyz[temp_count++] = temp_xyz[d_dim * (*(iter)-1) + j];
			}
		}

		cur_bi.rxyz = local_xyz;

		cur_bi.nnodes = local_nnodes;

		std::insert_iterator<std::set<int, std::less<int>>> hereOBCInserter(*hereObcNodes, hereObcNodes->begin());
		std::insert_iterator<std::set<int, std::less<int>>> hereWBCInserter(*hereWbcNodes, hereWbcNodes->begin());

		std::set_intersection(hereNodes->begin(), hereNodes->end(), overNodes->begin(), overNodes->end(), hereOBCInserter);
		std::set_intersection(hereNodes->begin(), hereNodes->end(), wallNodes->begin(), wallNodes->end(), hereWBCInserter);

		delete[] hereNodes;
		delete[] overNodes;
		delete[] wallNodes;

		int *obc = NULL;
		int *wbc = NULL;
		if (hereObcNodes->size() > 0)
		{
			std::cout << "size of obc on proc " << cur_proc << " meshtag " << m + 1 << " is " << hereObcNodes->size() << '\n';
			obc = new int[hereObcNodes->size()];
		}
		if (hereWbcNodes->size() > 0)
		{
			std::cout << "size of wbc on proc " << cur_proc << " meshtag " << m + 1 << " is " << hereWbcNodes->size() << '\n';
			wbc = new int[hereWbcNodes->size()];
		}

		temp_count = 0;
		for (auto iter = hereObcNodes->begin(); iter != hereObcNodes->end(); iter++)
		{
			obc[temp_count] = hereMapping.global2local(*iter);
			++temp_count;
		}
		temp_count = 0;
		for (auto iter = hereWbcNodes->begin(); iter != hereWbcNodes->end(); iter++)
		{
			wbc[temp_count] = hereMapping.global2local(*iter);

			++temp_count;
		}
		cur_bi.nwbc = hereWbcNodes->size();
		cur_bi.nobc = hereObcNodes->size();
		delete[] hereWbcNodes;
		delete[] hereObcNodes;
		std::cout << "meshtag " << m + 1 << " proc " << cur_proc << " nwbc and nobc are " << cur_bi.nwbc << " " << cur_bi.nobc << '\n';
		cur_bi.wbc = wbc;
		cur_bi.obc = obc;
		recorder.close();
		delete[] temp_xyz;
	}
}

void CobaltReader::readOversetBoundary()
{
	if (cur_proc == 0)
	{
		std::cout << "Entering readOversetBoundary() " << cur_proc << "\n";
	}

	for (int m = 0; m < bis.size(); m++)
	{
		auto &cur_bi = bis[m];
		std::string filename = d_filenames[m];

		std::ifstream finGrid;

		// First read cellfile to fill the nodes

		// We need to first read how many cells there are.

		finGrid.open(filename, std::ios::in);
		if (!finGrid.is_open())
		{
			throw std::runtime_error("Open grid file failure\n");
		}
		int temp_nnodes, temp_nedges, temp_ncells;
		int to_throw_away;
		finGrid >> to_throw_away >> to_throw_away >> to_throw_away;
		finGrid >> temp_nnodes >> temp_nedges >> temp_ncells >> to_throw_away >> to_throw_away;
		// Now we have got vconn and nc and nv and ntypes all set, it's time to get node informations.
		double *temp_xyz = new double[temp_nnodes * d_dim];
		int temp_count = 0;
		for (temp_count = 0; temp_count < temp_nnodes; temp_count++)
		{
			for (int j = 0; j < d_dim; j++)
			{
				finGrid >> temp_xyz[d_dim * temp_count + j];
			}
		}
		int edgeNodes[10];
		int edgeCells[2];
		int temp_int = 0; // Used to record ptr

		int wall_temp_int = 0;

		for (int i = 0; i < temp_nedges; i++)
		{

			finGrid >> temp_count;

			for (int j = 0; j < temp_count; j++)
			{
				finGrid >> edgeNodes[j];
			}
			finGrid >> edgeCells[0] >> edgeCells[1];

			if (edgeCells[1] == signWall)
			{
				meshes_wallPtr[m].push_back(wall_temp_int);
				wall_temp_int += temp_count;
				for (int j = 0; j < temp_count; j++)
				{
					meshes_wallInd[m].push_back(edgeNodes[j]);
				}
			}
			else if (edgeCells[1] == signOver)
			{

				meshes_overPtr[m].push_back(temp_int); // Starting from 0
				temp_int += temp_count;
				for (int j = 0; j < temp_count; j++)
				{
					meshes_overInd[m].push_back(edgeNodes[j]); // Starting from 1
				}
			}
			else if (edgeCells[1] > 0)
			{
			}
			else
			{
				std::cout << "Pre is " << edgeNodes[0] << " " << edgeNodes[1] << '\n';
				std::cout << "proc and meshtag is " << cur_proc << " " << m + 1 << " " << i << " out of " << temp_nedges << " : " << edgeCells[0] << " " << edgeCells[1] << '\n';
				std::cout << "Something wrong,stop\n";
			}
		}
		finGrid.close();
		meshes_overPtr[m].push_back(temp_int);
		meshes_wallPtr[m].push_back(wall_temp_int);
		meshes_xyz[m] = temp_xyz;
	}
}

void CobaltReader::writeFiles()
{
	int mtag = 1;
	auto itera = d_filenames.begin();

	for (auto iter = bis.begin(); iter != bis.end(); iter++)
	{
		iter->setDim(d_dim);
		iter->make_iblank();
		iter->write_grd(*itera, mtag, cur_proc);
		iter->write_scatter_wbc(*itera + "wbc", mtag, cur_proc);
		iter->write_scatter_obc(*itera + "obc", mtag, cur_proc);
		iter->write_scatter_gwbc(*itera + "gwbc", mtag, cur_proc);
		iter->write_scatter_gobc(*itera + "gobc", mtag, cur_proc);

		itera++;
		mtag++;
	}

	for (int i = 0; i < d_nmesh; i++)
	{
		UBs2[i].write_grd(d_filenames[i] + "unstruct", i + 1, cur_proc);
		UBs2[i].write_grdSendRecv(d_filenames[i] + "UnstructSRFlags", i + 1, cur_proc);
	}
	for (int i = 0; i < d_filenames.size(); i++)
	{
		std::cout << i << " " << d_filenames[i] << " " << meshes_overInd[i].size() << " " << meshes_overPtr[i].size() - 1 << '\n';
	}
}
void CobaltReader::performCommunication()
{
}

void CobaltReader::readForTioga()
{
	selectReadFile();
	MPI_Barrier(MPI_COMM_WORLD);
	std::cout << "Exiting readForTioga\n";
	for (auto iter = bis.begin(); iter != bis.end(); iter++)
	{
		iter->make_iblank();
	}
}
void CobaltReader::readForSamrai()
{
	readOversetBoundary();
}

void CobaltReader::readForUnstruct()
{
	selectReadForUnstruct();
}