//#include<iostream>
#include<vector>
#include<string>
#include"codetypes.h"
#include"ADT.h"

class CartGrid;
class parallelComm;

class MeshBlock
{
public:
//private:
	int d_dim = 2;

	int nnodes = 0;
	int ncells = 0;
	int ntypes = 0;
	int* nv = NULL;
	int* nc = NULL;


	int nobc = 0;
	int nwbc = 0;

	double* x = NULL;
	
	int** vconn = NULL;
	int* wbcnode = NULL;
	int* obcnode = NULL;

	double* nodeRes = NULL;
	double* nodeResWeight = NULL;

	double* cellRes = NULL;
	double* cellResWeight = NULL;


	double* elementBbox = NULL;
	int* elementList = NULL;
	//int* elementList = NULL;

	ADT* adt = NULL;

	DONORLIST** donorList = NULL;
	
	int ninterp = 0;
	int interpListSize = 0;
	INTERPLIST* interpList = NULL;

	int* interp2donor = NULL;
	int* cancelArray = NULL;
	INTEGERLIST* cancelList = NULL;
	int ncancel = 0;

	void (*get_nodes_per_cell)(int*, int*);
	void (*get_receptor_nodes)(int*, int*, double*);
	void (*donor_inclusion_test)(int*, double*, int*, double*);
	void (*donor_frac)(int*, double*, int*, int*, double*, double*, int*);
	void (*convert_to_modal)(int*, int*, double*, int*, int*, double*);

	int nreceptorCells = 0;
	int* ctag = NULL;
	int* pointsPerCell = NULL;
	int maxPointsPerCell = 0;
	double* rxyz = NULL;

	int ipoint = 0;
	int* picked = NULL;

	int nreceptorCellsCart = 0;
	int* ctag_cart = NULL;
	int* pickedCart = NULL;
	int uniform_hex = 0;
	double dx[3];
	double xlow[3];
	int idims[3];
	int* uindx = NULL;
	int* invmap = NULL;     // inverse map
	int* mapmask = NULL;    // mask
	int* icft = NULL; 	   // frequency table for nodal containment
	int mapdims[3];  // dimensions of the map
	double mapdx[3]; // sides of the map
public:
	int* iblank = NULL;      /** < iblank value for each grid node */
	int* iblank_reduced = NULL;
	int* iblank_cell = NULL;
	int ntotalPointsCart;
	double* rxyzCart = NULL;
	int* donorIdCart = NULL;
	int donorListLength = 0;

	int nfringe = 1;
	int mexclude = 3;
	int meshtag = -1; /** < tag of the mesh that this block belongs to */
	int check_uniform_hex_flag = 0;
	double resolutionScale = 1.0;
	//
	// oriented bounding box of this partition
	// 
	OBB* obb = NULL;
	OBB* obh = NULL;
	//
	int nsearch = 0;        /** < number of query points to search in this block */
	int* isearch = NULL;       /** < index of query points in the remote process */
	int* tagsearch = NULL;       /** < index of query points in the remote process */
	double* res_search = NULL;   /** < resolution of search points */
	int* xtag = NULL;            /** < hash to determine if there are duplicates */
	double* xsearch = NULL;    /** < coordinates of the query points */
	double* rst = NULL;            /**  natrural coordinates */
	int* donorId = NULL;       /** < donor indices for those found */
	
	
	int donorCount = 0;
	int myid = -1;
	//double* cellRes = NULL;  /** < resolution for each cell */
	int ntotalPoints = 0;        /**  total number of extra points to interpolate */
	int ihigh = 0;
	int ninterp2 = 0;            /** < number of interpolants for high-order points */
	int interp2ListSize = 0;
	INTERPLIST* interpList2 = NULL; /** < list for high-interpolation points */
	int ninterpCart = 0;
	int interpListCartSize = 0;
	INTERPLIST* interpListCart = NULL;
	int* receptorIdCart = NULL;
	//User specified
	uint64_t* cellGID = NULL;
	uint64_t* nodeGID = NULL;

	double* userSpecifiedNodeRes = NULL;
	double* userSpecifiedCellRes = NULL;
	std::vector<uint64_t> gid_search; /**< Global node ID for the query points */

	//
  // call back functions to use p4est to search
  // its own internal data
  //
	void (*p4estsearchpt) (double*, int*, int*, int*);
	void (*check_intersect_p4est) (int*, int*);

	MeshBlock() {};

	~MeshBlock();
	void setData(int dim, int btag, int nnodesi, double* xyzi, int* ibli, int nwbci, int nobci,
		int* wbcnodei, int* obcnodei,
		int ntypesi, int* nvi, int* nci, int** vconni,
		uint64_t* cell_gid = NULL, uint64_t* node_gid = NULL);

	void preprocess();

	void tagBoundary();

	void writeGridFile(const std::string& filename);
	void writeCellFile(const std::string& filename);
	void writeOBB(const std::string& filename);


	void writeFlowFile(int bid, double* q, int nvar, int type);

	


	//void setResolutions(double* nres, double* cres);

	void search();
	//void search_uniform_hex();
	
	void writeOBB2(const std::string& filename, OBB* obc);
	//void writeOBB2(OBB* obc, int bid);

	void updateSolnData(int inode, double* qvar, double* q, int nvar, int interptype);

	int getNinterp(void) { return ninterp; };

	void getInterpolatedSolution(int* nints, int* nreals, int** intData, double** realData, double* q,
		int nvar, int interptype);

	void getInterpolatedSolutionAMR(int* nints, int* nreals, int** intData, double** realData, double* q,
		int nvar, int interptype);

	void checkContainment(int* cellIndex, int adtElement, double* xsearch);
	void checkContainment(int* cellIndex, int adtElement, double* xsearch,bool indicator);
	void getWallBounds(int* mtag, int* existWall, double wbox[6]);

	void markWallBoundary(int* sam, int nx[3], double extents[6]);

	void getQueryPoints2(OBB* obb, int* nints, int** intData, int* nreals,
		double** realData);
	void getQueryPoints(OBB* obb, int* nints, int** intData, int* nreals,
		double** realData);


	/** routines that do book keeping */

	void getDonorPacket(PACKET* sndPack, int nsend);

	void initializeDonorList();

	void insertAndSort(int pointid, int senderid, int meshtag, int remoteid, double donorRes, double receptorRes);

	void processDonors(HOLEMAP* holemap, int nmesh, int** donorRecords, double** receptorResolution,
		int* nrecords);

	void initializeInterpList(int ninterp_input);

	void findInterpData(int* recid, int irecord, double receptorRes);


	void findInterpListCart();

	void set_ninterp(int);

	void getCancellationData(int* nints, int** intData);

	void getCancellationReduce(int* nints, int** intData);

	void cancelDonor(int irecord);

	void getInterpData(int* nrecords, int** intData);

	void clearIblanks(void);

	void getStats(int mstat[2]);

	void setIblanks(int inode);

	void getDonorCount(int* dcount, int* fcount);

	void getDonorInfo(int* receptors, int* indices, double* frac);

	void getReceptorInfo(int* receptors);

	void getReducedOBB(OBB*, double*);
	void getReducedOBB2(OBB*, double*);

	void resetCoincident();
	//
	// routines for high order connectivity and interpolation
	//
	void getCellIblanks(void);
	void getCellIblanks2(void);
	void set_cell_iblank(int* iblank_cell_input)
	{
		iblank_cell = iblank_cell_input;
	}
	void set_node_iblank(int* iblank_node_input)
	{
		iblank = iblank_node_input;
	}
	void setcallback(void (*f1)(int*, int*),
		void (*f2)(int*, int*, double*),
		void (*f3)(int*, double*, int*, double*),
		void (*f4)(int*, double*, int*, int*, double*, double*, int*),
		void (*f5)(int*, int*, double*, int*, int*, double*))
	{
		get_nodes_per_cell = f1;
		get_receptor_nodes = f2;
		donor_inclusion_test = f3;
		donor_frac = f4;
		convert_to_modal = f5;
	}

	void setp4estcallback(void (*f1)(double*, int*, int*, int*),
		void (*f2)(int*, int*))
	{
		p4estsearchpt = f1;
		check_intersect_p4est = f2;
	}

	
	void getInternalNodes(void);
	void getExtraQueryPoints(OBB* obb, int* nints, int** intData, int* nreals,
		double** realData);
	void processPointDonors(void);
	void getInterpolatedSolutionAtPoints(int* nints, int* nreals, int** intData,
		double** realData,
		double* q,
		int nvar, int interptype);
	void updatePointData(double* q, double* qtmp, int nvar, int interptype);
	void outputOrphan(FILE* fp, int i)
	{
		fprintf(fp, "%f %f %f\n", rxyz[3 * i], rxyz[3 * i + 1], rxyz[3 * i + 2]);
	}
	void clearOrphans(HOLEMAP* holemap, int nmesh, int* itmp);
	void getUnresolvedMandatoryReceptors();
	void getCartReceptors(CartGrid* cg, parallelComm* pc);
	void setCartIblanks();

	// Getters
	inline int getMeshTag() const { return meshtag + (1 - BASE); }

	/**
	 * Get donor packet for multi-block/partition setups
	 *
	 */
	void getMBDonorPktSizes(std::vector<int>&,
		std::vector<int>&);

	void getMBDonorPackets(std::vector<int>&,
		std::vector<int>&,
		PACKET*);

	/** Reset interpolation list data structure
	 *
	 *  Reset the data structures in situations where the performConnectivity
	 *  method is invoked at every timestep when meshes undergo relative motion.
	 */
	void resetInterpData() {
		if (interpList) {
			for (int i = 0; i < interpListSize; i++) {
				if (interpList[i].inode)
					TIOGA_FREE(interpList[i].inode);
				if (interpList[i].weights)
					TIOGA_FREE(interpList[i].weights);
			}
			TIOGA_FREE(interpList);
		}
		ninterp = 0;
		interpListSize = 0;
	}
	void reduce_fringes();
private:
	void outPutSearch(const std::string& filename);
	void check_for_uniform_hex();

	void create_hex_cell_map();
	void checkContainment2(int* cellIndex, int adtElement, double* x2search);
	
public:
	

	void writeCellFile2(const std::string& filename, int* cellInfo);
	void writeCellFile2(const std::string& filename, double* cellInfo);
	void writeGridFile2(const std::string& filename, int* nodeInfo);
	void writeGridFile2(const std::string& filename, double* nodeInfo);

	void writeOBBPair(const std::string& filename,OBB* obc);
	void set_Iblank(int pid, int value)
	{
		iblank[pid] = value;
	}
	void writeMandatoryReceptor(const std::string& filename);
	void writeWbc(const std::string& filename);
	void writeObc(const std::string& filename);
	void get_distances();
	void search_distance(ADT* adt, double* distance, int WorO);
	void set_global_bc(int inwbc, int inobc, double* iwbcXYZ, double* iobcXYZ);
	int ngWbc = 0;
	int ngObc = 0;
	double* gWbcXYZ = NULL;
	double* gObcXYZ = NULL;
	double* wallDistance = NULL;
	double* overDistance = NULL;
	void reResolution();

	void write_distance_files(const std::string& filename)
	{
		if (wallDistance)
		{
			writeGridFile2(filename+"Wall",wallDistance);
		}
		if (overDistance)
		{
			writeGridFile2(filename+"Over",overDistance);
		}
		
	}
	void write_rxyz(const std::string& filename);
	//int* gWbcNode;//Do i need the global index when no global index is maintained here?
	//int* gObcNode;
	void writeGridHDF(const std::string& filename);
	void writeGridFile2(const std::string& filename, double* nodeInfo, int nvar);
	void writeCellFile2(const std::string& filename,int nvar, double* cellInfo);
};