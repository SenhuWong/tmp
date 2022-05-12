#include "Euler2D.h"
#include "HLLCFluxStrategy.h"

double computeEdgeNormal(int dim, double xv[4][3], int nvert, double xnorm[3], bool counter_clock = false);
#include <math.h>
#include "mpi.h"
#include <hdf5.h>
#include <algorithm>
class PressureHolder2D
{
    public:

    GeomElements::vector3d<2,double> center;
    double normal_y;
    double cp;
    PressureHolder2D(double x,double y,double normal,double cp_in)
    {
        center[0] = x;
        center[1] = y;
        normal_y = normal;
        cp = cp_in;
    }
    PressureHolder2D(PressureHolder2D& another)
    {
        center[0] = (another.center)[0];
        center[1] = (another.center)[1];
        normal_y = another.normal_y;
        cp = another.cp;
    }
    PressureHolder2D(const PressureHolder2D& another)
    {
        center = another.center;
        normal_y = another.normal_y;
        cp = another.cp;
    }

};
//Increasing ranking: e1<e2
bool rankWithXY(PressureHolder2D &e1,PressureHolder2D &e2)
{
    //If e1's normaly is negative(downward, upper face) and e2's normaly is positive (upward lower face)
    //return false means that upper face is e1 > e2
    if(e1.normal_y<0 and e2.normal_y>=0) return false;
    else if(e1.normal_y>=0 and e2.normal_y < 0) return true;
    else
    {
        if(e1.normal_y>0)//Both upper
        {
            if(e1.center[0] < e2.center[0])
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        else
        {
            if(e1.center[0] < e2.center[0])
            {
                return false;
            }
            else
            {
                return true;
            }
        }
    }
}

Euler2D::Euler2D(UnstructTopologyHolder *hder)
    : TopologyHolderStrategy(hder, 4)
{
    // Communication set
    cur_proc = d_hder->curproc;
    num_proc = d_hder->numproc;
    // Topology holder set
    d_limiter = new Vankatakrishnan_Limiter<2>(d_hder, this);
#ifdef DEBUG
    d_limiter->setComm(cur_proc, num_proc);

#endif // DEBUG
    d_fluxComputer = new HLLCFluxStrategy(d_hder, this);
    d_fluxComputer->setComm(num_proc, cur_proc);

    d_boundaryHandler = new BoundaryStrategy(d_hder,this);
    // Allocate Data Storage for this problem
    registerTopology();
    // This should be done with an input_db
    TopologyHolderStrategy::set_fs_variable(0.75, 3, 1.225, 101325, 1);
    //
    TopologyHolderStrategy::update_nondim_variable();
}

void Euler2D::registerTopology()
{
    // std::cout<<"Topology registered\n";
    // Initialize Conservative Variables
    U = new double *[d_nmesh];
    // Update Edge Variable
    U_edge = new double *[d_nmesh];
    // With Edge Variables compute gradient
    gradU = new GeomElements::vector3d<2, double> *[d_nmesh];
    Spectrum_cell = new double *[d_nmesh];
    // Compute left and right value of Conservative Variables
    UL = new double *[d_nmesh];
    UR = new double *[d_nmesh];
    Spectrum_edge = new double *[d_nmesh];
    // Compute flux and fluxSum
    Flux_edge = new double *[d_nmesh];
    Residual = new double *[d_nmesh];

    // Compute Local dt
    // dt = new double **[d_nmesh];

    // Allocate space for sndBuffer
    // One for each communication requiring proc
    
    // normal_vector = new double *[d_nmesh];
    for (int i = 0; i < d_nmesh; i++)
    {
        U[i] = new double [d_NEQU*d_hder->nCells(i)];
        Residual[i] = new double [d_NEQU*d_hder->nCells(i)];
        gradU[i] = new GeomElements::vector3d<2, double> [d_NEQU*d_hder->nCells(i)];
        Spectrum_cell[i] = new double[d_hder->nCells(i)];

        U_edge[i] = new double [d_NEQU*d_hder->nEdges(i)];
        UL[i] = new double [d_NEQU*d_hder->nEdges(i)];
        UR[i] = new double [d_NEQU*d_hder->nEdges(i)];
        Spectrum_edge[i] = new double[d_hder->nEdges(i)];

        Flux_edge[i] = new double [d_NEQU*d_hder->nEdges(i)];
    }
}

// Initialize Data to t0
void Euler2D::initializeData()
{
    if (cur_proc == 0)
        std::cout << "Initializing Data\n";
    for (int i = 0; i < d_nmesh; i++)
    {
        for (int k = 0; k < d_hder->nCells(i); k++)
        {
            for (int j = 0; j < d_NEQU; j++)
            {
                U[i][j+d_NEQU*k] = fs_primVar[j];
            }
        }
    }
    for (int k = 0; k < d_NEQU; k++)
    {
        std::cout << fs_primVar[k] << '\n';
    }
}

void Euler2D::ReconstructionInitial()
{
    // Grad should be set to 0
    for (int i = 0; i < d_nmesh; i++)
    {

        for (int j = 0; j < d_NEQU; j++)
        {
            for (int k = 0; k < d_hder->nCells(i); k++)
            {
                gradU[i][j+d_NEQU*k].reset();
            }
        }
    }
    // Spectrum of cell should be set to 0
    for (int i = 0; i < d_nmesh; i++)
    {
        for (int k = 0; k < d_hder->nCells(i); k++)
        {
            Spectrum_cell[i][k] = 0;
            for (int j = 0; j < d_NEQU; j++)
            {
                Residual[i][j+d_NEQU*k] = 0;
            }
        }
    }
    // Limiter related should be set
    d_limiter->init();
}

// Boundary and normal edge's conservative update
void Euler2D::updateEdgeValues()
{
    d_boundaryHandler->HandleBoundary();
}
// Gradient computation for each cell and LeftRight Value Computation for each edge
// Since there ain't no cell to edge ptr, we have to make two loops

void Euler2D::solveGradient()
{
    for (int i = 0; i < d_nmesh; i++)
    {
        auto &curBlk = d_hder->blk2D[i];
        for (int k = 0; k < d_hder->nEdges(i); k++)
        {
            auto &curEdge = curBlk.d_localEdges[k];
            int lC = curEdge.lCInd();
            int rC = curEdge.rCInd();
            for (int j = 0; j < d_NEQU; j++)
            {
                gradU[i][j+d_NEQU*lC] = gradU[i][j+d_NEQU*lC] + curEdge.normal_vector() * (curEdge.area() * U_edge[i][j+d_NEQU*k]);
            }
            if (rC >= 0)
            {
                for (int j = 0; j < d_NEQU; j++)
                {
                    gradU[i][j+d_NEQU*rC] = gradU[i][j+d_NEQU*rC] - curEdge.normal_vector() * (curEdge.area() * U_edge[i][j+d_NEQU*k]);
                }
            }
        }
        for (int k = 0; k < d_hder->nCells(i); k++)
        {
            auto &curCell = curBlk.d_localCells[k];
            for (int j = 0; j < d_NEQU; j++)
            {
                gradU[i][j+d_NEQU*k] = gradU[i][j+d_NEQU*k] / curCell.volume();
            }
        }
    }
}

void Euler2D::solveSpectralRadius()
{
    GeomElements::vector3d<2, double> velocity_edge;
    for (int i = 0; i < d_nmesh; i++)
    {
        auto &curBlk = d_hder->blk2D[i];
        for (int k = 0; k < d_hder->nEdges(i); k++)
        {
            auto &curEdge = curBlk.d_localEdges[k];
            velocity_edge[0] = U_edge[i][2+d_NEQU*k];
            velocity_edge[1] = U_edge[i][3+d_NEQU*k];
            double c = sqrt(Gamma * U_edge[i][1+d_NEQU*k] / U_edge[i][0+d_NEQU*k]);
            double Vn = velocity_edge.dot_product(curEdge.normal_vector());
            int lC = curEdge.lCInd();
            int rC = curEdge.rCInd();
            Spectrum_edge[i][k] = (std::abs<double>(Vn) + c) * curEdge.area();
            Spectrum_cell[i][lC] += Spectrum_edge[i][k];
            if (rC >= 0)
            {
                Spectrum_cell[i][rC] += Spectrum_edge[i][k];
            }
        }
    }
}

void Euler2D::updateEdgeLeftRightValues()
{
    const double EPSILON = 1.0e-5;
    // Compute Limiter for each cell
    // This should be done by a Limiter class
    d_limiter->computeLimiter();
    bool hasTrouble = false;
    // Compute left and right value for each Edge
    for (int i = 0; i < d_nmesh; i++)
    {
        auto &curBlk = d_hder->blk2D[i];
        for (int k = 0; k < curBlk.nEdges(); k++)
        {

            auto &curEdge = curBlk.d_localEdges[k];
            int lc = curEdge.lCInd();
            int rc = curEdge.rCInd();

            if (lc >= 0 and rc >= 0)
            {
                for (int j = 0; j < d_NEQU; j++)
                {
                    UL[i][j+d_NEQU*k] = UL[i][j+d_NEQU*k] * d_limiter->getLimiter(i, j, lc) + U[i][j+d_NEQU*lc];
                }
                for (int j = 0; j < d_NEQU; j++)
                {
                    UR[i][j+d_NEQU*k] = UR[i][j+d_NEQU*k] * d_limiter->getLimiter(i, j, rc) + U[i][j+d_NEQU*rc];
                }
            }
            else
            {
                if (lc < 0)
                {
                    throw std::runtime_error("Left Cell can't be NULL\n");
                }
                if (rc == GeomElements::edge3d<2>::BoundaryType::WALL) // Wall, set to be the cell variables
                {
                    for (int j = 0; j < d_NEQU; j++)
                    {
                        UL[i][j+d_NEQU*k] = UR[i][j+d_NEQU*k] = U[i][j+d_NEQU*lc];
                    }
                }
                else if (rc == GeomElements::edge3d<2>::BoundaryType::FARFIELD) // Far_field, set to be the face variables
                {
                    for (int j = 0; j < d_NEQU; j++)
                    {
                        UL[i][j+d_NEQU*k] = UR[i][j+d_NEQU*k] = U_edge[i][j+d_NEQU*k];
                    }
                }
                else
                {
                    std::cout << lc << "," << rc << '\n';
                    std::cout << "SOmething wrong with left and right cell\n";
                    std::cin.get();
                }
            }
        }
    }
    if (hasTrouble)
    {
        throw std::runtime_error("Limiter Error");
    }
}
// Update flux with the corresponding left and right value
void Euler2D::updateFlux()
{
    d_fluxComputer->computeFlux();
    // First clear FLuxsum from last iteration
    for (int i = 0; i < d_nmesh; i++)
    {
        for (int j = 0; j < d_NEQU; j++)
        {
            for (int k = 0; k < d_hder->nCells(i); k++)
            {
                Residual[i][j+d_NEQU*k] = 0;
            }
        }
    }
    for (int i = 0; i < d_nmesh; i++)
    {
        auto &curBlk = d_hder->blk2D[i];

        for (int k = 0; k < d_hder->nEdges(i); k++)
        {

            auto &curEdge = curBlk.d_localEdges[k];
            int lc = curEdge.lCInd();
            int rc = curEdge.rCInd();

            if (lc >= 0) // Left plus right minus
            {
                for (int j = 0; j < d_NEQU; j++)
                {
                    Residual[i][j+d_NEQU*lc] += Flux_edge[i][j+d_NEQU*k];
                }
            }
            else
            {
                std::cout << "left cell ind can't be negative\n";
                std::cin.get();
            }
            if (rc >= 0)
            {
                for (int j = 0; j < d_NEQU; j++)
                {
                    Residual[i][j+d_NEQU*rc] -= Flux_edge[i][j+d_NEQU*k];
                }
            }
        }
    }
}


void Euler2D::outPutCp(std::string& filename, int mesh_ind)
{
    auto& curBlk = d_hder->blk2D[mesh_ind];
    std::set<int,std::less<int>>* recvCellSet = curBlk.getRecvIndSet();
    std::vector<double> raw_data;
    for(int k = 0;k<d_hder->nEdges(mesh_ind);k++)
    {
        auto& curEdge = curBlk.d_localEdges[k];
        int lC = curEdge.lCInd();
        int rC = curEdge.rCInd();
        if(rC==GeomElements::edge3d<2>::BoundaryType::WALL)
        {
            //Check if lC is inside the recvCell
            if(recvCellSet->find(lC)==recvCellSet->end())
            {
                double cp = (U[mesh_ind][1+d_NEQU*lC] - fsnd_pressure)*2/(fsnd_density*(fsnd_velocity_components[0]*fsnd_velocity_components[0]+fsnd_velocity_components[1]*fsnd_velocity_components[1]));
                raw_data.push_back(curEdge.center()[0]);
                raw_data.push_back(curEdge.center()[1]);
                raw_data.push_back(curEdge.normal_vector()[1]);
                raw_data.push_back(cp);
            }
        }
    }
    delete[] recvCellSet;
    int nlocalWall = raw_data.size()/4;
    int* nEachLocalWall = new int[num_proc+1];

    MPI_Allgather(&nlocalWall,1,MPI_INT,nEachLocalWall+1,1,MPI_INT,MPI_COMM_WORLD);

    nEachLocalWall[0] = 0;
    for(int i = 0;i<num_proc;i++)
    {
        nEachLocalWall[i+1] = nEachLocalWall[i] + nEachLocalWall[i+1];
    }
    
    //I'm using HDF5 Parallel IO feature for this shit.
    
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    



    hsize_t dims1[2];

    hid_t acc_tpl = H5Pcreate(H5P_FILE_ACCESS);

    herr_t status = H5Pset_fapl_mpio(acc_tpl,comm,info);
    
    hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC,H5P_DEFAULT,acc_tpl);

    status = H5Pclose(acc_tpl);
    
    dims1[0] = nEachLocalWall[num_proc];
    dims1[1] = 4;

    hid_t dataspace_id = H5Screate_simple(2,dims1,NULL);

    hid_t dataset_id = H5Dcreate(file_id,"Cp",H5T_NATIVE_DOUBLE,dataspace_id,
                                H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

    //Now each Proc has its own slice to fill in.
    
    if(nEachLocalWall[cur_proc+1]-nEachLocalWall[cur_proc]!=0)
    {
        hsize_t count[2] = {nlocalWall,4};
        hsize_t offset[2] = {nEachLocalWall[cur_proc],0};
        hsize_t stride[2] = {1,1};
        hsize_t block[2] = {1,1};
        //Need to write in

        hid_t memspace_id = H5Screate_simple(2,count,NULL);

        status = H5Sselect_hyperslab(dataspace_id,H5S_SELECT_SET, offset, stride, count, block);

        status = H5Dwrite(dataset_id,H5T_NATIVE_DOUBLE,memspace_id,dataspace_id,H5P_DEFAULT,raw_data.data());

        status = H5Sclose(memspace_id);
    }
    status = H5Sclose(dataspace_id);
    status = H5Dclose(dataset_id);
    status = H5Fclose(file_id);

    MPI_Barrier(MPI_COMM_WORLD);
    if(cur_proc==0)
    {
        file_id = H5Fopen(filename.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
        dataset_id = H5Dopen(file_id,"Cp",H5P_DEFAULT);

        dataspace_id =H5Dget_space(dataset_id);
        double* buffer = new double[nEachLocalWall[num_proc]*4];
        status = H5Dread(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,dataspace_id,H5P_DEFAULT,buffer);
        std::vector<PressureHolder2D> holderVector;
        for(int i = 0;i<nEachLocalWall[num_proc];i++)
        {
            holderVector.push_back(PressureHolder2D(buffer[4*i],buffer[4*i+1],buffer[4*i+2],buffer[4*i+3]));
        }
        std::sort(holderVector.begin(),holderVector.end(),rankWithXY);
        for(int i = 0;i<nEachLocalWall[num_proc];i++)
        {
            buffer[4*i] = holderVector[i].center[0];
            buffer[4*i+1] = holderVector[i].center[1];
            buffer[4*i+2] = holderVector[i].normal_y;
            buffer[4*i+3] = holderVector[i].cp;
        }
        status = H5Dwrite(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,dataspace_id,H5P_DEFAULT,buffer);

        status = H5Sclose(dataspace_id);
        status = H5Dclose(dataset_id);
        status = H5Fclose(file_id);
        delete[] buffer;
    }

    delete[] nEachLocalWall;

}