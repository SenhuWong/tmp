#include "ViscousFlow2D.h"
#include "HLLCFluxStrategy.h"

double computeEdgeNormal(int dim, double xv[4][3], int nvert, double xnorm[3], bool counter_clock = false);
#include <math.h>
#include "mpi.h"
#include <hdf5.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
class PressureHolder2D;
//Increasing ranking: e1<e2
bool rankWithXY(PressureHolder2D &e1,PressureHolder2D &e2);
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
void ViscousFlow2D::initializeSubStrategy()
{
    // Communication set
    cur_proc = d_hder->curproc;
    num_proc = d_hder->numproc;
    // Topology holder set
    auto limiter = new Vankatakrishnan_Limiter<2>(d_hder, this);
#ifdef DEBUG
    limiter->setComm(cur_proc, num_proc);

#endif // DEBUG
    d_cfluxComputer = new HLLCFluxStrategy(d_hder, this, limiter);
    d_cfluxComputer->setComm(num_proc, cur_proc);
    d_vfluxComputer = new LaminarFluxStrategy(d_hder,this);
    // Allocate Data Storage for this problem
    d_boundaryHandler = new BoundaryStrategy(d_hder,this);
    registerTopology();

    d_turbulentModel = new SSTkomegaModel<2>(d_hder,this,d_vfluxComputer);
    turbModelTimeStrategy = new LUSGSIntegrator(d_hder,d_turbulentModel);

}
ViscousFlow2D::ViscousFlow2D(UnstructTopologyHolder *hder)
    : TopologyHolderStrategy(hder,4)
{
    // Communication set
    cur_proc = d_hder->curproc;
    num_proc = d_hder->numproc;
    // Topology holder set
    auto limiter = new Vankatakrishnan_Limiter<2>(d_hder, this);
#ifdef DEBUG
    limiter->setComm(cur_proc, num_proc);

#endif // DEBUG
    TopologyHolderStrategy::set_fs_variable3(0.73, 2.79, 288.0, 101325, 6.5e6);
    d_cfluxComputer = new HLLCFluxStrategy(d_hder, this, limiter);
    d_cfluxComputer->setComm(num_proc, cur_proc);
    d_vfluxComputer = new LaminarFluxStrategy(d_hder,this);
    // Allocate Data Storage for this problem
    d_boundaryHandler = new BoundaryStrategy(d_hder,this);
    registerTopology();

    d_turbulentModel = new SSTkomegaModel<2>(d_hder,this,d_vfluxComputer);
    turbModelTimeStrategy = new LUSGSIntegrator(d_hder,d_turbulentModel);
    
}

void ViscousFlow2D::registerTopology()
{
    // std::cout<<"Topology registered\n";
    // Initialize Conservative Variables
    U = new double *[d_nmesh];
    // Update Edge Variable
    U_edge = new double *[d_nmesh];
    // With Edge Variables compute gradient
    gradU = new GeomElements::vector3d<2, double> *[d_nmesh];
    gradT = new GeomElements::vector3d<2, double> *[d_nmesh];
    Spectrum_cell_c = new double *[d_nmesh];
    Spectrum_cell_v = new double *[d_nmesh];
    // Compute left and right value of Conservative Variables
    Residual = new double *[d_nmesh];
    d_stableDt = new double*[d_nmesh];
    // normal_vector = new double *[d_nmesh];
    for (int i = 0; i < d_nmesh; i++)
    {
        d_stableDt[i] = new double[d_hder->nCells(i)];
        U[i] = new double [d_NEQU*d_hder->nCells(i)];
        Residual[i] = new double [d_NEQU*d_hder->nCells(i)];
        gradU[i] = new GeomElements::vector3d<2, double> [d_NEQU*d_hder->nCells(i)];
        gradT[i] = new GeomElements::vector3d<2, double> [d_hder->nCells(i)];
        Spectrum_cell_c[i] = new double[d_hder->nCells(i)];
        Spectrum_cell_v[i] = new double[d_hder->nCells(i)];

        U_edge[i] = new double [d_NEQU*d_hder->nEdges(i)];
    }
}

// Initialize Data to t0
void ViscousFlow2D::initializeData()
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
        std::cout<<"fs_primVars:" << fs_primVar[k] << '\n';
    }
    turbModelTimeStrategy->initializeData();
}

void ViscousFlow2D::ReconstructionInitial()
{
    // Grad should be set to 0
    for (int i = 0; i < d_nmesh; i++)
    {
        for (int k = 0; k < d_hder->nCells(i); k++)
        {
            gradT[i][k].reset();
            for (int j = 0; j < d_NEQU; j++)
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
            Spectrum_cell_c[i][k] = 0;
            Spectrum_cell_v[i][k] = 0;
            for (int j = 0; j < d_NEQU; j++)
            {
                Residual[i][j+d_NEQU*k] = 0;
            }
        }
    }
    // Limiter is left for FLuxStrategy to handle.
}

// Boundary and normal edge's conservative update
void ViscousFlow2D::updateEdgeValues()
{
    d_boundaryHandler->HandleBoundary();
}
// Gradient computation for each cell and LeftRight Value Computation for each edge
// Since there ain't no cell to edge ptr, we have to make two loops

void ViscousFlow2D::solveGradient(int istage)
{
    for (int i = 0; i < d_nmesh; i++)
    {
        auto &curBlk = d_hder->blk2D[i];
        for (int k = 0; k < d_hder->nEdges(i); k++)
        {
            auto &curEdge = curBlk.d_localEdges[k];
            int lC = curEdge.lCInd();
            int rC = curEdge.rCInd();
            //Compute T at left and right
            double edgeT = U_edge[i][1+d_NEQU*k]/U_edge[i][0+d_NEQU*k];
            
            gradT[i][lC] = gradT[i][lC] + curEdge.normal_vector() * (curEdge.area() * edgeT);
            for (int j = 0; j < d_NEQU; j++)
            {
                gradU[i][j+d_NEQU*lC] = gradU[i][j+d_NEQU*lC] + curEdge.normal_vector() * (curEdge.area() * U_edge[i][j+d_NEQU*k]);
            }
            if (rC >= 0)
            {
                gradT[i][rC] = gradT[i][rC] - curEdge.normal_vector() * (curEdge.area() * edgeT);
                for (int j = 0; j < d_NEQU; j++)
                {
                    gradU[i][j+d_NEQU*rC] = gradU[i][j+d_NEQU*rC] - curEdge.normal_vector() * (curEdge.area() * U_edge[i][j+d_NEQU*k]);
                }
            }
        }
        for (int k = 0; k < d_hder->nCells(i); k++)
        {
            auto &curCell = curBlk.d_localCells[k];
            gradT[i][k] = gradT[i][k] / curCell.volume();
            for (int j = 0; j < d_NEQU; j++)
            {
                gradU[i][j+d_NEQU*k] = gradU[i][j+d_NEQU*k] / curCell.volume();
            }
        }
    }
}

void ViscousFlow2D::solveSpectralRadius()
{
    double Spectrum_edge;
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
            Spectrum_edge = (fabsf64(Vn) + c) * curEdge.area();
            Spectrum_cell_c[i][lC] += Spectrum_edge;
            if (rC >= 0)
            {
                Spectrum_cell_c[i][rC] += Spectrum_edge;
            }
            auto& leftCell = curBlk.d_localCells[lC];
            Spectrum_edge = std::max<double>(4/(3.0*U_edge[i][0+d_NEQU*k]),Gamma/U_edge[i][0+d_NEQU*k])*getMuOverPrEdge(i,k)*curEdge.area()*curEdge.area();
            Spectrum_cell_v[i][lC] += Spectrum_edge / leftCell.volume();
            if(rC >= 0)
            {
                auto& rightCell = curBlk.d_localCells[rC];
                Spectrum_cell_v[i][rC] += Spectrum_edge /  rightCell.volume();
            }
        }
    }
}


void ViscousFlow2D::updateFlux(int istage)
{
    d_cfluxComputer->computeFlux();

    d_vfluxComputer->computeFlux(istage);
}

void ViscousFlow2D::outPutCp(std::string& filename, int mesh_ind, int step)
{
    auto& curBlk = d_hder->blk2D[mesh_ind];
    // int
    double* p_average = new double[d_hder->nPoints(mesh_ind)];
    int* count_average = new int[d_hder->nPoints(mesh_ind)];
    for(int i = 0;i<d_hder->nPoints(mesh_ind);i++)
    {
        count_average[i] = 0;
        p_average[i] = 0;
    }
    for(int i = 0;i<d_hder->nCells(mesh_ind);i++)
    {
        auto& curCell = curBlk.d_localCells[i];

        for(int k = 0;k<curCell.size();k++)
        {
            int nodeInd = curCell.pointInd(k);
            
            ++count_average[nodeInd];
            p_average[nodeInd] += U[mesh_ind][1+d_NEQU*i];

        }
    }
    double fs_rhoV2 = fsnd_density*(fsnd_velocity_components[0]*fsnd_velocity_components[0]+fsnd_velocity_components[1]*fsnd_velocity_components[1]);
    
    for(int i = 0;i<d_hder->nPoints(mesh_ind);i++)
    {
        p_average[i] = p_average[i]/count_average[i];
        // std::cout<<count_average[i]<<'\n';
        p_average[i] = (p_average[i] - fsnd_pressure)*2/fs_rhoV2;
    }

    std::set<int,std::less<int>>* wallNodeInds = new std::set<int,std::less<int>>[1];
    

    for(int k = 0;k<d_hder->nEdges(mesh_ind);k++)
    {
        auto& curEdge = curBlk.d_localEdges[k];
        int rC = curEdge.rCInd();
        if(rC==GeomElements::edge3d<2>::BoundaryType::WALL)
        {
            for(int j = 0;j<curEdge.size();j++)
            {
                
                int ptId = curEdge.pointInd(j);
                if(curEdge.normal_vector()[1]>0)//Upper half
                {
                    //std::cout<<"Proc "<<cur_proc<<" found one lower half"<<'\n';
                    count_average[ptId] = -std::abs(count_average[ptId]);
                }
                wallNodeInds->insert(ptId);
            }
        }
    }
    int nlocalWall = wallNodeInds->size();
    //Xloc,Yloc,Count,
    double* realData;
    if(nlocalWall>0)
    {
        realData = new double[nlocalWall*4];
    }
    int m = 0;
    for(auto iter = wallNodeInds->begin();iter!=wallNodeInds->end();iter++)
    {
        int nodeId = *iter;
        auto curPoint = curBlk.d_localPoints[nodeId].d_pos;
        realData[m++] = curPoint[0];
        realData[m++] = curPoint[1];
        realData[m++] = count_average[nodeId];
        realData[m++] = p_average[nodeId];
    }

    // std::ofstream fout0;
    // fout0.open("Cp0"+std::to_string(cur_proc));
    // for(int i = 0;i<m/4;i++)
    // {
    //     fout0<<realData[4*i]<<"\t"<<realData[4*i+1]<<"\t"<<realData[4*i+2]<<"\t"<<realData[4*i+3]<<'\n';
    // }
    // fout0.close();

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

    hid_t dataset_id = H5Dcreate(file_id,("Cp"+std::to_string(step)).c_str(),H5T_NATIVE_DOUBLE,dataspace_id,
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

        status = H5Dwrite(dataset_id,H5T_NATIVE_DOUBLE,memspace_id,dataspace_id,H5P_DEFAULT,realData);

        status = H5Sclose(memspace_id);
    }
    status = H5Sclose(dataspace_id);
    status = H5Dclose(dataset_id);
    status = H5Fclose(file_id);
    MPI_Barrier(MPI_COMM_WORLD);
    if(cur_proc==0)
    {
        file_id = H5Fopen(filename.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
        dataset_id = H5Dopen(file_id,("Cp"+std::to_string(step)).c_str(),H5P_DEFAULT);

        dataspace_id =H5Dget_space(dataset_id);
        double* buffer = new double[nEachLocalWall[num_proc]*4];
        status = H5Dread(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,dataspace_id,H5P_DEFAULT,buffer);
        std::vector<PressureHolder2D> holderVector;
        for(int i = 0;i<nEachLocalWall[num_proc];i++)
        {
            holderVector.push_back(PressureHolder2D(buffer[4*i],buffer[4*i+1],buffer[4*i+2],buffer[4*i+3]));   
        }
        std::sort(holderVector.begin(),holderVector.end(),rankWithXY);
        //Remove the redundant ones.
        bool TheSame = true;
        int cur,next;
        int theSameCount = 0;
        m = 0;
        for(int i = 0;i<holderVector.size();i++)
        {
            cur = i;
            //Add cur to the buffer
            buffer[m++] = holderVector[cur].center[0];
            buffer[m++] = holderVector[cur].center[1];
            buffer[m++] = holderVector[cur].normal_y;
            buffer[m++] = holderVector[cur].cp;
            next = (i+1)%(holderVector.size());
            //Compare their x_loc,If the same 
            TheSame=true;
            while(TheSame)
            {
                if(fabsf64(holderVector[cur].center[0]-holderVector[next].center[0])>0.0001)
                {
                    TheSame =false;
                }
                else
                {
                    if(fabsf64(holderVector[next].normal_y)>fabsf64(holderVector[cur].normal_y))
                    {
                        buffer[m-1] = holderVector[next].cp;
                        buffer[m-2] = holderVector[next].normal_y;
                        cur = next;
                        next = next+1;
                    }
                    else
                    {
                        next = next+1;
                    }
                    i++;
                }
            }
        }
        int nCp = m/4;
        std::ofstream fout;
        fout.open("Cp_"+std::to_string(step));
        fout<<"Title = \"Dussin_exp\"\n";
        fout<<"VARIABLES = \"X\",\"Cp\"\n";
        fout<<"ZONE T=\"Dussin_exp\", I = "<<nCp+1<<", F=POINT\n";

        for(int i = 0;i<nCp;i++)
        {
            fout<<buffer[4*i]<<"\t"<<buffer[4*i+3]<<'\n';//<<"\t"<<buffer[4*i+2]<<"\t"<<buffer[4*i+3]<<'\n';
        }
        fout<<buffer[4*0]<<"\t"<<buffer[4*0+3]<<'\n';
        fout.close();

        

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
    delete[] wallNodeInds;
    delete[] p_average;
    delete[] count_average;
    if(nlocalWall>0)
    {
        delete[] realData;
        
    }
    delete[] nEachLocalWall;
}