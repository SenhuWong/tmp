#include "InterpStrategy.h"
#include "PureGeoIntegrator.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include <numeric>

InterpStrategy::InterpStrategy(const std::string &object_name,
                               const SAMRAI::tbox::Dimension &dim,
                               std::shared_ptr<SAMRAI::tbox::Database> input_db,
                               std::shared_ptr<SAMRAI::geom::CartesianGridGeometry> grid_geometry)
    : d_object_name(object_name),
      d_dim(dim),
      d_grid_geometry(grid_geometry),
      d_allocator(SAMRAI::tbox::AllocatorDatabase::getDatabase()->getDefaultAllocator())

{
    getFromInput(input_db, false);
    const SAMRAI::tbox::SAMRAI_MPI &mpi = SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld();
    cur_proc = mpi.getRank();
    num_proc = mpi.getSize();
}
void InterpStrategy::getFromInput(const std::shared_ptr<SAMRAI::tbox::Database> &input_db, bool is_from_restart)
{
    least_refine = input_db->getInteger("geo_fine_level");
    qnodesize = input_db->getInteger("qnodeSize");
    qnode = new double[qnodesize];
    input_db->getDoubleArray("qnode", qnode, qnodesize);
}

void InterpStrategy::registerModelVariables(PureGeometricIntegrator *integrator)
{
    std::list<std::shared_ptr<SAMRAI::hier::Variable>> &ZipVariables = integrator->d_field_variable;
    //This requires that PureGeoStrategy is initialized first.
    ZipDepth = ZipVariables.size();
    ZipIndexes = new int[ZipDepth];
    SAMRAI::hier::VariableDatabase *var_db = SAMRAI::hier::VariableDatabase::getDatabase();
    int count = 0;
    for (std::list<std::shared_ptr<SAMRAI::hier::Variable>>::iterator iter(ZipVariables.begin()); iter != ZipVariables.end(); iter++)
    {
        int index = var_db->mapVariableAndContextToIndex(*iter, workspace);
        ZipIndexes[count++] = index;
    }
    //Now we get all the indexes of the variables that we need to Zip;
    zip_variable = std::make_shared<SAMRAI::pdat::CellVariable<double>>(d_dim, "V_dot_Zip", ZipDepth);
    integrator->registerVariable(
        zip_variable,
        SAMRAI::hier::IntVector::getZero(d_dim),
        d_grid_geometry,
        "NO_COARSEN",
        "NO_REFINE",
        PureGeometricIntegrator::CURRENT_ONLY);
}
//For now there is just copy and paste
double *InterpStrategy::make_zip(SAMRAI::hier::Patch &patch, double zip_data_time)
{
    //First allocate all the data that is needed
    SAMRAI::hier::VariableDatabase *var_db = SAMRAI::hier::VariableDatabase::getDatabase();
    int zip_ind = var_db->mapVariableAndContextToIndex(zip_variable, workspace);
    patch.allocatePatchData(zip_ind);
    std::shared_ptr<SAMRAI::pdat::CellData<double>> zip_data(
        SAMRAI_SHARED_PTR_CAST<SAMRAI::pdat::CellData<double>, SAMRAI::hier::PatchData>(patch.getPatchData(zip_ind)));
    int count = 0;
    const SAMRAI::hier::Box &pbox = patch.getBox();
    int box_size = pbox.size();
    for (int i = 0; i < ZipDepth; i++)
    {
        double *dep_i = zip_data->getPointer(i);
        std::shared_ptr<SAMRAI::pdat::CellData<double>> to_zip_data(
            SAMRAI_SHARED_PTR_CAST<SAMRAI::pdat::CellData<double>, SAMRAI::hier::PatchData>(patch.getPatchData(ZipIndexes[i])));
        double *data_i = to_zip_data->getPointer();

        for (int j = 0; j < box_size; j++)
        {
            dep_i[j] = data_i[j];
        }
    }
    return zip_data->getPointer();
    //Shall I pack all these{q,ibl} and {idata,rdata} inside
}

int *InterpStrategy::make_ibl(SAMRAI::hier::Patch &patch, PureGeometricIntegrator *integrator, double make_block_time)
{

    std::shared_ptr<SAMRAI::pdat::CellVariable<int>> ibl_var = integrator->d_patch_strategy->d_iblank;
    SAMRAI::hier::VariableDatabase *var_db = SAMRAI::hier::VariableDatabase::getDatabase();
    int ibl_ind = var_db->mapVariableAndContextToIndex(ibl_var, workspace);
    patch.allocatePatchData(ibl_ind);
    std::shared_ptr<SAMRAI::pdat::CellData<int>> ibl_data(
        SAMRAI_SHARED_PTR_CAST<SAMRAI::pdat::CellData<int>, SAMRAI::hier::PatchData>(patch.getPatchData(ibl_ind)));
    int *ibl_ptr = ibl_data->getPointer();
    return ibl_ptr;
}
void InterpStrategy::formCartBlock(SAMRAI::hier::Patch &patch, PureGeometricIntegrator *integrator, double make_block_time)
{
    double *q = make_zip(patch, make_block_time);
    int *ibl = make_ibl(patch, integrator, make_block_time);
    const SAMRAI::hier::Box &pbox = patch.getBox();
    std::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> patch_geom(
        SAMRAI_SHARED_PTR_CAST<SAMRAI::geom::CartesianPatchGeometry, SAMRAI::hier::PatchGeometry>(patch.getPatchGeometry()));
    const SAMRAI::hier::Index &ifirst = pbox.lower();
    const SAMRAI::hier::Index &ilast = pbox.upper();
    const double *dx = patch_geom->getDx();
    const double *xlo = patch_geom->getXLower();

    related_patches.push_back(DataPack());
    DataPack &p = related_patches.back();
    p.idata[0] = -1;
    p.idata[1] = patch.getPatchLevelNumber();
    p.idata[2] = cur_proc;
    p.idata[3] = 0; //This could be replaced by a function so that porder is dynamically decided
    p.idata[4] = related_patches.size() - 1;
    for (int i = 0; i < d_dim.getValue(); i++)
    {
        p.idata[5 + i] = ifirst(i);
        p.idata[5 + d_dim.getValue() + i] = ilast(i);
        p.rdata[i] = xlo[i];
        p.rdata[i + d_dim.getValue()] = dx[i];
    }
    p.ibl = ibl;
    p.q = q;
    
}

void InterpStrategy::formCartGrid(int *nf, int *qstride, double **qnodein,
 int *ngridsin, int *qnodsize, int **idata, double **rdata,int**block_global_id)
{
    int loc_ngrids = related_patches.size();
    int dim = d_dim.getValue();

    int loc_intSize = (2 * dim + 5) * loc_ngrids;
    int loc_realSize = (2 * dim) * loc_ngrids;
    int *loc_intData = new int[loc_intSize];
    double *loc_realData = new double[loc_realSize];

    std::vector<int> nBlocksPerProc(num_proc);
    MPI_Allgather(&loc_ngrids, 1, MPI_INT, nBlocksPerProc.data(), 1, MPI_INT, MPI_COMM_WORLD);

    int ntotalBlocks = std::accumulate(nBlocksPerProc.begin(), nBlocksPerProc.end(), 0);

    int glo_intSize = (2 * dim + 5) * ntotalBlocks;
    int glo_realSize = (2 * dim) * ntotalBlocks;
    int *glo_intData = new int[glo_intSize];
    double *glo_realData = new double[glo_realSize];

    //Now fill in the
    int icount = 0;
    int rcount = 0;
    for (auto iter = related_patches.begin(); iter != related_patches.end(); iter++)
    {
        for (int i = 0; i < (2 * dim + 5); i++)
        {
            loc_intData[icount++] = iter->idata[i];
            if (cur_proc == 0)
            {
                std::cout << iter->idata[i] << '\t';
            }
        }
        if (cur_proc == 0)
        {
            std::cout << '\n';
        }

        for (int i = 0; i < 2 * dim; i++)
        {
            loc_realData[rcount++] = iter->rdata[i];
        }
    }
    std::vector<int> displs(num_proc + 1);
    std::vector<int> dataSizePerProc(num_proc);
    displs[0] = 0;
    dataSizePerProc[0] = (2 * dim + 5) * nBlocksPerProc[0];
    for (int i = 1; i < num_proc; i++)
    {
        displs[i] = displs[i - 1] + dataSizePerProc[i-1];
        
        if (i < num_proc)
            dataSizePerProc[i] = (2 * dim + 5) * nBlocksPerProc[i];
    }
    
    MPI_Allgatherv(loc_intData, loc_intSize, MPI_INT, glo_intData, dataSizePerProc.data(), displs.data(), MPI_INT, MPI_COMM_WORLD);
    
    //MPI_Barrier(MPI_COMM_WORLD);
    
    displs[0] = 0;
    dataSizePerProc[0] = (2 * dim) * nBlocksPerProc[0];
    for (int i = 1; i < num_proc; i++)
    {
        displs[i] = displs[i - 1] + dataSizePerProc[i - 1];
        if (i < num_proc)
            dataSizePerProc[i] = (2 * dim) * nBlocksPerProc[i];
    }
    MPI_Allgatherv(loc_realData, loc_realSize, MPI_DOUBLE, glo_realData, dataSizePerProc.data(), displs.data(), MPI_DOUBLE, MPI_COMM_WORLD);
    //MPI_Barrier(MPI_COMM_WORLD);
    //Now we are free to provide this to other API.
    int count = 0;
    for (int i = 0; i < ntotalBlocks; i++)
    {
        glo_intData[i * (2 * dim + 5)] = i;
        if(glo_intData[i * (2 * dim + 5)+2]==cur_proc)
        {
            (*block_global_id)[count++] = i;
        }
    }

    if (cur_proc == 0)
    {
        std::cout<<"GLobal"<<'\n';
        for (int i = 0; i < ntotalBlocks; i++)
        {
            for (int j = 0; j < 2 * dim + 5; j++)
            {
                std::cout << glo_intData[(2 * dim + 5) * i + j]<< '\t' ;
            }
            std::cout << '\n';
        }
    }


    *nf = 0;
    *qstride = ZipDepth;
    *qnodein = qnode;
    *ngridsin = ntotalBlocks;
    *idata = glo_intData;
    *rdata = glo_realData;
    *qnodsize = qnodesize;

    //I am not sure whether this will be needed somewhere else.
    delete[] loc_intData;
    delete[] loc_realData;
}

void InterpStrategy::takeCartBlock(double ***qs, int ***ibls, int *nblocks)
{
    *nblocks = related_patches.size();
    *qs = new double *[*nblocks];
    *ibls = new int *[*nblocks];
    int count = 0;
    for (auto iter = related_patches.begin(); iter != related_patches.end(); iter++)
    {
        (*ibls)[count] = iter->ibl;
        (*qs)[count] = iter->q;
        count++;
    }
}