#include "SamraiWrapper.h"
#include "SAMRAI/tbox/BalancedDepthFirstTree.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "memory"
#include "numeric"

WrapSAMRAI::WrapSAMRAI()
{
}
WrapSAMRAI::WrapSAMRAI(int argc, char *argv[])
{
    if ((argc != 2) && (argc != 4))
    {
        SAMRAI::tbox::pout << "USAGE: " << argv[0] << "<input filename>"
                           << "<restartdir><restore number>[options]\n"
                           << " options:\n"
                           << " none at this time"
                           << std::endl;
        SAMRAI::tbox::SAMRAI_MPI::abort();
    }
    else
    {
        input_filename = argv[1];
        if (argc == 4)
        {
            restart_read_dirname = argv[2];
            restore_num = atoi(argv[3]);
            is_from_restart = true;
        }
    }
}
void WrapSAMRAI::getFromInput(const std::shared_ptr<SAMRAI::tbox::Database> &input_db, bool is_from_restart)
{
    //Set dimension for tioga and reader.
    const SAMRAI::tbox::Dimension dim(static_cast<unsigned short>(input_db->getInteger("dim")));

    //Logging stuff for SAMRAI
    const std::string base_name = input_db->getStringWithDefault("base_name", " unnamed");
    const std::string log_filename = input_db->getStringWithDefault("log_filename", base_name + ".log");

    bool log_all_nodes = false;
    if (input_db->keyExists("log_all_nodes"))
    {
        log_all_nodes = input_db->getBool("log_all_nodes");
    }
    if (log_all_nodes)
    {
        SAMRAI::tbox::PIO::logAllNodes(log_filename);
    }
    else
    {
        SAMRAI::tbox::PIO::logOnlyNodeZero(log_filename);
    }

    int viz_dump_interval = 0;
    if (input_db->keyExists("viz_dump_interval"))
    {
        viz_dump_interval = input_db->getInteger("viz_dump_interval");
    }

    const std::string viz_dump_dirname =
        input_db->getStringWithDefault("viz_dump_dirname", base_name + ".visit");

    int visit_number_procs_per_file = 1;
    if (viz_dump_interval > 0)
    {
        if (input_db->keyExists("visit_number_procs_per_file"))
        {
            visit_number_procs_per_file =
                input_db->getInteger("visit_number_procs_per_file");
        }
    }

    viz_dump_data = (viz_dump_interval > 0);

    int restart_interval = 0;
    if (input_db->keyExists("restart_interval"))
    {
        restart_interval = input_db->getInteger("restart_interval");
    }

    const std::string restart_write_dirname =
        input_db->getStringWithDefault("restart_write_dirname", base_name + ".restart");
#if (TESTING == 1) && !defined(HAVE_HDF5)
    is_from_restart = false;
    restart_interval = 0;
#endif
    int BUFFER_SIZE = input_db->getInteger("buffer_size");

    const bool write_restart = (restart_interval > 0) && !(restart_write_dirname.empty());
    SAMRAI::tbox::RestartManager *restart_manager = SAMRAI::tbox::RestartManager::getManager();
    const SAMRAI::tbox::SAMRAI_MPI &mpi(SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld());
    if (is_from_restart)
    {
        restart_manager->openRestartFile(restart_read_dirname, restore_num, mpi.getSize());
    }

    //We need a struct here to pack all that SAMRAI needs
    grid_geometry = std::make_shared<SAMRAI::geom::CartesianGridGeometry>(
        dim,
        "CartesianGeometry",
        input_db->getDatabase("CartesianGeometry"));

    patch_hierarchy = std::make_shared<SAMRAI::hier::PatchHierarchy>(
        "PatchHierarchy",
        grid_geometry,
        input_db->getDatabase("PatchHierarchy"));

    geo_model = new PureGeoStrategy(
        "PureGeoStrategy",
        dim,
        input_db->getDatabase("PureGeoStrategy"),
        grid_geometry);

    geo_overset = new OversetStrategy(
        "OversetStrategy",
        dim,
        grid_geometry);

    geo_interp = new InterpStrategy(
        "InterpStrategy",
        dim,
        input_db->getDatabase("InterpStrategy"),
        grid_geometry);

    geo_integrator = std::make_shared<PureGeometricIntegrator>(
        "PureGeometircIntegrator",
        dim,
        input_db->getDatabase("PureGeometircIntegrator"),
        geo_model,
        geo_overset,
        geo_interp);

    error_detector = std::make_shared<SAMRAI::mesh::StandardTagAndInitialize>(
        "StandardTagAndInitialize",
        geo_integrator.get(),
        input_db->getDatabase("StandardTagAndInitialize"));
    box_generator = std::make_shared<SAMRAI::mesh::BergerRigoutsos>(
        dim,
        input_db->getDatabase("BergerRigoutsos"));
    load_balancer = std::make_shared<SAMRAI::mesh::TreeLoadBalancer>(
        dim,
        "LoadBalancer",
        input_db->getDatabase("LoadBalancer"),
        std::shared_ptr<SAMRAI::tbox::RankTreeStrategy>(new SAMRAI::tbox::BalancedDepthFirstTree));

    load_balancer->setSAMRAI_MPI(SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld());

    gridding_algorithm = std::make_shared<SAMRAI::mesh::GriddingAlgorithm>(
        patch_hierarchy,
        "GriddingAlgorithm",
        input_db->getDatabase("GriddingAlgorithm"),
        error_detector,
        box_generator,
        load_balancer);

#ifdef HAVE_HDF5
    visit_data_writer = std::make_shared<SAMRAI::appu::VisItDataWriter>(
        dim,
        "PureGeo Visit Writer",
        viz_dump_dirname,
        visit_number_procs_per_file);
    geo_model->registerVisItDataWriter(visit_data_writer);
    //geo_overset->registerVisItDataWriter(visit_data_writer);
#endif
    tec_writer = std::make_shared<TecplotWriter>(
        "IamBored",
        dim);
    tec_writer->init();
    tec_writer->attachHierarchy(patch_hierarchy);
    geo_model->registerTecplotWriter(tec_writer);
    geo_overset->registerTecplotWriter(tec_writer);

    SAMRAI::tbox::plog << "\n Check input data and variables before simulation:"
                       << std::endl;
    SAMRAI::tbox::plog << "Input database..." << std::endl;
    input_db->printClassData(SAMRAI::tbox::plog);
    SAMRAI::tbox::plog << "\n Variable database..." << std::endl;
    SAMRAI::hier::VariableDatabase::getDatabase()->printClassData(SAMRAI::tbox::plog);
    geo_integrator->initializeIntegrator(gridding_algorithm);

    tag_buffer_array.resize(patch_hierarchy->getMaxNumberOfLevels());
    regrid_start_time.resize(patch_hierarchy->getMaxNumberOfLevels());
    std::fill(tag_buffer_array.begin(), tag_buffer_array.end(), BUFFER_SIZE);
}

void WrapSAMRAI::init(SamraiFeeder *feeder)
{
    //Get geo_overset and geo_interp the oversetGrid data that they might be interested in.
    geo_overset->registerFeeder(feeder);
    geo_interp->registerFeeder(feeder);
    //Now as init() being called, the grid is
}

void WrapSAMRAI::make_geometric_refine()
{
    double loop_time = 0;
    int loop_circle = 0;
    gridding_algorithm->makeCoarsestLevel(loop_time);
    bool done = false;
    bool initial_circle = true;
    for (int ln = 0; patch_hierarchy->levelCanBeRefined(ln) and !done; ++ln)
    {
        gridding_algorithm->makeFinerLevel(
            tag_buffer_array[ln],
            initial_circle,
            loop_circle,
            loop_time);
        done = !(patch_hierarchy->finerLevelExists(ln));
        std::cout << ln << '\t' << done << '\t' << patch_hierarchy->levelCanBeRefined(ln) << '\n';
    }
}
void WrapSAMRAI::dump_data(double time)
{
    if (viz_dump_data)
    {
#ifdef HAVE_HDF5
        visit_data_writer->writePlotData(patch_hierarchy, 10, time);
#endif
    }
    //tec_writer->writeAll();
}

void WrapSAMRAI::pack_patch_for_overset(
    int *nf, int *qstride, double **qnodein,
    int **idata, double **rdata,
    int *ngridsin, int *qnodesize,
    int ***ibl, double ***q, int **block_global_id, int *nblockin)
{
    //int ln = geo_interp->geo_fine_level();
    for (int ln = 0; ln < patch_hierarchy->getNumberOfLevels(); ln++)
    {
        geo_integrator->findOversetCartBlocks(patch_hierarchy, ln, 0);
    }

    geo_integrator->extractCartGrid(
        nf, qstride, qnodein,
        idata, rdata,
        ngridsin, qnodesize,
        ibl, q, block_global_id, nblockin);
}

WrapSAMRAI::~WrapSAMRAI()
{
    // if (geo_model)
    //     delete[] geo_model;
    // if (geo_overset)
    //     delete[] geo_overset;
    // grid_geometry.reset();
    // patch_hierarchy.reset();
    // geo_integrator.reset();
    // error_detector.reset();
    // box_generator.reset();
    // load_balancer.reset();
    // gridding_algorithm.reset();
    // visit_data_writer.reset();
    tec_writer->finalize();
    tec_writer.reset();
    MPI_Barrier(MPI_COMM_WORLD);
}

void WrapSAMRAI::extractGrid(
    int *nf, int *qstride, double **qnodein,
    int **idata, double **rdata,
    int *ngridsout, int *qnodesize)
{
    int local_ngrids = patch_hierarchy->getNumberBlocks();
    //We need to specify an nf for all variables
    int dim = grid_geometry->getDim().getValue();
    //In CURRENT context, there is no ghost cell
    *nf = 0;
    //Number of Variables is set to be 1,We might need to make a depth related Variable.
    *qstride = 1;
    *qnodein = new double[3];
    //These should be put in the init() part
    *qnodesize = 3; //1+2
    *qnodein[0] = 0;
    *qnodein[1] = -1;
    *qnodein[2] = 1;
    int local_intSize = (2 * dim + 5) * local_ngrids;
    int *local_intData = new int[local_intSize];
    int local_realSize = (2 * dim) * local_ngrids;
    double *local_realData = new double[local_realSize];

    int num_proc = SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld().getSize();

    fillCartGridData(*idata, *rdata, local_ngrids);
    //*idata = new int[(2 * dim + 4) * (*ngrids)];
    //*rdata = new double[2 * dim * (*ngrids)];
    std::vector<int> localBlockCount(num_proc);
    MPI_Allgather(&local_ngrids, 1, MPI_INT, localBlockCount.data(), 1, MPI_INT, MPI_COMM_WORLD);
    int ntotalCarts = std::accumulate(localBlockCount.begin(), localBlockCount.end(), 0);
    int global_intSize = (2 * dim + 5) * ntotalCarts;
    int *global_intData = new int[global_intSize];
    int global_realSize = (2 * dim) * ntotalCarts;
    double *global_realData = new double[global_realSize];

    std::vector<int> displs(num_proc + 1);
    std::vector<int> iDataSizePerProc(num_proc);
    std::vector<int> rDataSizePerProc(num_proc);

    displs[0] = 0;
    iDataSizePerProc[0] = localBlockCount[0] * (2 * dim + 5);

    for (int i = 1; i <= num_proc; i++)
    {
        displs[i] = displs[i - 1] + iDataSizePerProc[i - 1];
        if (i < num_proc)
            iDataSizePerProc[i] = localBlockCount[i] * (2 * dim + 5);
    }
    MPI_Allgatherv(local_intData, local_intSize, MPI_INT, global_intData, iDataSizePerProc.data(), displs.data(), MPI_INT, MPI_COMM_WORLD);
    displs[0] = 0;
    rDataSizePerProc[0] = localBlockCount[0] * (2 * dim);
    for (int i = 1; i <= num_proc; i++)
    {
        displs[i] = displs[i - 1] + rDataSizePerProc[i - 1];
        if (i < num_proc)
            rDataSizePerProc[i] = localBlockCount[i] * (2 * dim);
    }
    MPI_Allgatherv(local_realData, local_realSize, MPI_DOUBLE, global_realData, rDataSizePerProc.data(), displs.data(), MPI_DOUBLE, MPI_COMM_WORLD);
    //We might want to do an ordering here?
    //Not now
    //
    *rdata = global_realData;
    *idata = global_intData;
    for (int i = 0; i < ntotalCarts; i++)
    {
        (*idata)[i * (2 * dim + 5)] = i;
    }
    delete[] local_intData;
    delete[] local_realData;
    delete[] global_intData;
    delete[] global_realData;
}

void WrapSAMRAI::fillCartGridData(int *idata, double *rdata, int ngrids)
{
    int dim = patch_hierarchy->getDim().getValue();
    int icount = 0;
    int rcount = 0;
    int local_id = 0;
    int cur_proc = SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld().getRank();

    for (int ln = 0; ln < patch_hierarchy->getMaxNumberOfLevels(); ln++)
    {
        std::shared_ptr<SAMRAI::hier::PatchLevel> level(
            patch_hierarchy->getPatchLevel(ln));
        for (SAMRAI::hier::PatchLevel::iterator ip(level->begin()); ip != level->end(); ip++)
        {
            const std::shared_ptr<SAMRAI::hier::Patch> &patch = *ip;
            const std::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> patch_geom(
                SAMRAI_SHARED_PTR_CAST<SAMRAI::geom::CartesianPatchGeometry, SAMRAI::hier::PatchGeometry>(patch->getPatchGeometry()));
            const double *dx = patch_geom->getDx();
            const double *xlo = patch_geom->getXLower();
            const double *xhi = patch_geom->getXUpper();

            const SAMRAI::hier::Box &pbox = patch->getBox();
            const SAMRAI::hier::Index ifirst = pbox.lower();
            const SAMRAI::hier::Index ilast = pbox.upper();
            icount++; //Leave a space for global id.
            idata[icount++] = ln;
            idata[icount++] = cur_proc;
            idata[icount++] = 0; //porder is set to 0 for now, later we might need to dynamically assign porder?//This needs another strategy
            idata[icount++] = local_id;
            for (int i = 0; i < dim; i++)
            {
                idata[icount++] = ifirst(i);
            }
            for (int i = 0; i < dim; i++)
            {
                idata[icount++] = ilast(i);
            }
            for (int i = 0; i < dim; i++)
            {
                rdata[rcount++] = xlo[i];
            }
            for (int i = 0; i < dim; i++)
            {
                rdata[rcount++] = dx[i];
            }
        }
    }
}