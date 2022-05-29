#include "samrai/PureGeoIntegrator.h"
//#include"SAMRAI/include/tbox/Dimension.h"
#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/tbox/RankTreeStrategy.h"
#include "SAMRAI/tbox/BalancedDepthFirstTree.h"

#include "UnstructSolver/UnstructIntegrator.h"
#include "UnstructSolver/Euler2D.h"
#include "UnstructSolver/ViscousFlow2D.h"
#include "UnstructSolver/RungeKuttaStrategy.h"
#include "UnstructSolver/LU_SGS_Strategy.h"
#include "reader/CobaltReader.h"
#include "Orchestra.h"
#include "UniTioga/tioga.h"
#include <fstream>
#include "reader/UnstructReader.h"
#include "toolBox/Timer.h"
void test_amr_index(int dim);
void test_list();
int tioga_main()
{
    // test_list();
    // return 0;
    // test_amr_index(2);
    // return 0;
    MPI_Init(NULL, NULL);
    int num_proc;
    int cur_proc;
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &cur_proc);
    std::cout << "curproc " << cur_proc << " numproc " << num_proc << '\n';
    CobaltReader *rd = new CobaltReader[1];
    // rd->addFile("part1.grd");
    // rd->addFile("part2.grd");
    rd->addFile("fuselage_part.grd");
    rd->addFile("rotor1.grd");
    rd->addFile("rotor2.grd");
    // rd->addFile("background.grd");
    // rd->addFile("fuselage.grd");
    // rd->addFile("rotor.grd");
    rd->setDim(3);
    rd->setComm(cur_proc, num_proc);
    // rd->setOverSign(-2147483647);
    // rd->setWallSign(-2);
    rd->setWallSign(-1);
    rd->setOverSign(-2);
    rd->readAll();
    // rd->writeFiles();
    TIOGA::tioga *tg = new TIOGA::tioga[1];
    tg->setDim(3);
    tg->setCommunicator(MPI_COMM_WORLD, cur_proc, num_proc);
    int mexclude = 1;
    tg->setMexclude(&mexclude);
    bool backgrounded = false;
    bool backgrounded2D = false;
    bool backgrounded3D = false;
    if (backgrounded)
    {
        if(backgrounded2D)
        {
        std::ofstream fout;
        std::string local_name = "CartInfo" + std::to_string(cur_proc);
        fout.open(local_name);
        int nf = 1;
        int qstride = 1;
        double *qnodein = new double[3];

        qnodein[0] = 0;
        qnodein[1] = -0.5;
        qnodein[2] = 0.5;
        int *idata = new int[4 * 9];
        double *rdata = new double[4 * 2 * 2];
        for (int i = 0; i < 4; i++)
        {
            idata[9 * i] = i;                            //global id
            idata[9 * i + 1] = 0;                        //level is 0
            idata[9 * i + 2] = i;                        //proc id
            idata[9 * i + 3] = 0;                        //porder set to 0;
            idata[9 * i + 4] = 0;                        //local id
            idata[9 * i + 5] = (i % 2 == 0 ? 1 : 161);   //ilow0
            idata[9 * i + 6] = (i % 3 == 0 ? 1 : 161);   //ilow1
            idata[9 * i + 7] = (i % 2 == 0 ? 160 : 320); //ihigh0
            idata[9 * i + 8] = (i % 3 == 0 ? 160 : 320); //ihigh1
        }

        for (int i = 0; i < 4; i++)
        {
            rdata[4 * i] = (i % 2 == 0 ? -4 : 0);
            rdata[4 * i + 1] = (i % 3 == 0 ? -4 : 0);
            rdata[4 * i + 2] = 0.025;
            rdata[4 * i + 3] = 0.025;
        }
        if (cur_proc == 0)
        {
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 9; j++)
                    fout << idata[9 * i + j] << '\t';
                fout << '\n';
                for (int j = 0; j < 4; j++)
                    fout << rdata[4 * i + j] << '\t';
                fout << '\n';
            }
        }
        int ngridsin = 4;
        int qnodesize = 3; //1+2=3
        int *q_ibl = new int[322 * 322];
        double *q_quan = new double[322 * 322];
        for (int i = 0; i < 322 * 322; i++)
        {
            q_ibl[i] = 1;
            q_quan[i] = 1;
        }
        tg->register_amr_global_data(nf, qstride, qnodein, idata, rdata, ngridsin, qnodesize);

        tg->set_amr_patch_count(1);

        tg->register_amr_local_data(0, cur_proc, q_ibl, q_quan);
        }
        if(backgrounded3D)
        {
        std::ofstream fout;
        std::string local_name = "CartInfo" + std::to_string(cur_proc);
        fout.open(local_name);
        int nf = 1;
        int qstride = 1;
        double *qnodein = new double[3];

        qnodein[0] = 0;
        qnodein[1] = -0.5;
        qnodein[2] = 0.5;
        int *idata = new int[4 * 11];
        double *rdata = new double[4 * 3 * 2];
        for (int i = 0; i < 4; i++)
        {
            idata[11 * i] = i;                            //global id
            idata[11 * i + 1] = 0;                        //level is 0
            idata[11 * i + 2] = i;                        //proc id
            idata[11 * i + 3] = 0;                        //porder set to 0;
            idata[11 * i + 4] = 0;                        //local id
            idata[11 * i + 5] = (i % 2 == 0 ? 1 : 301);   //ilow0
            idata[11 * i + 6] = (i % 3 == 0 ? 1 : 301);   //ilow1
            idata[11 * i + 7] = 1;                        //ilow2
            idata[11 * i + 8] = (i % 2 == 0 ? 301 : 601); //ihigh0
            idata[11 * i + 9] = (i % 3 == 0 ? 301 : 601); //ihigh1
            idata[11 * i + 10] = 1501;                     //ihigh2
        }

        for (int i = 0; i < 4; i++)
        {
            rdata[6 * i] = (i % 2 == 0 ?  -5 : 0);
            rdata[6 * i + 1] = (i % 3 == 0 ? -5 : 0);
            rdata[6 * i + 2] = -10;//-12 to 27 is 39
            rdata[6 * i + 3] = 1/6;
            rdata[6 * i + 4] = 1/6;
            rdata[6 * i + 5] = 1/6;
        }
        if (cur_proc == 0)
        {
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 9; j++)
                    fout << idata[9 * i + j] << '\t';
                fout << '\n';
                for (int j = 0; j < 4; j++)
                    fout << rdata[4 * i + j] << '\t';
                fout << '\n';
            }
        }
        int ngridsin = 4;
        int qnodesize = 3; //1+2=3
        int *q_ibl = new int[303 * 303* 1503];
        double *q_quan = new double[303 * 303 * 1503];
        for (int i = 0; i < 303 * 303 * 1503; i++)
        {
            q_ibl[i] = 1;
            q_quan[i] = 1;
        }
        tg->register_amr_global_data(nf, qstride, qnodein, idata, rdata, ngridsin, qnodesize);

        tg->set_amr_patch_count(1);

        tg->register_amr_local_data(0, cur_proc, q_ibl, q_quan);

        }
    }
    
    tg->set_distanced(false);
    tg->registerFromFeeder(rd);
    
    tg->profile();
    MPI_Finalize();
    return 1;
    tg->reResolution();
    
    tg->performConnectivity();
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "Connectivity Performed\n";
    if (backgrounded)
    {
        tg->performConnectivityAMR();
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // tg->performConnectivity();
    // delete[] tg;
    // delete[] rd;
    MPI_Finalize();
}

int samrai_main(int argc, char *argv[])
{
    //1:Iitialize MPI and SAMRAI.
    SAMRAI::tbox::SAMRAI_MPI::init(&argc, &argv);
    SAMRAI::tbox::SAMRAIManager::initialize();

    //2:StartUp SAMRAI.
    SAMRAI::tbox::SAMRAIManager::startup();

    const SAMRAI::tbox::SAMRAI_MPI &mpi(SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld());
    int num_proc;
    int cur_proc;
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &cur_proc);
    std::cout << "curproc " << cur_proc << " numproc " << num_proc << '\n';
    CobaltReader *rd = new CobaltReader[1];
    rd->addFile("part1.grd");
    rd->addFile("part2.grd");
    rd->setDim(2);
    rd->setComm(cur_proc, num_proc);
    rd->setOverSign(-2147483647);
    rd->setWallSign(-2);
    rd->readAll();
    //rd->writeFiles();
    //delete[] rd;

    int num_failures = 0;

    if (true)
    {
        std::string input_filename;
        std::string restart_read_dirname;
        int restore_num = 0;
        bool is_from_restart = false;
        if ((argc != 2) && (argc != 4))
        {
            SAMRAI::tbox::pout << "USAGE: " << argv[0] << "<input filename>"
                               << "<restartdir><restore number>[options]\n"
                               << " options:\n"
                               << " none at this time"
                               << std::endl;
            SAMRAI::tbox::SAMRAI_MPI::abort();
            return -1;
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

        SAMRAI::tbox::plog << "input_filename = " << input_filename << std::endl;
        SAMRAI::tbox::plog << "restart_read_dirname = " << restart_read_dirname << std::endl;
        SAMRAI::tbox::plog << "restore_num = " << restore_num << std::endl;

        //3.Create input database and parse input file
        std::shared_ptr<SAMRAI::tbox::InputDatabase> input_db(new SAMRAI::tbox::InputDatabase("input_db"));
        SAMRAI::tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

        //Reading the input database
        if (input_db->keyExists("GlobalInputs"))
        {
            std::shared_ptr<SAMRAI::tbox::Database> global_db(input_db->getDatabase("GlobalInputs"));
            if (global_db->keyExists("call_abort_in_serial_instead_of_exit"))
            {
            }
            bool flag = global_db->getBool("call_abort_in_serial_instead_of_exit");
            SAMRAI::tbox::SAMRAI_MPI::setCallAbortInParallelInsteadOfMPIAbort(flag);
        }

        std::shared_ptr<SAMRAI::tbox::Database> main_db(input_db->getDatabase("Main"));

        const SAMRAI::tbox::Dimension dim(static_cast<unsigned short>(main_db->getInteger("dim")));
        const std::string base_name = main_db->getStringWithDefault("base_name", " unnamed");
        const std::string log_filename = main_db->getStringWithDefault("log_filename", base_name + ".log");

        bool log_all_nodes = false;
        if (main_db->keyExists("log_all_nodes"))
        {
            log_all_nodes = main_db->getBool("log_all_nodes");
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
        if (main_db->keyExists("viz_dump_interval"))
        {
            viz_dump_interval = main_db->getInteger("viz_dump_interval");
        }

        const std::string viz_dump_dirname =
            main_db->getStringWithDefault("viz_dump_dirname", base_name + ".visit");

        int visit_number_procs_per_file = 1;
        if (viz_dump_interval > 0)
        {
            if (main_db->keyExists("visit_number_procs_per_file"))
            {
                visit_number_procs_per_file =
                    main_db->getInteger("visit_number_procs_per_file");
            }
        }

        const bool viz_dump_data = (viz_dump_interval > 0);

        int restart_interval = 0;
        if (main_db->keyExists("restart_interval"))
        {
            restart_interval = main_db->getInteger("restart_interval");
        }

        const std::string restart_write_dirname =
            main_db->getStringWithDefault("restart_write_dirname", base_name + ".restart");
#if (TESTING == 1) && !defined(HAVE_HDF5)
        is_from_restart = false;
        restart_interval = 0;
#endif
        int BUFFER_SIZE = main_db->getInteger("buffer_size");
        const bool write_restart = (restart_interval > 0) && !(restart_write_dirname.empty());

        SAMRAI::tbox::RestartManager *restart_manager = SAMRAI::tbox::RestartManager::getManager();
        if (is_from_restart)
        {
            restart_manager->openRestartFile(restart_read_dirname, restore_num, mpi.getSize());
        }
        //There suppose to be a recorder for main but we ignore it for now

        //
        std::shared_ptr<SAMRAI::geom::CartesianGridGeometry> grid_geometry(
            new SAMRAI::geom::CartesianGridGeometry(
                dim,
                "CartesianGeometry",
                input_db->getDatabase("CartesianGeometry")));

        std::shared_ptr<SAMRAI::hier::PatchHierarchy> patch_hierarchy(
            new SAMRAI::hier::PatchHierarchy(
                "PatchHierarchy",
                grid_geometry,
                input_db->getDatabase("PatchHierarchy")));

        PureGeoStrategy *geo_model = new PureGeoStrategy(
            "PureGeoStrategy",
            dim,
            input_db->getDatabase("PureGeoStrategy"),
            grid_geometry);
        OversetStrategy *geo_overset = new OversetStrategy(
            "OversetStrategy",
            dim,
            grid_geometry,
            rd);

        //geo_model->set_geometry(rd->files.size(),)

        std::shared_ptr<PureGeometricIntegrator> geo_integrator(
            new PureGeometricIntegrator(
                "PureGeometircIntegrator",
                dim,
                input_db->getDatabase("PureGeometircIntegrator"),
                geo_model,
                geo_overset));
        //For debugging
        geo_integrator->set_max_level(patch_hierarchy->getMaxNumberOfLevels());

        std::shared_ptr<SAMRAI::mesh::StandardTagAndInitialize> error_detector(
            new SAMRAI::mesh::StandardTagAndInitialize(
                "StandardTagAndInitialize",
                geo_integrator.get(),
                input_db->getDatabase("StandardTagAndInitialize")));
        std::shared_ptr<SAMRAI::mesh::BergerRigoutsos> box_generator(
            new SAMRAI::mesh::BergerRigoutsos(
                dim,
                input_db->getDatabase("BergerRigoutsos")));
        std::shared_ptr<SAMRAI::mesh::TreeLoadBalancer> load_balancer(
            new SAMRAI::mesh::TreeLoadBalancer(
                dim,
                "LoadBalancer",
                input_db->getDatabase("LoadBalancer"),
                std::shared_ptr<SAMRAI::tbox::RankTreeStrategy>(new SAMRAI::tbox::BalancedDepthFirstTree)));

        load_balancer->setSAMRAI_MPI(SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld());

        std::shared_ptr<SAMRAI::mesh::GriddingAlgorithm> gridding_algorithm(
            new SAMRAI::mesh::GriddingAlgorithm(
                patch_hierarchy,
                "GriddingAlgorithm",
                input_db->getDatabase("GriddingAlgorithm"),
                error_detector,
                box_generator,
                load_balancer));
#ifdef HAVE_HDF5
        std::shared_ptr<SAMRAI::appu::VisItDataWriter> visit_data_writer(
            new SAMRAI::appu::VisItDataWriter(
                dim,
                "PureGeo Visit Writer",
                viz_dump_dirname,
                visit_number_procs_per_file));
        geo_model->registerVisItDataWriter(visit_data_writer);
        geo_overset->registerVisItDataWriter(visit_data_writer);
#endif
        std::shared_ptr<TecplotWriter> tec_writer(
            new TecplotWriter(
                "IamBored",
                dim));
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

        //Now let's do some regrid

        std::vector<int> tag_buffer_array(patch_hierarchy->getMaxNumberOfLevels(), BUFFER_SIZE);
        std::vector<double> regrid_start_time(
            patch_hierarchy->getMaxNumberOfLevels());
        double loop_time = 0;
        int loop_circle = 0;
        if (SAMRAI::tbox::RestartManager::getManager()->isFromRestart())
        {
        }
        else
        {
            std::cout << "Not from restart\n";
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
        //End of restart if there is one, close the restart manager
        SAMRAI::tbox::RestartManager::getManager()->closeRestartFile();
        //std::cout<<patch_hierarchy->getNumberOfLevels()<<'\n';
        gridding_algorithm->regridAllFinerLevels(0,
                                                 tag_buffer_array,
                                                 0, 2, regrid_start_time);

        geo_integrator->checking(patch_hierarchy);
        MPI_Barrier(MPI_COMM_WORLD);
        std::cout << "ANother one bites the dust\n";
        geo_integrator->checking(patch_hierarchy);
        if (viz_dump_data)
        {
#ifdef HAVE_HDF5
            visit_data_writer->writePlotData(patch_hierarchy, 2, 0);
#endif
        }
        tec_writer->writeAll();
        std::cin.get();
        //End of the program, deallocating

        box_generator.reset();
        load_balancer.reset();
        gridding_algorithm.reset();
#ifdef HAVE_HDF5
        visit_data_writer.reset();
#endif
        tec_writer->finalize();
        error_detector.reset();
        geo_integrator.reset();
        if (geo_model)
            delete geo_model;
        patch_hierarchy.reset();
        grid_geometry.reset();

        main_db.reset();
        input_db.reset();
    }
    if (num_failures == 0)
    {
        SAMRAI::tbox::pout << "\n PASSED" << std::endl;
    }
    SAMRAI::tbox::SAMRAIManager::shutdown();
    SAMRAI::tbox::SAMRAIManager::finalize();
    SAMRAI::tbox::SAMRAI_MPI::finalize();
    return num_failures;
}

void unstruct_serial_v()
{
    UnstructReader* rd = new UnstructReader[1];
    // rd->setOverSign(-2);
    // rd->setWallSign(-1);
    // rd->addFile("naca0012.grd");
    rd->setOverSign(-2);
    rd->setWallSign(-1);
    rd->addFile("mix_element_laminar.grd");
    rd->setDim(2);
    std::cin.get();
    rd->readAll();

    std::cout<<"ReadAll done\n";
    std::cin.get();
    rd->writeFiles();
    int* nc;
    int** indC;
    rd->takeBoundary(-1,&nc,&indC);
    UnstructTopologyHolder* integrator = new UnstructTopologyHolder(rd);
    std::cout<<"Holder initailization done\n";
    std::cin.get();
    ViscousFlow2D* flower = new ViscousFlow2D(integrator);
    flower->test_unwantedSweep(indC,nc);
    std::cout<<"Flower initialization done\n";
    std::cin.get();
    //RungeKuttta
    bool use_implicit = true;
    if(!use_implicit)
    {
        RungeKuttaStrategy* rk_integrator = new RungeKuttaStrategy(integrator,flower,5);
        rk_integrator->initialize();
        for(int i = 0;i<40001;i++)
        {
            rk_integrator->singleStepSerial(i);
            if(i%5000==0)
            {
                std::cout<<i<<"\n";
                flower->writeCellData("cellDataSerial_"+std::to_string(i),0,0,flower->getU());
            }
        }
    }
    else
    {
        LUSGSStrategy* lusgs_integrator = new LUSGSStrategy(integrator,flower);
        lusgs_integrator->initialize();
        for(int i = 0;i<5001;i++)
        {
            lusgs_integrator->singleStepSerial(i);
            // std::cout<<i<<"Implicit Serial\n";
            // // if(i%500==0)
            // // {
            //     // std::cout<<i<<"Implicit Serial"<<'\n';
            //     flower->writeCellData("cellDataSerial_LUSGS_"+std::to_string(i),0,0,flower->getU());
            // // }
        }
    }
}

void unstruct_parallelv()
{
    MPI_Init(NULL,NULL);
    int num_proc;
    int cur_proc;
    MPI_Comm_size(MPI_COMM_WORLD,&num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD,&cur_proc);
    CobaltReader *rd = new CobaltReader[1];
    rd->setOverSign(-2);
    rd->setWallSign(-1);
    // rd->addFile("1.grd");
    rd->addFile("mix_element_laminar.grd");
    // rd->addFile("naca0012.grd");
    rd->setDim(2);
    rd->setComm(cur_proc,num_proc);
    // rd->setOverSign(-2);
    // rd->setWallSign(-1);
    //rd->setOverSign(-2147483647);
    //rd->setWallSign(-2);
    rd->readAll();

    
    
    
    MPI_Barrier(MPI_COMM_WORLD);
    rd->writeFiles();
    MPI_Barrier(MPI_COMM_WORLD);
    rd->performCommunication();
    UnstructTopologyHolder* integrator = new UnstructTopologyHolder(rd);
    integrator->arrange_communication();
    MPI_Barrier(MPI_COMM_WORLD);
    integrator->write("Meta");
    //
    //std::cout<<integrator->relatedProcs.size()<<"--------------------------"<<'\n';
    ViscousFlow2D* flower = new ViscousFlow2D(integrator);

    flower->test_communication();
    flower->test_partialcomm();
    
    // flower->outPutNondim("NOndim",cur_proc);
    
    bool use_implicit = true;
    if(!use_implicit)
    {
        RungeKuttaStrategy* rk_integrator = new RungeKuttaStrategy(integrator,flower,5);
        rk_integrator->initialize();
        simpleTimer timeIs;
        for(int i = 0;i<10;i++)
        {
            
            rk_integrator->singleStep(i);
            if(i%500==0)
            {

            }
            if(i%5000==0)
            {
                std::cout<<i<<"\n";
                flower->AllCellCommunication(flower->getU());
                flower->writeCellData("cellDatParallel"+std::to_string(i+1),cur_proc,0,flower->getU());   
            }
        }
        std::string Cp_name = "OneAndOnlyLegendaryCp.h5";
        flower->outPutCp(Cp_name,0,10000);
    }
    else
    {
        simpleTimer timeIs;
        LUSGSStrategy* lusgs = new LUSGSStrategy(integrator,flower);
        lusgs->initialize();
        for(int i = 0;i<1000001;i++)
        {
            lusgs->singleStep(i);
            if(i%5000==0)
            {
                std::cout<<i<<'\n';
                lusgs->postprocessUpdate();
                flower->AllCellCommunication(flower->getU());
                flower->writeCellData("cellDataParallelLUSGS"+std::to_string(i+1),cur_proc,0,flower->getU());
                std::string Cp_name = "OneAndOnlyLegendaryCp.h5";
                flower->outPutCp(Cp_name,0,i);
            }
        }
        std::cout<<"Finished?\n";
        
        
        
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

}

void unstruct_serial()
{
    UnstructReader* rd = new UnstructReader[1];
    // rd->setOverSign(-2);
    // rd->setWallSign(-1);
    // rd->addFile("naca0012.grd");
    rd->setOverSign(-1);
    rd->setWallSign(-4);
    rd->addFile("rae2822_euler.grd");
    rd->setDim(2);
    std::cin.get();
    rd->readAll();

    std::cout<<"ReadAll done\n";
    std::cin.get();
    rd->writeFiles();
    int* nc;
    int** indC;
    rd->takeBoundary(-1,&nc,&indC);
    UnstructTopologyHolder* integrator = new UnstructTopologyHolder(rd);
    std::cout<<"Holder initailization done\n";
    std::cin.get();
    Euler2D* euler = new Euler2D(integrator);
    euler->test_unwantedSweep(indC,nc);
    std::cout<<"Euler initialization done\n";
    std::cin.get();
    //RungeKuttta
    bool use_implicit = true;
    if(!use_implicit)
    {
        RungeKuttaStrategy* rk_integrator = new RungeKuttaStrategy(integrator,euler,5);
        rk_integrator->initialize();
        for(int i = 0;i<20000;i++)
        {
            rk_integrator->singleStepSerial(i);
            if(i%500==0)
            {
                std::cout<<i<<"\n";
                euler->writeCellData("cellDataSerial_"+std::to_string(i),0,0,euler->getU());
            }
        }
    }
    else
    {
        LUSGSStrategy* lusgs_integrator = new LUSGSStrategy(integrator,euler);
        lusgs_integrator->initialize();
        lusgs_integrator->checkEdgeCountSerial("Serial_edgeCount");
        simpleTimer timeIs;
        for(int i = 0;i<10001;i++)
        {
            lusgs_integrator->singleStepSerial(i);
            if(i%500==0)
            {
                std::cout<<i<<"Using Implicit"<<'\n';
                euler->writeCellData("cellDataSerial_LUSGSOn"+std::to_string(i+1),0,0,euler->getU());
            }
        }


    }
    std::cout<<"MPI BARRIER?\n";
}
    

void unstruct_main()
{
    MPI_Init(NULL,NULL);
    int num_proc;
    int cur_proc;
    MPI_Comm_size(MPI_COMM_WORLD,&num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD,&cur_proc);
    CobaltReader *rd = new CobaltReader[1];
    // rd->addFile("naca0012.grd");
    // rd->setOverSign(-2);
    // rd->setWallSign(-1);
    

    rd->setOverSign(-1);
    rd->setWallSign(-4);
    rd->addFile("rae2822_euler.grd");

    // rd->addFile("naca0012.grd");
    rd->setDim(2);
    rd->setComm(cur_proc,num_proc);
    //rd->setOverSign(-2147483647);
    //rd->setWallSign(-2);
    rd->readAll();

    
    
    
    MPI_Barrier(MPI_COMM_WORLD);
    rd->writeFiles();
    MPI_Barrier(MPI_COMM_WORLD);
    rd->performCommunication();
    UnstructTopologyHolder* integrator = new UnstructTopologyHolder(rd);
    integrator->arrange_communication();
    MPI_Barrier(MPI_COMM_WORLD);
    integrator->write("Meta");
    //
    //std::cout<<integrator->relatedProcs.size()<<"--------------------------"<<'\n';
    Euler2D* euler = new Euler2D(integrator);

    euler->test_communication();
    euler->test_partialcomm();
    
    // euler->outPutNondim("NOndim",cur_proc);
    
    bool use_implicit = true;
    if(!use_implicit)
    {
        RungeKuttaStrategy* rk_integrator = new RungeKuttaStrategy(integrator,euler,5);
        rk_integrator->initialize();
        
        for(int i = 0;i<20001;i++)
        {
            rk_integrator->singleStep(i);
            if(i%5000==0)
            {
                std::cout<<i<<"\n";
            

            }
                
            if(i%5000==0)
            {
                euler->AllCellCommunication(euler->getU());
                euler->writeCellData("cellDatParallel"+std::to_string(i+1),cur_proc,0,euler->getU());
                std::string Cp_name = "OneAndOnlyLegendaryCp.h5";
                euler->outPutCp(Cp_name,0,i);
            }
        }
        
    }
    else
    {
        simpleTimer timeIs;
        LUSGSStrategy* lusgs = new LUSGSStrategy(integrator,euler);
        lusgs->initialize();
        for(int i = 0;i<5001;i++)
        {
            lusgs->singleStep(i);
            if(i%500==0)
            {
                std::cout<<i<<'\n';
                euler->AllCellCommunication(euler->getU());
                euler->writeCellData("cellDataParallelLUSGS"+std::to_string(i+1),cur_proc,0,euler->getU());
                std::string Cp_name = "OneAndOnlyLegendaryCp.h5";
                euler->outPutCp(Cp_name,0,i);
            }
        }
        std::cout<<"Finished?\n";
        
        
    }
    std::cout<<"MPI BARRIER?\n";
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout<<"MPI BARRIERED\n";
    MPI_Finalize();
    std::cout<<"MPI FINALIZED?\n";
}

void orchestra_main(int argc, char *argv[])
{
    Orchestra *orchestra = new Orchestra(argc, argv);

    orchestra->init();
    MPI_Barrier(MPI_COMM_WORLD);
    orchestra->finalize();

    //delete orchestra;
}
int main(int argc, char *argv[])
{
    bool tioga = false;
    if (tioga)
    {
        tioga_main();
    }
    bool samrai = false;
    if (samrai)
    {
        samrai_main(argc, argv);
    }
    bool orchestra_used = false;
    if (orchestra_used)
    {
        orchestra_main(argc, argv);
    }
    bool unstruct_used = true;
    bool serial_used = false;
    
    if(unstruct_used)
    {
        if(serial_used)
        {
            unstruct_serial_v();
        }
        else
        {
            unstruct_parallelv();
        }
    }
}