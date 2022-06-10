#include "Orchestra.h"
//#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/BalancedDepthFirstTree.h"

#include <fstream>
//In the constructor, reader, samrai and tioga's init is done and
Orchestra::Orchestra(int argc, char *argv[])
{
    //Do SAMRAI environment initialize first
    SAMRAI::tbox::SAMRAI_MPI::init(&argc, &argv);
    SAMRAI::tbox::SAMRAIManager::initialize();
    SAMRAI::tbox::SAMRAIManager::startup();
    const SAMRAI::tbox::SAMRAI_MPI &mpi(SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld());
    num_proc = mpi.getSize();
    cur_proc = mpi.getRank();

    //Create 3 main part of orchestra
    rd = new CobaltReader[1];
    tg = new TIOGA::tioga[1];
    samrai = new WrapSAMRAI[1];
    solver = new WrapSolver[1];

    //Get input to initialize SAMRAI Strategy together with rd and tioga

    std::string input_filename;
    std::string restart_read_dirname;
    int restore_num = 0;
    bool is_from_restart = false;
    //For example,in a command like mpirun -np 4 ./myAPP test.input
    //argc is 2 [./myAPP] [test.input]
    // in a command like mpirun -np 4 ./myAPP test.input restart.data 8
    //argc is 4 [./myAPP] [test.input] [restart.data] [8]
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

    //3.Create input database and parse input file
    std::shared_ptr<SAMRAI::tbox::InputDatabase> input_db(new SAMRAI::tbox::InputDatabase("input_db"));
    //std::cout << input_filename << '\n';
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

    //Initialize for READER
    rd->setComm(cur_proc, num_proc);
    tg->setCommunicator(MPI_COMM_WORLD, cur_proc, num_proc);
    int dim;
    if (input_db->keyExists("Reader"))
    {
        std::shared_ptr<SAMRAI::tbox::Database> reader_db(input_db->getDatabase("Reader"));
        dim = reader_db->getInteger("dim");

        int nmesh = reader_db->getInteger("nmesh");
        int overTag = reader_db->getInteger("overTag");
        int wallTag = reader_db->getInteger("wallTag");
        std::string mesh_n;
        std::string mesh_name_n;
        for (int i = 1; i < nmesh + 1; i++)
        {
            mesh_n = "mesh_" + std::to_string(i);
            mesh_name_n = reader_db->getString(mesh_n);
            std::cout << mesh_name_n << '\n';
            rd->addFile(mesh_name_n);
        }
        rd->setDim(dim);
        rd->setOverSign(overTag);
        rd->setWallSign(wallTag);
        reader_db.reset();
        //rd->readALl() is moved into at the init() part
    }
    //Initialize for TIOGA
    if (input_db->keyExists("TIOGA"))
    {
        std::shared_ptr<SAMRAI::tbox::Database> tioga_db(input_db->getDatabase("TIOGA"));
        int dim = tioga_db->getInteger("dim");
        tg->setDim(dim);
        tg->setCommunicator(MPI_COMM_WORLD,cur_proc,num_proc);
        int mexclude = tioga_db->getInteger("mexclude");
        int nfringe = tioga_db->getInteger("nfringe");
        bool uses_weight = tioga_db->getBool("using_weight");
        tg->setMexclude(&mexclude);
        tg->setNfringe(&nfringe);
        tg->set_distanced(false);

        tioga_db.reset();
    }

    if (input_db->keyExists("SAMRAI"))
    {
        std::shared_ptr<SAMRAI::tbox::Database> samrai_db(input_db->getDatabase("SAMRAI"));
        samrai->getFromInput(samrai_db, false);
        samrai_db.reset();
    }

    if (input_db->keyExists("SOLVER"))
    {
        std::cout<<"Solver exist\n";
        std::shared_ptr<SAMRAI::tbox::Database> solver_db(input_db->getDatabase("SOLVER"));
        solver->setTopology(rd);
        
        std::string flow_type = solver_db->getString("FlowType");
        
        std::string advance_type = solver_db->getString("AdvanceType");
        
        std::string solve_type = solver_db->getString("SolveType");

        int nMaxStep = solver_db->getInteger("MaximumTimeStep");
        double minimalResidual = solver_db->getDouble("MinimumResidual");
        double ma;
        double re;
        double aoa;
        std::shared_ptr<SAMRAI::tbox::Database> flowParam_db(solver_db->getDatabase("Params"));
        {
            ma = flowParam_db->getDouble("Ma");
            re =  flowParam_db->getDouble("Re");
            if(dim==2)
            {
                aoa = flowParam_db->getDouble("AOA");//2D for now.
            }
            else 
            {
                throw std::runtime_error("3D not implemented yet\n");
            }
        }
        solver->setParams(flow_type,advance_type,solve_type,nMaxStep,minimalResidual,ma,re,&aoa);
        solver_db.reset();
    }
    input_db.reset();
}
//In the destructor, samrai reader and tioga's destructor is called
Orchestra::~Orchestra()
{
    delete[] samrai;
}
//Something Constructor wouldn't do should be put here.
void Orchestra::performGeometricRefineOverset()
{
    int dim = samrai->geo_model->d_dim.getValue();
    std::string local_filename = "aaa_rd_CartGrids" + std::to_string(cur_proc) + ".dat";
    std::string global_filename = "aaa_rd_CartBlocks" + std::to_string(cur_proc) + ".dat";

    //Looks like orchestra need a recorder for time
    samrai->make_geometric_refine();
    MPI_Barrier(MPI_COMM_WORLD);
    //We need to specify which are to be destroyed by SAMRAI and which are released right here.
    int nf;
    int qstride;
    double *qnodein;
    int *idata;
    double *rdata;
    int ngridsin;
    int qnodesize;
    int *global_id;
    int **ibl;
    double **q;
    int nblocksin;
    samrai->pack_patch_for_overset(&nf, &qstride, &qnodein, &idata, &rdata, &ngridsin, &qnodesize, &ibl, &q, &global_id, &nblocksin);
    MPI_Barrier(MPI_COMM_WORLD);

    //NOw let's connect tioga

    tg->register_amr_global_data(nf, qstride, qnodein, idata, rdata, ngridsin, qnodesize);
    tg->set_amr_patch_count(nblocksin);
    for (int i = 0; i < nblocksin; i++)
    {
        tg->register_amr_local_data(i, global_id[i], ibl[i], q[i]);
        std::cout<<"Cur proc is "<<cur_proc<<" this block id is "<<global_id[i]<<'\n';
    }
    tg->performConnectivity();
    MPI_Barrier(MPI_COMM_WORLD);
    tg->performConnectivityAMR();
    MPI_Barrier(MPI_COMM_WORLD);
}
void Orchestra::init()
{
    rd->readAll();
    rd->make_movement(1,-0.7,0.15,0);


    solver->init(rd);
    // 
    // for(int i = 0;i<rd->d_nmesh;i++)
    // {
    //     auto& curBi = rd->bis[i];
    //     curBi.ibl = solver->d_flow_strategy->getIblankCell()[i];
    // }
    
    //std::cout<<"End of readaLL\n";
    MPI_Barrier(MPI_COMM_WORLD);
    tg->registerFromFeeder(rd);
    MPI_Barrier(MPI_COMM_WORLD);
    
    
    tg->profile();
    MPI_Barrier(MPI_COMM_WORLD);
    tg->reResolution();
    //Set iblank_cell from solver to tioga
    for(int i = 0;i<rd->d_nmesh;i++)
    {
        tg->set_cell_iblank(i+1,solver->d_flow_strategy->getIblankCell()[i]);
    }
    tg->performConnectivity();
    MPI_Barrier(MPI_COMM_WORLD);
    //At here,we perform a simple assignment of each value in flow strategy to be its mesh ind
    
    double** U = solver->d_flow_strategy->getU();
    int nequ = solver->d_flow_strategy->getNEquation();
    for(int i = 0;i<rd->d_nmesh;i++)
    {
        int value = i+1;
        for(int k = 0;k<solver->d_flow_strategy->getNEquation()*solver->d_topology_holder->nCells(i);k++)
        {
            U[i][k] = value;
        }
    }

    double** nodeUs = new double*[rd->d_nmesh];
    for(int i = 0;i<rd->d_nmesh;i++)
    {
        int value = i+1;
        nodeUs[i] = new double[nequ*solver->d_topology_holder->nPoints(i)];
        for(int k = 0;k<nequ*solver->d_topology_holder->nPoints(i);k++)
        {
            nodeUs[i][k] = value;
        }
    }
    for(int i = 0;i<rd->d_nmesh;i++)
    {
        tg->mblocks[i]->writeGridFile2("beforeNode",nodeUs[i],nequ);
    }
    //The dataUpdate between nodeUs is done.
    for(int i = 0;i<rd->d_nmesh;i++)
    {
        tg->registerSolution(i+1,nodeUs[i]);
    }

    tg->dataUpdate(nequ,0,0);
    for(int i = 0;i<rd->d_nmesh;i++)
    {
        tg->mblocks[i]->writeGridFile2("afterrNode",nodeUs[i],nequ);
    }
    for(int i = 0;i<rd->d_nmesh;i++)
    {
        int cellId = 0;
        for(int k = 0;k<tg->mblocks[i]->ntypes;k++)
        {
            int nvert = tg->mblocks[i]->nv[k];
            for(int typeC = 0;typeC < tg->mblocks[i]->nc[k];typeC++)
            {
                for(int l = 0;l< nequ; l++)
                {
                    U[i][cellId*nequ+l] = 0;
                    for(int j = 0;j<nvert;j++)
                    {
                        int nodeId = tg->mblocks[i]->vconn[k][nvert*typeC+j];
                        
                        U[i][cellId*nequ+l] += nodeUs[i][nodeId*nequ+l];
                    }
                    U[i][cellId*nequ+l] /= nvert;

                }

                ++cellId;
            }
        }
        // for(int k = 0;k<solver->d_topology_holder->nPoints(i);k++)
        // {
        //     for(int j = 0;j<nequ;j++)
        //     {
        //         int offset = k*nequ+j;
            
        //         if(fabsf64(nodeUs[i][offset] - value) >0.5)
        //         {
        //             std::cout<<nodeUs[i][offset]<<" should be "<<value<<'\n';
        //             fout<< solver->d_topology_holder->blk2D[i].d_localPoints[k][0]
        //             <<" "<<solver->d_topology_holder->blk2D[i].d_localPoints[k][1]<<"\n";
        //         // std::cin.get();
        //         }
        //     }
        //     // nodeUs[i][k] = value;
        // }
        // fout.close();
    }
    //Now it is time to take it back to cellUs
    // for(int i = 0;i<rd->d_nmesh;i++)
    // {
    //     int value = i+1;
    //     auto& curBlk = solver->d_topology_holder->blk2D[i];
        
    //     for(int k = 0;k<solver->d_topology_holder->nCells(i);k++)
    //     {
    //         // std::cout<<k<<" Checking\n";
    //         // std::cout<<"solver->d_topology_holder->nCells(k)"<<solver->d_topology_holder->nCells(k)<<'\n';
    //         auto& curCell = curBlk.d_localCells[k];
            
    //         for(int l = 0;l< nequ;l++)
    //         {    
    //             U[i][k*nequ+l] = 0;
                
    //             for(int j = 0;j<curCell.size();j++)
    //             {
    //                 int nodeId = curCell.pointInd(j);
                    
    //                 if(fabsf64(nodeUs[i][nodeId*nequ+l] - value)>0.5)
    //                 {
    //                     std::cout<<"Updating cell should have some change\n";
    //                 }
    //                 U[i][k*nequ+l] += nodeUs[i][nodeId*nequ+l];
    //             }
    //             U[i][k*nequ+l] = U[i][k*nequ+l]/curCell.size();
    //         }
    //     }
    // }

    tg->mblocks[0]->writeCellFile2("finalTioga0",nequ,U[0]);
    tg->mblocks[1]->writeCellFile2("finalTioga1",nequ,U[1]);

    solver->d_flow_strategy->writeNOdeData("finalnode0",cur_proc,0,nodeUs);
    solver->d_flow_strategy->writeNOdeData("finalnode1",cur_proc,1,nodeUs);

    solver->d_flow_strategy->writeCellData("final0",cur_proc,0,U);
    solver->d_flow_strategy->writeCellData("final1",cur_proc,1,U);
    

    // samrai->init(rd);
    // MPI_Barrier(MPI_COMM_WORLD);
    

    return;
   
    // solver->Loop();

    performGeometricRefineOverset();
    // samrai->dump_data(10);
    // dumpCartesian("I",0.0);
    return;

    //Take all the CartGrid out of SAMRAI into TIOGA

    //extractGrid(&nf,&qstride,qnodein,idata,rdata,&ngridsin,&qnodesize);
}
//Time advancing and Regridding if necessary
void Orchestra::singleStep(int &cur_step, double &cur_time)
{
}

void Orchestra::run(int step, double time)
{
    //Do a while loop here
}

//Something Destructor wouldn't do should be put here.
void Orchestra::finalize()
{
    std::cout<<"Calling finalizing Orchestra\n";
    delete[] samrai;
}

void Orchestra::dumpAll(double current_time)
{
    dumpCartesian("Cart",current_time);
    dumpUnstruct("Unstruct",current_time);
}
//This should dump the data in TIOGA(iblank flag and all its variables) at time current
void Orchestra::dumpCartesian(const std::string& filename,double current_time)
{
    samrai->tec_writer->writeAll();
}
//This should dump the data in SAMRAI(iblank flag and all its variables) at time current
void Orchestra::dumpUnstruct(const std::string& filename,double current_time)
{

}