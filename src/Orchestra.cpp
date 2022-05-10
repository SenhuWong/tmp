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

    if (input_db->keyExists("Reader"))
    {
        std::shared_ptr<SAMRAI::tbox::Database> reader_db(input_db->getDatabase("Reader"));
        int dim = reader_db->getInteger("dim");

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
    std::ofstream fout;
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
    // return;
    MPI_Barrier(MPI_COMM_WORLD);
    tg->performConnectivityAMR();

    // fout.open(local_filename);
    // fout << nblocksin << '\n';
    // for (int i = 0; i < nblocksin; i++)
    // {
    //     fout << i << '\t' << q[i] << '\t' << ibl[i] << '\n';
    // }
    // fout.close();
    // fout.open(global_filename);
    // fout << ngridsin << '\t' << nf << '\t' << qstride << '\t' << qnodesize << '\t' << '\n';
    // int iOffset = 5 + 2 * dim;
    // int rOffset = 2 * dim;
    // for (int i = 0; i < ngridsin; i++)
    // {
    //     fout << i << '\t';
    //     for (int j = 0; j < iOffset; j++)
    //     {
    //         fout << idata[i * iOffset + j] << '\t';
    //     }
    //     fout << '\n';
    // }
    // fout.close();
    MPI_Barrier(MPI_COMM_WORLD);
}
void Orchestra::init()
{
    rd->readAll();
    // rd->make_movement(1,-0.7,0.15,0);
    //std::cout<<"End of readaLL\n";
    MPI_Barrier(MPI_COMM_WORLD);
    tg->registerFromFeeder(rd);
    MPI_Barrier(MPI_COMM_WORLD);
    
    
    tg->profile();
    MPI_Barrier(MPI_COMM_WORLD);
    tg->reResolution();
    
    samrai->init(rd);
    MPI_Barrier(MPI_COMM_WORLD);

    //For now, let's just perform a pureGeometric refine here
    //This need to be implemented in samraiWrapper as something like
    // samrai->make_geometric_refine()
    // samrai->find_patch_for_overset()
    // samrai->extract patch_for_overset()
    // tioga->register_amr_global_data()
    // tioga->register_amr_local_data()
    performGeometricRefineOverset();
    samrai->dump_data(10);
    dumpCartesian("I",0.0);
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