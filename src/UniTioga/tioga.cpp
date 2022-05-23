#include "tioga.h"
#include <iostream>
using namespace TIOGA;
tioga::~tioga()
{
}
void tioga::setCommunicator(MPI_Comm communicator, int id_proc, int nprocs)
{
    scomm = communicator;
    myid = id_proc;
    numprocs = nprocs;
    sendCount = new int[numprocs];
    recvCount = new int[numprocs];
    //
    // only one mesh block per process for now
    // this can be changed at a later date
    // but will be a fairly invasive change
    //
    // nblocks=0;
    // mb=new MeshBlock[1];
    //
    // instantiate the parallel communication class
    //
    pc = new parallelComm[1];
    pc->myid = myid;
    pc->scomm = scomm;
    pc->numprocs = numprocs;

    // instantiate the parallel communication class
    //
    pc_cart = new parallelComm[1];
    pc_cart->myid = myid;
    pc_cart->scomm = scomm;
    pc_cart->numprocs = numprocs;
}

void tioga::registerGridData(int btag, int nnodes, double *xyz, int *ibl, int nwbc, int nobc,
                             int *wbcnode, int *obcnode, int ntypes, int *nv, int *nc, int **vconn,
                             uint64_t *cell_gid, uint64_t *node_gid)
{
    int iblk;

    auto idxit = tag_iblk_map.find(btag);
    if (idxit == tag_iblk_map.end())
    {
        mtags.push_back(btag);
        mytag.push_back(btag);
        mblocks.push_back(std::unique_ptr<MeshBlock>(new MeshBlock));
        nblocks = mblocks.size();
        iblk = nblocks - 1;
        tag_iblk_map[btag] = iblk;
    }
    else
    {
        iblk = idxit->second;
    }

    auto &mb = mblocks[iblk];
    mb->setData(d_dim, btag, nnodes, xyz, ibl, nwbc, nobc, wbcnode, obcnode, ntypes, nv,
                nc, vconn, cell_gid, node_gid);
    mb->myid = myid;
}

void tioga::registerGridData(int btag, BlockInfo *bi)
{
    int nnodes = bi->nnodes;
    double *xyz = bi->rxyz;
    int *ibl = bi->ibl;
    int nwbc = bi->nwbc;
    int nobc = bi->nobc;
    int *wbcnode = bi->wbc;
    int *obcnode = bi->obc;
    int ntypes = bi->ntypes;
    int *nv = bi->eachNodeCount;
    int *nc = bi->eachCellCount;
    int **vconn = bi->vconn;
    registerGridData(btag, nnodes, xyz, ibl, nwbc, nobc, wbcnode, obcnode, ntypes, nv, nc, vconn);
}

void tioga::registerGlobalBoundary(int btag, BlockInfo *bi)
{
    auto idxit = tag_iblk_map.find(btag);
    if (idxit == tag_iblk_map.end())
    {
        std::cout << "Something wrong with setting global BC\n";
        return;
    }

    auto &mb = mblocks[idxit->second];
    mb->set_global_bc(bi->gNwbc, bi->gNobc, bi->gWxyz, bi->gOxyz);
}
void tioga::set_distanced(bool distancing)
{
    d_uses_distancing = distancing;
}
void tioga::registerFromFeeder(TiogaFeeder *fder)
{
    int mtag = 1;
    for (auto iter = fder->bis.begin(); iter != fder->bis.end(); iter++)
    {
        registerGridData(mtag, &(*iter));
        mtag++;
    }
    if (d_uses_distancing)
    {
        mtag = 1;
        for (auto iter = fder->bis.begin(); iter != fder->bis.end(); iter++)
        {
            registerGlobalBoundary(mtag, &(*iter));
            mtag++;
        }
    }
}

void tioga::register_amr_global_data(int nf, int qstride, double *qnodein, int *idata, double *rdata, int ngridsin, int qnodesize)
{
    if (cg)
    {
        TIOGA_FREE(cg);
    }
    cg = new CartGrid[1];
    cg->myid = myid;
    std::cout<<ngridsin<<"registering  global data\n";
    cg->registerData(d_dim,nf, qstride, qnodein, idata, rdata, ngridsin, qnodesize);
}

void tioga::set_amr_patch_count(int npatchesin)
{
    ncart = npatchesin;
    if(cb) TIOGA_FREE(cb);
    cb = new CartBlock[ncart];
}

void tioga::register_amr_local_data(int ipatch, int global_id, int *iblank, double *q)
{
    cb[ipatch].registerData(d_dim,ipatch, global_id, iblank, q);
}

void tioga::profile()
{
    for (int ib = 0; ib < nblocks; ib++)
    {
        auto &mb = mblocks[ib];
        mb->mexclude = mexclude;
        mb->nfringe = nfringe;
        mb->preprocess();
        mblocks[ib]->writeGridFile2("nodeRes", mblocks[ib]->nodeRes);
            mblocks[ib]->writeCellFile2("cellRes", mblocks[ib]->cellRes);
        //mb->writeMandatoryReceptor("AfterProfile");
    }
}
#include <iostream>
void tioga::performConnectivity()
{
    getHoleMap();
    outputHoleMap("holemaps");
    std::cout << "End of getting Holemap on proc " << myid << "\n";
    for (int ib = 0; ib < nblocks; ib++)
    {
        auto &mb = mblocks[ib];
        // mb->writeGridFile("beforeExchange");
    }
    exchangeBoxes();

    std::cout << "obblist size on proc " << myid << "is " << obblist.size() << '\n';

    if (true)
    {
        for (int ib = 0; ib < nblocks; ib++)
        {
            // mblocks[ib]->writeGridFile("grd");
            // mblocks[ib]->writeGridFile2("nodeRes", mblocks[ib]->nodeRes);
            // mblocks[ib]->writeCellFile2("cellRes", mblocks[ib]->cellRes);
            // mblocks[ib]->writeGridHDF("OneAndOnlyGrid");
            mblocks[ib]->writeOBB("boundingbox");
            for (int i = 0; i < obblist.size(); i++)
            {
                mblocks[0]->writeOBB2("obbLists" + std::to_string(myid) + " " + std::to_string(i) + " ", &obblist[i]);
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    exchangeSearchData();
    // outPutQuery("queries");
    //The code here can't be replaced with auto& mb = mblocks[ib];
    //Why?
    for (int ib = 0; ib < nblocks; ib++)
    {
        mblocks[ib]->resetInterpData();
        mblocks[ib]->search();
        //auto &mb = mblocks[ib];
        //mb->ihigh = 0;
        //mb->resetInterpData();
        //mb->writeMandatoryReceptor("AfterExchangeSearch");
        //mb->search();
        //mb->writeCellFile("LordShowMeTheWay");
    }

    std::cout << "Entering exchangeDonors\n";
    exchangeDonors(false);

    reduce_fringes();
    for (int ib = 0; ib < nblocks; ib++)
    {
        //auto &mb = mblocks[ib];
        if (ihighGlobal)
        {
            std::cout << "Calling getCellIblanks2" << '\n';
            mblocks[ib]->getCellIblanks2();
            //mb->getCellIblanks2();
        }
        else
        {
            std::cout << "Calling getCellIblanks" << '\n';
            mblocks[ib]->getCellIblanks();
            //mb->getCellIblanks();
        }
        mblocks[ib]->writeCellFile2("iblank_cell", mblocks[ib]->iblank_cell);
    }
}

void tioga::performConnectivityAMR()
{

    int iamr;

    iamr = (ncart > 0) ? 1 : 0;
    MPI_Allreduce(&iamr, &iamrGlobal, 1, MPI_INT, MPI_MAX, scomm);
    cg->setDim(d_dim);
    cg->preprocess();
    
    for (int i = 0; i < ncart; i++)
    {
        cb[i].setDim(d_dim);
        cb[i].preprocess(cg);
    }
    
    std::cout<<"Step before if\n";
    
    if (nblocks > 0)
    {
        for (int ib = 0; ib < nblocks; ib++)
        {
            std::cout<<"Stop before getCartReceptors\n";

            auto &mb = mblocks[ib];
            mb->getCartReceptors(cg, pc_cart); //From cartBlock
            mb->ihigh = ihigh;                 //TODO::Where is the ihigh of tioga being set?
            std::cout<<"Stop before mb->search\n";
            mb->search();
            std::cout<<"Stop before getUnresolvedMandatory receptor\n";
            mb->getUnresolvedMandatoryReceptors(); //At meshBlock
            std::cout<<"Stop before cg->search\n";
            cg->search(mb->rxyzCart, mb->donorIdCart, mb->ntotalPointsCart);
        }
    }
    

    //cb[0].writeCellFile("cartBlock");
    // outPutQuery("amrQueries");
    // outPutCartQuery("amrCartQueries");
    std::cout<<"Queries being output\n";
    

    // return;
    exchangeAMRDonors();
    //MPI_Barrier(scomm);
    std::cout<<"End of performAMRCOnnectivity\n";
    for (int ib = 0; ib < nblocks; ib++)
    {
        auto &mb = mblocks[ib];
        mb->getCellIblanks();
        mb->writeGridFile2("IBLANKS_NODE",mb->iblank);
        mb->writeCellFile2("IBLANKS_mb",mb->iblank_cell);
    }
    for(int ib=0;ib<ncart;ib++)
    {
        cb[ib].writeCellFile("IBLANKS_cb");
    }
    
}

//void tioga::test_from_search(double* xsearches, int nsearches)
//{
//    auto& mb = mblocks[0];
//    mb->isolated_search(nsearches, xsearches);
//}

void tioga::reduce_fringes()
{
    int nsend, nrecv;
    int *sndMap;
    int *rcvMap;
    PACKET *sndPack, *rcvPack;
    pc->getMap(&nsend, &nrecv, &sndMap, &rcvMap);
    for (int ib = 0; ib < nblocks; ib++)
    {
        auto &mb = mblocks[ib];
        mb->reduce_fringes();
    }
    //Repeat what the process in exchangeDonors would do

    //Find cancellation data.
    std::vector<int> nrecords(nblocks, 0);
    int **donorRecords = new int *[nblocks];
    for (int i = 0; i < nblocks; i++)
    {
        donorRecords[i] = NULL;
        nrecords[i] = 0;
        mblocks[i]->getCancellationReduce(&(nrecords[i]), &(donorRecords[i]));
    }
    sndPack = new PACKET[nsend];
    rcvPack = new PACKET[nrecv];

    pc->initPackets(sndPack, rcvPack);

    std::vector<int> nintsSend(nsend, 0);
    std::vector<int> ixOffset(nsend, 0);
    std::fill(nintsSend.begin(), nintsSend.end(), 0);
    std::fill(ixOffset.begin(), ixOffset.end(), 0);

    for (int n = 0; n < nblocks; n++)
    {
        for (int i = 0; i < nrecords[n]; i++)
        {
            int k = donorRecords[n][3 * i];
            sndPack[k].nints += 2;
        }
    }

    for (int k = 0; k < nsend; k++)
    {
        sndPack[k].intData = new int[sndPack[k].nints];
    }

    for (int n = 0; n < nblocks; n++)
    {
        for (int i = 0; i < nrecords[n]; i++)
        {
            int k = donorRecords[n][3 * i];
            sndPack[k].intData[ixOffset[k]++] = donorRecords[n][3 * i + 1];
            sndPack[k].intData[ixOffset[k]++] = donorRecords[n][3 * i + 2];
        }
    }

    pc->sendRecvPackets(sndPack, rcvPack);
    //Cancel donors
    for (int k = 0; k < nrecv; k++)
    {
        int m = 0;
        for (int j = 0; j < rcvPack[k].nints / 2; j++)
        {
            int recid = rcvPack[k].intData[m++];
            int ib = tag_iblk_map[rcvPack[k].intData[m++]];
            mblocks[ib]->cancelDonor(recid);
        }
    }
    for (int ib = 0; ib < nblocks; ib++)
    {
        auto &mb = mblocks[ib];
        mb->resetCoincident();
    }

    pc->clearPackets(sndPack, rcvPack);
    //Find remaining donors

    for (int i = 0; i < nblocks; i++)
    {
        if (donorRecords[i])
        {
            TIOGA_FREE(donorRecords[i]);
            donorRecords[i] = NULL;
        }
        nrecords[i] = 0;
        mblocks[i]->getInterpData(&(nrecords[i]), &(donorRecords[i]));
    }
    std::fill(nintsSend.begin(), nintsSend.end(), 0);
    std::fill(ixOffset.begin(), ixOffset.end(), 0);
    for (int n = 0; n < nblocks; n++)
    {
        for (int i = 0; i < nrecords[n]; i++)
        {
            int k = donorRecords[n][3 * i];
            sndPack[k].nints += 2;
        }
    }
    for (int k = 0; k < nsend; k++)
        sndPack[k].intData = new int[sndPack[k].nints];
    for (int n = 0; n < nblocks; n++)
    {
        for (int i = 0; i < nrecords[n]; i++)
        {
            int k = donorRecords[n][3 * i];
            sndPack[k].intData[ixOffset[k]++] = donorRecords[n][3 * i + 1];
            sndPack[k].intData[ixOffset[k]++] = donorRecords[n][3 * i + 2];
        }
    }
    pc->sendRecvPackets(sndPack, rcvPack);
    //Find all valid iblanks

    for (int ib = 0; ib < nblocks; ib++)
    {
        mblocks[ib]->clearIblanks();
    }

    for (int k = 0; k < nrecv; k++)
    {
        int m = 0;
        for (int j = 0; j < rcvPack[k].nints / 2; j++)
        {
            int recid = rcvPack[k].intData[m++];
            int ib = rcvPack[k].intData[m++];
            //std::cout << mblocks[ib]->iblank << '\n';
            //std::cout << "pointid " << pointid << " , " << "ib " << ib << '\n';
            mblocks[ib]->setIblanks(recid);
        }
    }
    for (int ib = 0; ib < nblocks; ib++)
    {
        auto &mb = mblocks[ib];
        mb->writeGridFile("FindReduced");
    }
    pc->clearPackets(sndPack, rcvPack);
    TIOGA_FREE(sndPack);
    TIOGA_FREE(rcvPack);
    if (donorRecords)
    {
        for (int i = 0; i < nblocks; i++)
        {
            if (donorRecords[i])
                TIOGA_FREE(donorRecords[i]);
        }
        TIOGA_FREE(donorRecords);
    }
}