#include "tioga.h"

using namespace TIOGA;

void tioga::exchangeAMRDonors()
{

    int nsend_sav;
    int nrecv_sav;
    int nsend;
    int nrecv;
    int *sndMap;
    int *rcvMap;
    int *sndMapAll;
    int *rcvMapAll;
    //Last time pc is set at getCartReceptor
    //We set it to be the ones with cartGrid's amr node inside OBB;
    //But that is not helpful, so we
    pc_cart->getMap(&nsend_sav, &nrecv_sav, &sndMap, &rcvMap);
    sndMapAll = new int[pc_cart->numprocs];
    rcvMapAll = new int[pc_cart->numprocs];
    nsend = nrecv = pc_cart->numprocs;

    int *imap = new int[pc_cart->numprocs];
    int *icount = new int[pc_cart->numprocs];

    for (int i = 0; i < pc_cart->numprocs; i++)
    {
        sndMapAll[i] = rcvMapAll[i] = imap[i] = i;
        icount[i] = 0;
    }
    pc_cart->setMap(nsend, nrecv, sndMapAll, rcvMapAll);

    PACKET *sndPack = new PACKET[nsend];
    PACKET *rcvPack = new PACKET[nrecv];
    int *obdonors = new int[nsend];
    int *obreceptors = new int[nsend];

    pc_cart->initPackets(sndPack, rcvPack);

    for (int i = 0; i < nsend; i++)
    {
        obdonors[i] = obreceptors[i] = 0;
    }

    if (nblocks > 0)
    {
        for (int ib = 0; ib < nblocks; ib++)
        {
            auto &mb = mblocks[ib];
            for (int i = 0; i < mb->ntotalPointsCart; i++)
            {
                if (mb->donorIdCart[i] != -1) //Finds a valid donor on cb
                {
                    int gid = mb->donorIdCart[i];
                    obdonors[imap[cg->proc_id[gid]]]++;
                    //nobdonors is the number of valid receptor on Carts side
                }
            }

            for (int i = 0; i < mb->nsearch; i++)
            {
                if (mb->donorId[i] != -1)
                {
                    obreceptors[imap[mb->isearch[3 * i]]]++;
                    //nobreceptor is the number of valid donor on meshBlocks
                }
            }
        }
    }
    //Allocate data packets

    for (int i = 0; i < nsend; i++)
    {
        sndPack[i].nints = obdonors[i] * 3 + obreceptors[i] * 5 + 2; //obreceptor counts for donors on mb
        sndPack[i].nreals = obdonors[i] * d_dim + obreceptors[i];
        sndPack[i].intData = new int[sndPack[i].nints];
        sndPack[i].realData = new double[sndPack[i].nreals];
        sndPack[i].intData[0] = obdonors[i]; //number of
        sndPack[i].intData[1] = obreceptors[i];
    }

    //pack the data

    int *intcount = new int[nsend];
    int *realcount = new int[nsend];
    for (int i = 0; i < nsend; i++)
    {
        intcount[i] = 2;
        realcount[i] = 0;
    }
    if (nblocks > 0)
    {
        for (int ib = 0; ib < nblocks; ib++)//First allocate obdonors(intData[0]) for all meshBLocks
        {
            auto &mb = mblocks[ib];
            for (int i = 0; i < mb->ntotalPointsCart; i++)
            {
                if (mb->donorIdCart[i] != -1)
                {
                    int gid = mb->donorIdCart[i];
                    int procid = imap[cg->proc_id[gid]];
                    int localid = cg->local_id[gid];
                    sndPack[procid].intData[intcount[procid]++] = localid;
                    sndPack[procid].intData[intcount[procid]++] = i;
                    sndPack[procid].intData[intcount[procid]++] = ib;
                    for (int j = 0; j < d_dim; j++)
                    {
                        sndPack[procid].realData[realcount[procid]++] = mb->rxyzCart[d_dim * i + j];
                    }
                }
            }
        }
        for(int ib=0;ib<nblocks;ib++)//Then allocate obreceptors for all meshBlocks.
        {
            auto &mb = mblocks[ib];
            for (int i = 0; i < mb->nsearch; i++)
            {
                if (mb->donorId[i] != -1)
                {
                    int procid = imap[mb->isearch[3 * i]];
                    sndPack[procid].intData[intcount[procid]++] = mb->isearch[3 * i + 1];
                    sndPack[procid].intData[intcount[procid]++] = mb->isearch[3 * i + 2];
                    sndPack[procid].intData[intcount[procid]++] = mb->meshtag;
                    sndPack[procid].intData[intcount[procid]++] = i;
                    sndPack[procid].intData[intcount[procid]++] = ib;
                    sndPack[procid].realData[realcount[procid]++] = mb->cellRes[mb->donorId[i]];
                }
            }
        }
    }
    TIOGA_FREE(realcount);

    pc_cart->sendRecvPackets(sndPack, rcvPack);
    
    //decode the data

    int *bcount = new int[ncart];
    for (int i = 0; i < ncart; i++)
    {
        cb[i].clearLists();
        cb[i].initializeLists();
        bcount[i] = 0;
    }

    for (int i = 0; i < nrecv; i++)
    {
        if (rcvPack[i].nreals > 0)
        {
            int m = 2;
            int n = 0;
            int interpCount = rcvPack[i].intData[0]; //interp to meshblock from cb
            int donorCount = rcvPack[i].intData[1];  //donor on meshblock for cb
            icount[i] = interpCount + donorCount;
            

            for (int j = 0; j < interpCount; j++)
            {
                int local_id = rcvPack[i].intData[m++];
                int remote_id = rcvPack[i].intData[m++];
                int remote_blockid = rcvPack[i].intData[m++];
                double xtemp[3];
                for (int k = 0; k < d_dim; k++)
                {
                    xtemp[k] = rcvPack[i].realData[n++];
                }
                cb[local_id].insertInInterpList(i, remote_id, remote_blockid, xtemp);
                bcount[local_id]++;
                std::cout<<"ncart vs local_id:"<<ncart<<" "<<local_id<<'\n';
            }

            for (int j = 0; j < donorCount; j++)
            {

                int local_id = rcvPack[i].intData[m++];
                int index = rcvPack[i].intData[m++];
                int meshtag = rcvPack[i].intData[m++];
                int remote_id = rcvPack[i].intData[m++];
                int remote_blockid = rcvPack[i].intData[m++];

                double cellRes = rcvPack[i].realData[n++];
                cb[local_id].insertInDonorList(i, index, meshtag, remote_id, remote_blockid, cellRes);
                bcount[local_id]++;
                // std::cout<<"ncart vs local_id:"<<ncart<<" "<<local_id<<'\n';
            }
         }
    }
    

    for (int i = 0; i < nsend; i++)
    {
        std::cout << i << " : "
                  << "obdonors " << obdonors[i] << "  obreceptors " << obreceptors[i] << '\n';
    }

    for (int i = 0; i < ncart; i++)
    {
        cb[i].writeCellFile("BeforeProcessDonors");
        cb[i].processDonors(holeMap, nmesh);
        cb[i].writeCellFile("AfterProcessDonors");
    }



    pc_cart->clearPackets2(sndPack, rcvPack);
    for (int i = 0; i < nsend; i++)
    {
        if (icount[i] > 0)
        {
            sndPack[i].nints = 3 * icount[i]; //INterpCount + DOnorCOunt
            sndPack[i].nreals = 0;
            sndPack[i].intData = new int[sndPack[i].nints];
        }
        intcount[i] = 0; //for later use
    }
    for (int i = 0; i < ncart; i++)
        std::cout << bcount[i] << "\n";

    int *cancelledData = NULL;
    for (int i = 0; i < ncart; i++)
    {
        if (cancelledData)
            TIOGA_FREE(cancelledData);
        int ncancel = bcount[i];
        if (ncancel > 0)
        {
            cancelledData = new int[4 * ncancel];
            cb[i].getCancellationData(cancelledData, &ncancel);
            for (int j = 0; j < ncancel; j++)
            {
                int procid = cancelledData[4 * j];
                int ctype = cancelledData[4 * j + 1];
                int remote_id = cancelledData[4 * j + 2];
                int remote_blockid = cancelledData[4 * j + 3];
                sndPack[procid].intData[intcount[procid]++] = ctype;
                sndPack[procid].intData[intcount[procid]++] = remote_id;
                sndPack[procid].intData[intcount[procid]++] = remote_blockid;
            }
        }
    }

    for (int i = 0; i < nsend; i++)
    {
        sndPack[i].nints = intcount[i];
    }
    pc_cart->sendRecvPackets(sndPack, rcvPack);
    //decode now
    for (int i = 0; i < nrecv; i++)
    {
        if (rcvPack[i].nints > 0)
        {
            int m = 0;
            for (int j = 0; j < rcvPack[i].nints / 3; j++)
            {
                int ctype = rcvPack[i].intData[m++];
                int id = rcvPack[i].intData[m++];
                int ib = rcvPack[i].intData[m++];
                if (ctype == 0)
                {
                    //mblocks[ib]->donorIdCart[id] = -1;
                }
                else
                {
                    mblocks[ib]->donorId[id] = -1;
                }
            }
        }
    }
    for (int ib = 0; ib < nblocks; ib++)
    {
        mblocks[ib]->setCartIblanks();
    }
    pc_cart->clearPackets2(sndPack, rcvPack);
    for (int ib = 0; ib < nblocks; ib++)
    {
        mblocks[ib]->findInterpListCart();
    }
    if (cancelledData)
        TIOGA_FREE(cancelledData);
    TIOGA_FREE(bcount);
    TIOGA_FREE(intcount);
    TIOGA_FREE(sndMapAll);
    TIOGA_FREE(rcvMapAll);
    TIOGA_FREE(imap);
    TIOGA_FREE(icount);
    TIOGA_FREE(sndPack);
    TIOGA_FREE(rcvPack);
}