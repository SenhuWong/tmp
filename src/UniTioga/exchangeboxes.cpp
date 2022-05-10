//
// This file is part of the Tioga software library
//
// Tioga  is a tool for overset grid assembly on parallel distributed systems
// Copyright (C) 2015 Jay Sitaraman
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

#include <vector>
#include <map>
#include <numeric>
#include <algorithm>
#include <utility>
#include<iostream>
#include <cassert>
#include "codetypes.h"
#include "tioga.h"
using namespace TIOGA;
int obbIntersectCheck(int dim, double vA[3][3], double xA[3], double dxA[3],
    double vB[3][3], double xB[3], double dxB[3]);
void tioga::exchangeBoxes(void)
{
    int* sndMap;
    int* rcvMap;
    int nsend;
    int nrecv;
    PACKET* sndPack, * rcvPack;

    std::vector<int> nbPerProc(numprocs); // Number of chunks per processor
    std::vector<int> obSizePerProc(numprocs); // Number of real data (OBBs) per processor

    MPI_Allgather(&nblocks, 1, MPI_INT, nbPerProc.data(), 1, MPI_INT, scomm);

    // Total number mesh chunks across all procs
    int ntotalblks = std::accumulate(nbPerProc.begin(), nbPerProc.end(), 0);

    std::vector<int> alltags(ntotalblks); // Mesh tags for all blocks across all procs
    std::vector<int> displs(numprocs + 1);  // Offsets for tags per proc
    std::vector<bool> sendFlag(numprocs, false); // Flag indicating send/recv from this proc

    displs[0] = 0;
    //%TODO 16 is a good number for OBB data of rank 3.
    obSizePerProc[0] = nbPerProc[0] * 16;
    for (int i = 1; i <= numprocs; i++) 
    {
        displs[i] = displs[i - 1] + nbPerProc[i - 1];
        if (i < numprocs) obSizePerProc[i] = nbPerProc[i] * 16;
    }

    MPI_Allgatherv(mytag.data(), nblocks, MPI_INT, alltags.data(),
        nbPerProc.data(), displs.data(), MPI_INT, scomm);

    int maxtag = -1;
    //for (auto itag: alltags)
    for (int i = 0; i < ntotalblks; i++) 
    {
        int itag = abs(alltags[i]);
        if (maxtag < itag) maxtag = itag;
    }
    int mxtgsqr = maxtag * maxtag;

    displs[0] = 0;
    for (int i = 1; i <= numprocs; i++) {
        displs[i] = displs[i - 1] + nbPerProc[i - 1] * 16;
    }

    std::vector<double> myOBBdata(nblocks * 16);
    std::vector<double> allOBBdata(ntotalblks * 16);

    int m = 0;
    for (int ib = 0; ib < nblocks; ib++) 
    {
        myOBBdata[m++] = (double)mytag[ib];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                myOBBdata[m++] = mblocks[ib]->obb->vec[i][j];
        for (int i = 0; i < 3; i++)
            myOBBdata[m++] = mblocks[ib]->obb->xc[i];
        for (int i = 0; i < 3; i++)
            myOBBdata[m++] = mblocks[ib]->obb->dxc[i];
    }

    MPI_Allgatherv(myOBBdata.data(), nblocks * 16, MPI_DOUBLE, allOBBdata.data(),
        obSizePerProc.data(), displs.data(), MPI_DOUBLE, scomm);

    // Determine total number of OBBs received 
    int nobb = ntotalblks;

    // Store all received OBBs in a temporary list
    std::vector<OBB> obbRecv(nobb);
    std::vector<int> obbID(nobb); // Mesh tag corresponding to the OBB
    std::vector<int> obbProc(nobb); // Proc ID corresponding to OBB
    m = 0;
    for (int k = 0, ix = 0; k < numprocs; k++) {
        for (int n = 0; n < nbPerProc[k]; n++)
        {
            obbProc[ix] = k;//Proc id of the OBB
            obbID[ix] = (int)(allOBBdata[m++] + 0.5);//Meshtag of the OBB
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    obbRecv[ix].vec[i][j] = allOBBdata[m++];

            for (int i = 0; i < 3; i++)
                obbRecv[ix].xc[i] = allOBBdata[m++];

            for (int i = 0; i < 3; i++)
                obbRecv[ix].dxc[i] = allOBBdata[m++];
            ix++;
        }
    }

    // Mapping of (local block_id, remote OBB block_id) for every intersected pair
    std::vector<std::pair<int, int>> intersectIDs;
    // Counter tracking number of intersected blocks per process
    std::vector<int> obPerProc(numprocs, 0);

    // Reset sendflags 
    std::fill(sendFlag.begin(), sendFlag.end(), false);
    //
    // Check for intersection of OBBs
    //
    nsend = nrecv = numprocs;
    for (int ob = 0; ob < nobb; ob++) 
    {
        for (int ib = 0; ib < nblocks; ib++) 
        {
            auto& mb = mblocks[ib];
            int meshtag = mb->getMeshTag();
            if (abs(obbID[ob]) == meshtag) continue;

            if (obbIntersectCheck(d_dim,
                mb->obb->vec, mb->obb->xc, mb->obb->dxc,
                obbRecv[ob].vec, obbRecv[ob].xc, obbRecv[ob].dxc)==1 or
                obbIntersectCheck(d_dim,
                    obbRecv[ob].vec, obbRecv[ob].xc, obbRecv[ob].dxc,
                    mb->obb->vec, mb->obb->xc, mb->obb->dxc)==1) 
            {
                int overlap_present = 1;
                //if (mb->meshtag == 3)//BackGround's pairs
                //{
                //    mb->writeOBBPair(std::to_string(ib) + " " + std::to_string(ob) + "OBBPAIR", &(obbRecv[ob]));
                //}
                if (obbID[ob] < 0 || mytag[ib] < 0) //obbID or mytag is negative means p4est routine.
                {
                    mb->check_intersect_p4est(&obbProc[ob], &overlap_present);
                }
                // If there is an intersection, store the index pair, increment number
                // of intersections for this processor, and activate send flag
                if (overlap_present) 
                {
                    intersectIDs.push_back(std::make_pair(ib, ob));//index of nblocks and nobbs
                    obPerProc[obbProc[ob]]++;//The remote proc intersection
                    sendFlag[obbProc[ob]] = true;//The remote proc need communication
                }
            }
        }
    }
    int new_send = std::count(sendFlag.begin(), sendFlag.end(), true);
    assert(new_send <= nsend);
    // Populate send and recv maps
    std::map<int, int> invMap;
    nsend = nrecv = new_send;

    // allocate sndMap and rcvMap
    sndMap = new int[nsend];
    rcvMap = new int[nrecv];
    sndPack = new PACKET[nsend];
    rcvPack = new PACKET[nrecv];
    
    for (int p = 0, ip = 0; p < numprocs; p++)
    {
        if (sendFlag[p]) 
        {
            sndMap[ip] = p;
            rcvMap[ip] = p;
            invMap[p] = ip;
            ip++;
        }
    }
    // clear packets before nsend and nrecv are modified in pc->setMap
    pc->setMap(nsend, nrecv, sndMap, rcvMap);
    pc->initPackets(sndPack, rcvPack);


    
    intBoxMap.clear();
    ibsPerProc.clear();
    ibProcMap.clear();
    obblist.clear();
    obblist.resize(intersectIDs.size());
    ibsPerProc.resize(nsend);
    ibProcMap.resize(nsend);

    // Determine packet sizes and reallocate arrays
    for (int k = 0; k < nsend; k++) 
    {
        sndPack[k].nints = 3 * obPerProc[sndMap[k]];
        sndPack[k].nreals = 6 * obPerProc[sndMap[k]];
        sndPack[k].intData = new int[sndPack[k].nints];
        sndPack[k].realData = new double[sndPack[k].nreals];
        ibsPerProc[k] = obPerProc[sndMap[k]];
        ibProcMap[k].resize(ibsPerProc[k]);
    }

    // Array tracking indices for populating reduced OBBs
    std::vector<int> idxOffset(nsend, 0);
    for (size_t i = 0; i < intersectIDs.size(); i++) 
    {
        auto ids = intersectIDs[i];
        int ib = ids.first;           // Block ID of the local mesh block
        int ob = ids.second;          // Index of the intersected block in OBB list 
        int k = invMap[obbProc[ob]];  // Index in sndMap for this proc ID
        auto& mb = mblocks[ib];       // Mesh block data object
        int ip = obbProc[ob];

        int ioff = idxOffset[k];      // Index to fill in sndPack
        int roff = idxOffset[k] * 6;

        int key_recv = mxtgsqr * ip + maxtag * (mtags[ib] - 1) + obbID[ob] - 1;
        int key_send = mxtgsqr * myid + maxtag * (obbID[ob] - 1) + (mtags[ib] - 1);
        intBoxMap[key_recv] = i;
        ibProcMap[k][ioff] = i;
        obblist[i].comm_idx = k;
        obblist[i].iblk_local = ib;
        // obblist[i].iblk_remote = obbID[ob];
        obblist[i].send_tag = key_send;
        obblist[i].recv_tag = key_recv;

        sndPack[k].intData[3 * ioff] = key_send; // mb->getMeshTag();
        sndPack[k].intData[3 * ioff + 1] = ib;
        sndPack[k].intData[3 * ioff + 2] = mb->getMeshTag();
        //mb->getReducedOBB2(&obbRecv[ob], &(sndPack[k].realData[roff]));
        mb->getReducedOBB(&obbRecv[ob], &(sndPack[k].realData[roff]));

        for (int ii = 0; ii < 3; ii++)
            for (int j = 0; j < 3; j++)
                obblist[i].vec[ii][j] = obbRecv[ob].vec[ii][j];

        // Increment index offset for next fill
        idxOffset[k]++;
        // std::cout << "# " << myid << " " << obbProc[ob] << " "
        //           << mtags[ib] << " " << obbID[ob] << " "
        //           << std::setw(3) << key_recv << " "
        //           << std::setw(3) << key_send << std::endl;
    }
    pc->sendRecvPackets(sndPack, rcvPack);
    for (int i = 0; i < nsend; i++)
    {
        std::cout << "sndpack " << i <<" on proc "<<myid<< " is:";
        for (int j = 0; j < sndPack[i].nreals; j++)
        {
            std::cout << sndPack[i].realData[j] << '\t';
        }
        std::cout << '\n';
    }
    //int count = 0;
    for (int k = 0; k < nrecv; k++) 
    {
        int m = 0;
        for (int n = 0; n < rcvPack[k].nints; n += 3) 
        {
            int key = rcvPack[k].intData[n];
            int ii = intBoxMap[key];
            obblist[ii].iblk_remote = rcvPack[k].intData[n + 1];
            obblist[ii].tag_remote = rcvPack[k].intData[n + 2];
            

            for (int i = 0; i < d_dim; i++)
            {
                obblist[ii].xc[i] = rcvPack[k].realData[m++];
                std::cout << "RAW XC " << obblist[ii].xc[i] << '\n';
            }
            for (int i = 0; i < d_dim; i++)
            {
                obblist[ii].dxc[i] = rcvPack[k].realData[m++];
                std::cout <<"RAW DXC " << obblist[ii].dxc[i] << '\n';
            }
            for (int i = 0; i < 3-d_dim; i++)
            {
                m += 2;
            }
            /*if (obblist[ii].tag_remote == 3 and obblist[ii].dxc[0]!=0)
            {
                mblocks[0]->writeOBB2("for3"+std::to_string(count), &obblist[ii]);
                count++;
            }*/
        }
    }
    for (int ii = 0; ii < obblist.size(); ii++)
    {
        int ib = obblist[ii].iblk_local;
        auto& mb = mblocks[ib];
        std::cout << "End of exchangeBoxes::" << ib << " versus " << obblist[ii].iblk_remote << "on proc " << myid << "out of obblist " << obblist.size() << "\n";
        std::cout << "End of exchangeBoxes::" << obblist[ii].dxc[0] << " " << obblist[ii].dxc[1] << '\n';
        std::cout << "End of exchangeBoxes::" << obblist[ii].xc[0] << " " << obblist[ii].xc[1] << '\n';
    }
    ////Do a cleansing by erasing all obblists with xc and dxc equal to 0
    //auto iter = obblist.begin();
    //bool signal = true;
    //while (signal)
    //{
    //    if (iter->dxc[0] == 0)
    //    {
    //        obblist.erase(iter);
    //        iter = obblist.begin();
    //    }
    //    else
    //    {
    //        signal = false;
    //    }
    //}
    //for (int i=0;i<obblist.size();i++)
    //{
    //    std::cout << i << " " << obblist.size()<<'\n';
    //    if ((iter+i)->dxc[0] == 0)
    //    {
    //        obblist.erase(iter+i);
    //        i--;
    //    }
    //}

    pc->clearPackets(sndPack, rcvPack);
    //
    // Free local memory
    //
    TIOGA_FREE(sndMap);
    TIOGA_FREE(rcvMap);
    TIOGA_FREE(sndPack);
    TIOGA_FREE(rcvPack);
}