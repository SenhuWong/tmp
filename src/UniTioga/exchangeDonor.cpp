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
#include<iostream>
#include <algorithm>
#include "codetypes.h"
#include "tioga.h"
using namespace TIOGA;
void tioga::exchangeDonors(bool reduce_fringe)
{
    int nsend, nrecv;
    int* sndMap;
    int* rcvMap;
    PACKET* sndPack, * rcvPack;
    //
    // get the processor map for sending
    // and receiving
    //
    pc->getMap(&nsend, &nrecv, &sndMap, &rcvMap);
    if (nsend == 0) return;
    //
    // create packets to send and receive
    // and initialize them to zero
    //
    sndPack = new PACKET[nsend];
    rcvPack = new PACKET[nrecv];

    //
    pc->initPackets(sndPack, rcvPack);

    std::vector<int> nintsSend(nsend, 0), nrealsSend(nsend, 0);

    // First pass through all mesh blocks and determine data sizes for each
    // intersection pair
    
    for (int ib = 0; ib < nblocks; ib++) 
    {
        auto& mb = mblocks[ib];
        mb->getMBDonorPktSizes(nintsSend, nrealsSend);
    }
    std::cout << nintsSend[0] <<'\t' << nrealsSend[0] << '\n';
    //
    // Allocate sndPack 
    //
    for (int k = 0; k < nsend; k++)
    {
        sndPack[k].nints = nintsSend[k];
        sndPack[k].nreals = nrealsSend[k];
        sndPack[k].intData = new int[sndPack[k].nints];
        sndPack[k].realData = new double[sndPack[k].nreals];
    }
    //
    // Populate send packets with data from each mesh block in this partition
    //
    std::vector<int> ixOffset(nsend, 0), rxOffset(nsend, 0);
    for (int ib = 0; ib < nblocks; ib++) 
    {
        auto& mb = mblocks[ib];
        mb->getMBDonorPackets(ixOffset, rxOffset, sndPack);
    }
    //
    // communicate donors (comm1)
    //
    pc->sendRecvPackets(sndPack, rcvPack);
    // Initialize linked lists and populate donor data from rcvPack
    for (int ib = 0; ib < nblocks; ib++) 
    {
        auto& mb = mblocks[ib];
        mb->initializeDonorList();
    }
    //
    for (int k = 0; k < nrecv; k++) 
    {
        int m = 0;
        int l = 0;
        for (int i = 0; i < rcvPack[k].nints / 4; i++)
        {
            int meshtag = rcvPack[k].intData[m++];// Meshtag on the other side(donor side).
            int pointid = rcvPack[k].intData[m++];// point id at this side(receptor)
            int remoteid = rcvPack[k].intData[m++];// local receptor's index in remote(donor side)'s nsearch
            int ib = rcvPack[k].intData[m++];// Block id on this(receptor) side
            double donorRes = rcvPack[k].realData[l++];
            double receptorRes = rcvPack[k].realData[l++];

            auto& mb = mblocks[ib];
            // mb->set_Iblank(pointid, -1);
            mb->insertAndSort(pointid, k, meshtag, remoteid, donorRes, receptorRes);
        }
    }
    
    //
    // Figure out the state of each point (i.e., if it is a hole, fringe, or a
    // field point)
    //
    std::vector<int> nrecords(nblocks, 0);
    int** donorRecords = new int* [nblocks];
    double** receptorResolution = new double* [nblocks];
    
    for (int ib = 0; ib < nblocks; ib++) {
        auto& mb = mblocks[ib];
        mb->writeGridFile("preProcessDonor");
        mb->processDonors(holeMap, nmesh, &(donorRecords[ib]),
            &(receptorResolution[ib]), &(nrecords[ib]));
        mb->writeGridFile("processDonor");
    }
    std::cout << "Donors processed\n";
    
    //
    // Reset all send/recv data structures
    //
    pc->clearPackets(sndPack, rcvPack);
    for (int i = 0; i < nsend; i++) {
        sndPack[i].nints = 0;
        sndPack[i].nreals = 0;
        nintsSend[i] = 0;
        nrealsSend[i] = 0;
        ixOffset[i] = 0;
        rxOffset[i] = 0;
    }

    for (int n = 0; n < nblocks; n++) {
        for (int i = 0; i < nrecords[n]; i++) {
            int k = donorRecords[n][3 * i];
            sndPack[k].nints += 2;
            sndPack[k].nreals++;
        }
    }

    for (int k = 0; k < nsend; k++)
    {
        sndPack[k].intData = new int[sndPack[k].nints];
        sndPack[k].realData = new double[sndPack[k].nreals];
    }

    for (int n = 0; n < nblocks; n++) {
        for (int i = 0; i < nrecords[n]; i++) {
            int k = donorRecords[n][3 * i];
            sndPack[k].intData[ixOffset[k]++] = donorRecords[n][3 * i + 1];//local receptor's index in remote nsearch
            sndPack[k].intData[ixOffset[k]++] = donorRecords[n][3 * i + 2];//meshtag of remote donor
            sndPack[k].realData[rxOffset[k]++] = receptorResolution[n][i];
        }
    }
    //
    // comm 2
    // notify donors that they are accepted (for now)
    //
    pc->sendRecvPackets(sndPack, rcvPack);

    std::vector<int> ninterp(nblocks, 0);
    for (int k = 0; k < nrecv; k++)
    {
        int m = 0;
        for (int j = 0; j < rcvPack[k].nints / 2; j++)
        {
            m++;   // skip over point id
            int ib = tag_iblk_map[rcvPack[k].intData[m++]];
            ninterp[ib]++;//How many receptor on each meshtag
        }
    }

    for (int i = 0; i < nblocks; i++) {
        mblocks[i]->initializeInterpList(ninterp[i]);
    }
    
    std::fill(ninterp.begin(), ninterp.end(), 0);

    for (int k = 0; k < nrecv; k++)
    {
        int l = 0;
        int m = 0;
        for (int j = 0; j < rcvPack[k].nints / 2; j++)
        {
            int recid = rcvPack[k].intData[m++];//receptor index in local nsearch
            int ib = tag_iblk_map[rcvPack[k].intData[m++]];
            double receptorRes = rcvPack[k].realData[l++];
            mblocks[ib]->findInterpData(&(ninterp[ib]), recid, receptorRes);
        }
    }
    std::cout << "InterpData Found\n";
    for (int ib = 0; ib < nblocks; ib++) {
        mblocks[ib]->set_ninterp(ninterp[ib]);
    }
    // //Output all the receptor
    // {
    //     std::ofstream fout;
    //     for(int ib=0;ib<nblocks;ib++)
    //     {
    //         fout.open("afterInterpReceptor_"+std::to_string(ib+1)+"_"+std::to_string(myid));
    //         auto& mb = mblocks[ib];
    //         for(int i = 0;i<mb->nsearch;i++)
    //         {
    //             int idnsearch = mb->interp2donor[i];
    //             if(idnsearch>=0)
    //             {
    //                 for(int j = 0;j<d_dim;j++)
    //                 {
    //                     fout<< mb->xsearch[d_dim*i+j]<< '\t';
    //                 }
    //                 fout<< '\n';
    //             }
    //         }
    //         fout.close();
    //     }
    

    // }

    pc->clearPackets(sndPack, rcvPack);
    //
    // Find cancellation data (based on donor quality)
    //
    for (int i = 0; i < nblocks; i++) {
        if (donorRecords[i]) {
            TIOGA_FREE(donorRecords[i]);
            donorRecords[i] = NULL;
        }
        nrecords[i] = 0;
        mblocks[i]->getCancellationData(&(nrecords[i]),
            &(donorRecords[i]));
    }
    std::fill(nintsSend.begin(), nintsSend.end(), 0);
    std::fill(ixOffset.begin(), ixOffset.end(), 0);

    for (int n = 0; n < nblocks; n++) {
        for (int i = 0; i < nrecords[n]; i++) {
            int k = donorRecords[n][3 * i];
            sndPack[k].nints += 2;
        }
    }
    for (int k = 0; k < nsend; k++)
        sndPack[k].intData = new int[sndPack[k].nints];

    for (int n = 0; n < nblocks; n++) {
        for (int i = 0; i < nrecords[n]; i++) {
            int k = donorRecords[n][3 * i];
            sndPack[k].intData[ixOffset[k]++] = donorRecords[n][3 * i + 1];
            sndPack[k].intData[ixOffset[k]++] = donorRecords[n][3 * i + 2];
        }
    }
    //
    // communciate cancellation data comm 3
    //
    pc->sendRecvPackets(sndPack, rcvPack);

    for (int k = 0; k < nrecv; k++) {
        int m = 0;
        for (int j = 0; j < rcvPack[k].nints / 2; j++) {
            int recid = rcvPack[k].intData[m++];
            int ib = tag_iblk_map[rcvPack[k].intData[m++]];
            mblocks[ib]->cancelDonor(recid);
        }
    }

    for (int ib = 0; ib < nblocks; ib++) {
        auto& mb = mblocks[ib];
        mb->resetCoincident();
    }

    //
    pc->clearPackets(sndPack, rcvPack);
    //
    // Find final interpolation data
    //
    for (int i = 0; i < nblocks; i++) {
        if (donorRecords[i]) {
            TIOGA_FREE(donorRecords[i]);
            donorRecords[i] = NULL;
        }
        nrecords[i] = 0;
        mblocks[i]->getInterpData(&(nrecords[i]),
            &(donorRecords[i]));
    }
    std::fill(nintsSend.begin(), nintsSend.end(), 0);
    std::fill(ixOffset.begin(), ixOffset.end(), 0);
    for (int n = 0; n < nblocks; n++) {
        for (int i = 0; i < nrecords[n]; i++) {
            int k = donorRecords[n][3 * i];
            sndPack[k].nints += 2;
        }
    }
    for (int k = 0; k < nsend; k++)
        sndPack[k].intData = new int[sndPack[k].nints];
    for (int n = 0; n < nblocks; n++) {
        for (int i = 0; i < nrecords[n]; i++) {
            int k = donorRecords[n][3 * i];
            sndPack[k].intData[ixOffset[k]++] = donorRecords[n][3 * i + 1];//point id in nnode
            sndPack[k].intData[ixOffset[k]++] = donorRecords[n][3 * i + 2];//block id
        }
    }
    
    //
    // comm 4
    // final receptor data to set iblanks
    //     
    pc->sendRecvPackets(sndPack, rcvPack);
    //
    for (int ib = 0; ib < nblocks; ib++)
        mblocks[ib]->clearIblanks();

    for (int k = 0; k < nrecv; k++) {
        int m = 0;
        for (int j = 0; j < rcvPack[k].nints / 2; j++)
        {
            int pointid = rcvPack[k].intData[m++];
            int ib = rcvPack[k].intData[m++];
            mblocks[ib]->setIblanks(pointid);
        }
    }
    for (int ib = 0; ib < nblocks; ib++)
    {
        auto& mb = mblocks[ib];
        mb->writeGridFile("FindInterp");
    }
    pc->clearPackets(sndPack, rcvPack);
    TIOGA_FREE(sndPack);
    TIOGA_FREE(rcvPack);

    if (donorRecords) {
        for (int i = 0; i < nblocks; i++) {
            if (donorRecords[i]) TIOGA_FREE(donorRecords[i]);
        }
        TIOGA_FREE(donorRecords);
    }
    if (receptorResolution) {
        for (int i = 0; i < nblocks; i++) {
            if (receptorResolution[i]) TIOGA_FREE(receptorResolution[i]);
        }
        TIOGA_FREE(receptorResolution);
    }
}

