#include"MeshBlock.h"
#include"mpi.h"
#include<iostream>
#include<fstream>
void deallocateLinkList(DONORLIST* temp);
void insertInList(DONORLIST** donorList, DONORLIST* temp1);
void deallocateLinkList2(INTEGERLIST* temp);
int checkHoleMap(int dim, double* x, int* nx, int* sam, double* extents);
void computeNodalWeights(int dim, double xv[8][3], double* xp, double* frac, int nvert);
void MeshBlock::getMBDonorPktSizes(std::vector<int>& nints, std::vector<int>& nreals)
{
    for (int i = 0; i < nsearch; i++) 
    {
        if (donorId[i] > -1) 
        {
            int ii = isearch[3 * i];
            nints[ii] += 4;
            nreals[ii] += 2;
        }
    }
}
void MeshBlock::getMBDonorPackets(std::vector<int>& ixOffset, std::vector<int>& rxOffset, PACKET* sndPack)
{
    for (int i = 0; i < nsearch; i++) 
    {
        if (donorId[i] < 0) continue;

        int k = isearch[3 * i];// proc ID
        int& ix = ixOffset[k];
        int& rx = rxOffset[k];

        sndPack[k].intData[ix++] = meshtag;           // Unique mesh tag
        sndPack[k].intData[ix++] = isearch[3 * i + 1];  // point ID
        sndPack[k].intData[ix++] = i;                 // point ID on donor side(index in nsearch)
        sndPack[k].intData[ix++] = isearch[3 * i + 2];  // receptor block ID
        sndPack[k].realData[rx++] = cellRes[donorId[i]];
        sndPack[k].realData[rx++] = res_search[xtag[i]];
    }
}

void MeshBlock::initializeDonorList(void)
{
    int i;
    if (donorList)
    {
        for (i = 0; i < donorListLength; i++) 
        {
            deallocateLinkList(donorList[i]);
            //printf("\n\t Deallocate (nnodes) i: %d %d ",nnodes,i);
        }
        TIOGA_FREE(donorList);
    }

    donorListLength = nnodes;
    donorList = new DONORLIST * [donorListLength];
    for (i = 0; i < donorListLength; i++)
        donorList[i] = NULL;
}

void MeshBlock::insertAndSort(int pointid, int senderid, int meshtagdonor, int remoteid,
    double donorRes, double receptorRes)
{
    DONORLIST* temp1;
    temp1 = new DONORLIST[1];

    temp1->donorData[0] = senderid;
    temp1->donorData[1] = meshtagdonor;
    temp1->donorData[2] = remoteid;
    temp1->donorRes = donorRes;
    temp1->receptorRes = receptorRes;
    //std::cout << pointid << '\n';
    insertInList(&donorList[pointid], temp1);
}

void MeshBlock::processDonors(HOLEMAP* holemap, int nmesh, int** donorRecords, double** receptorResolution,
    int* nrecords)
{
    
    int nvert;
    DONORLIST* temp;
    int* iflag;
    int meshtagdonor;
    //
    // first mark hole points
    //
    iflag = new int[nmesh];
    //
    
    for (int i = 0; i < nnodes; i++)
    {
        iblank[i] = 1;
        if (donorList[i] == NULL)//It has no donor, which means it is outside of every mesh(with holemap)'s reach, so every one must check its holemap; 
        {
            for (int j = 0; j < nmesh; j++)
                if (j != (meshtag - BASE) and holemap[j].existWall)
                {
                    if (checkHoleMap(d_dim,&x[d_dim * i], holemap[j].nx, holemap[j].sam, holemap[j].extents))
                    {
                        iblank[i] = 0;
                        break;
                    }
                }
        }
        else//Only need to check those who don't have a donor for this node
        {
            temp = donorList[i];
            for (int j = 0; j < nmesh; j++) iflag[j] = 0;
            while (temp != NULL)
            {
                meshtagdonor = temp->donorData[1] - BASE;
                iflag[meshtagdonor] = 1;
                //The nodeRes might be from 2 splitting and therefore vary, but a donor cell will have both's resolution when we make a unique octree in search()
                //, there fore we need to rewrite the nodeRes to be the bigger one of the two.
                nodeRes[i] = TIOGA_MAX(nodeRes[i], temp->receptorRes);
                temp = temp->next;
            }
            for (int j = 0; j < nmesh; j++)
            {
                if (j != (meshtag - BASE) and holemap[j].existWall)
                {
                    if (!iflag[j])
                        if (checkHoleMap(d_dim,&x[d_dim * i], holemap[j].nx, holemap[j].sam, holemap[j].extents))
                        {
                            iblank[i] = 0;
                            break;
                        }
                }
            }
        }
    }
    //not everyone has a wbc
    //std::cout << wbcnode[0] << '\n';
    for (int i = 0; i < nwbc; i++)
    {
        if (iblank[wbcnode[i] - BASE] == 0) 
        {
            std::cout << "--------------------------------------------------------------------\n";
            std::cout << "Alarm from process " << myid << ": wall node is being tagged as a hole " << wbcnode[i] - BASE << " " << donorList[wbcnode[i] - BASE] << "\n";
            
            int ii = wbcnode[i] - BASE;
           /* if (d_dim == 2)
            {
                std::cout << "xloc = (" << x[2 * ii] << "," << x[2 * ii + 1] << ")\n";
            }*/
            /*else if (d_dim == 3)
            {
                std::cout << "xloc = (" << x[3 * ii] << "," << x[3 * ii + 1] << "," << x[3 * ii + 2] << ")\n";
            }*/
            std::cout << "Computations will continue, but may suffer from accuracy problems\n";
            std::cout << "Please recheck positions of your grids\n";
            std::cout << "--------------------------------------------------------------------\n";
        }
    }
    
    //
    // mark mandatory fringes as neighbors (up to nfringe depth)
    // of hole points 
    //
    int* mtag = new int[nnodes];
    int* mtag1 = new int[nnodes];

    for (int i = 0; i < nnodes; i++)
    {
        mtag[i] = mtag1[i] = 0;
        if (iblank[i] == 0) mtag[i] = mtag1[i] = 1;
    }

    for (int iter = 0; iter < nfringe; iter++)
    {
        for (int n = 0; n < ntypes; n++)
        {
            nvert = nv[n];
            for (int i = 0; i < nc[n]; i++)
            {
                for (int m = 0; m < nvert; m++)
                {
                    if (mtag[(vconn[n][nvert * i + m] - BASE)] == 1)
                    {
                        for (int mm = 0; mm < nvert; mm++)
                            if (m != mm && mtag[vconn[n][nvert * i + mm] - BASE] != 1)
                                mtag1[vconn[n][nvert * i + mm] - BASE] = 1;
                    }
                }
            }
        }
        for (int i = 0; i < nnodes; i++) mtag[i] = mtag1[i];
    }
    for (int i = 0; i < nnodes; i++)
        if (mtag1[i] and iblank[i]) nodeRes[i] = BIGVALUE;//Set those with iblank 1 yet inside nfringe to mandatory receptor as they are close to some other's wall
    TIOGA_FREE(mtag);
    TIOGA_FREE(mtag1);
    //
    // now find fringes
    //
    *nrecords = 0;
    for (int i = 0; i < nnodes; i++)
    {
        if (donorList[i] != NULL && iblank[i] != 0)
        {
            //For debug
            //iblank[i] = -1;
            //continue;
            //End
            //Up to now is Okay.
            temp = donorList[i];
            while (temp != NULL)
            {
                if (temp->donorRes < nodeRes[i])
                {
                    iblank[i] = -temp->donorData[1];//If there is a valid donor, iblank should be negative meshtag
                    (*nrecords)++;//Don't forget to look at the sorting of insertInList
                    break;
                }
                temp = temp->next;
            }
        }
    }
   
    //
    // set the records to send back to the donor
    // process
    //
    (*donorRecords) = new int[(*nrecords) * 3];
    (*receptorResolution) = new double[(*nrecords)];

    int m = 0;
    int k = 0;
    for (int i = 0; i < nnodes; i++)
    {
        if (iblank[i] < 0)
        {
            temp = donorList[i];
            while (temp != NULL)
            {
                if (temp->donorRes < nodeRes[i])
                {
                    break;
                }
            }
            (*receptorResolution)[k++] = (resolutionScale > 1.0) ? -nodeRes[i] : nodeRes[i];//The receptorResolution should be negative if resolutionScale is not 1.0
            (*donorRecords)[m++] = temp->donorData[0];
            (*donorRecords)[m++] = temp->donorData[2];
            (*donorRecords)[m++] = temp->donorData[1];
            /*temp1->donorData[0] = senderid;
            temp1->donorData[2] = remoteid;
            temp1->donorData[1] = meshtagdonor;*/
        }
    }
    std::string filename = "mandatoryPostProcessD" + std::to_string(meshtag) + std::to_string(myid) + ".sct";
    std::ofstream fout;
    fout.open(filename);
    if (!fout.is_open())
    {
        std::cout << "Open failure\n";
    }
    else
    {
        for (int i = 0; i < nnodes; i++)
        {
            if (nodeRes[i] == BIGVALUE)
            {
                for (int j = 0; j < d_dim; j++)
                {
                    fout << x[d_dim * i + j] << "\t";
                }
                fout << '\n';
            }
        }
        fout.close();
    }
    
    //
    // release local memory
    //
    TIOGA_FREE(iflag);
}


void MeshBlock::initializeInterpList(int ninterp_input)
{
    if (interpList) {
        //for(i=0;i<ninterp;i++)
        for (int i = 0; i < interpListSize; i++)
        {
            if (interpList[i].inode) TIOGA_FREE(interpList[i].inode);
            if (interpList[i].weights) TIOGA_FREE(interpList[i].weights);
        }
        TIOGA_FREE(interpList);
    }
    ninterp = ninterp_input;
    interpListSize = ninterp_input;
    interpList = new INTERPLIST[interpListSize];
    for (int i = 0; i < interpListSize; i++) 
    {
        interpList[i].inode = NULL;
        interpList[i].weights = NULL;
    }
    if (cancelList) deallocateLinkList2(cancelList);
    cancelList = NULL;
    ncancel = 0;
    if (interp2donor) TIOGA_FREE(interp2donor);
    interp2donor = new int[nsearch];
    for (int i = 0; i < nsearch; i++)
    {
        interp2donor[i] = -1;
    }

}
static int cancel_count = 0;
void MeshBlock::findInterpData(int* recid, int irecord, double receptorRes2)
{
    int procid, pointid, blockid;
    double xv[8][3];
    double xp[3];
    double frac[8];
    int inode[8];
    int acceptFlag;
    double receptorRes;
    int meshtagrecv;

    INTEGERLIST* clist;
    //
    receptorRes = fabs(receptorRes2);
    procid = isearch[3 * irecord];
    pointid = isearch[3 * irecord + 1];
    blockid = isearch[3 * irecord + 2];
    meshtagrecv = tagsearch[irecord];
    
    int id = d_dim * irecord;
    for (int i = 0; i < d_dim; i++)
    {
        xp[i] = xsearch[id + i];
    }
    //
    int i;
    int isum = 0;
    int n;
    for (n = 0; n < ntypes; n++)
    {
        isum += nc[n];
        if (donorId[irecord] < isum)
        {
            i = donorId[irecord] - (isum - nc[n]);
            break;
        }
    }
    int nvert = nv[n];
    acceptFlag = 1;
    for (int m = 0; m < nvert; m++)
    {
        inode[m] = vconn[n][nvert * i + m] - BASE;
        id = d_dim * inode[m];
        if (iblank[inode[m]] <= 0 && receptorRes2 > 0.0)//If this donor has any node that is another's receptor and receptorRes2 >0(not from a scaled mesh)
        {
            if (nodeRes[inode[m]] == BIGVALUE) acceptFlag = 0;//If is mandatory or interpolated from where the receptor is from
            if (abs(iblank[inode[m]]) == meshtagrecv) acceptFlag = 0;
        }
        for (int j = 0; j < d_dim; j++)
            xv[m][j] = x[id + j];
    }

    if (receptorRes == BIGVALUE && resolutionScale == 1.0)//If the remote receptor is a mandatory ,then cancel all the unmandatory vertices of this donor
    {
        clist = cancelList;
        //
        // go to the end of the list 
        //
        if (clist != NULL) while (clist->next != NULL) clist = clist->next;
        //
        for (int m = 0; m < nvert; m++)
        {
            inode[m] = vconn[n][nvert * i + m] - BASE;
            if (iblank[inode[m]] <= 0 and nodeRes[inode[m]] != BIGVALUE)
            {
                if (iblank[inode[m]] < 0) iblank[inode[m]] = 1;
                if (clist == NULL)
                {
                    clist = new INTEGERLIST[1];
                    clist->inode = inode[m];
                    clist->next = NULL;
                    cancelList = clist;
                }
                else
                {
                    clist->next = new INTEGERLIST[1];
                    clist->next->inode = inode[m];
                    clist->next->next = NULL;
                    clist = clist->next;
                }
                ncancel++;
                cancel_count++;
            }
        }
    }
    //  
    computeNodalWeights(d_dim, xv, xp, frac, nvert);
    //
    interp2donor[irecord] = *recid;//map from index inside nsearch to interp index
    interpList[*recid].cancel = 0;
    interpList[*recid].nweights = nvert;
    interpList[*recid].receptorInfo[0] = procid;
    interpList[*recid].receptorInfo[1] = pointid;
    interpList[*recid].receptorInfo[2] = blockid;

    interpList[*recid].inode = new int[nvert + 1];
    interpList[*recid].weights = new double[nvert + 1];
    int m;
    for (m = 0; m < nvert; m++)
    {
        
        interpList[*recid].inode[m] = inode[m];
        interpList[*recid].weights[m] = frac[m];
        if (frac[m] < -0.2 || frac[m] > 1.2) 
        {
            throw std::runtime_error("Frac Problem Again? Oh come on");
            TRACEI(myid);
            TRACEI(irecord);
            TRACEI(meshtag);
            TRACEI(donorId[irecord]);
            TRACED(frac[m]);
            //MPI_Abort(MPI_COMM_WORLD, m);
        }
    }
    interpList[*recid].inode[m] = donorId[irecord];
    interpList[*recid].weights[m] = 0.0;
    if (acceptFlag == 0 && receptorRes != BIGVALUE) interpList[*recid].cancel = 1;//This interp is cancelled but the iblank is not set to 1 yet
    (*recid)++;
}

void MeshBlock::set_ninterp(int ninterp_input)
{
    ninterp = ninterp_input;
}

void MeshBlock::getCancellationData(int* nrecords, int** intData)
{
    int i;
    int inode;
    INTEGERLIST* clist;
    *nrecords = ncancel;
    if (ncancel > 0)
    {
        *intData = new int[3 * (*nrecords)];
        
        i = 0;
        for (clist = cancelList; clist != NULL; clist = clist->next)
        {
            inode = clist->inode;
            if (donorList[inode] != NULL) {
                (*intData)[i++] = donorList[inode]->donorData[0];
                (*intData)[i++] = donorList[inode]->donorData[2];
                (*intData)[i++] = donorList[inode]->donorData[1];
            }
        }
        *nrecords = i / 3;
    }
}

void MeshBlock::getCancellationReduce(int* nrecords, int** intData)
{
    int i;
    int inode;
    
    *nrecords = ncancel;
    if (ncancel > 0)
    {
        *intData = new int[3 * (*nrecords)];
        
        i = 0;
        for(int j = 0;j<ncancel;j++)
        {
            inode = cancelArray[j];
            if (donorList[inode] != NULL) {
                (*intData)[i++] = donorList[inode]->donorData[0];
                (*intData)[i++] = donorList[inode]->donorData[2];
                (*intData)[i++] = donorList[inode]->donorData[1];
            }
        }
        *nrecords = i / 3;
    }
}

void MeshBlock::cancelDonor(int irecord)
{
    int iptr;
    iptr = interp2donor[irecord];
    if (iptr > -1) interpList[iptr].cancel = 1;
}


void MeshBlock::resetCoincident()//Coincident means that one point shared by two splitting
{
    int i, iptr;
    int* ireset = new int[nsearch]; 
    for (i = 0; i < nsearch; i++) ireset[i] = 1;

    //int coincident_count = 0;
    for (i = 0; i < nsearch; i++)
    {
        iptr = interp2donor[i];
        if (iptr > -1)
        {
            ireset[xtag[i]] = TIOGA_MIN(ireset[xtag[i]], interpList[iptr].cancel);
        }
    }
    for (i = 0; i < nsearch; i++)
    {
        iptr = interp2donor[i];
        if (iptr > -1) 
        {
            ireset[xtag[i]] = TIOGA_MIN(ireset[xtag[i]], interpList[iptr].cancel);
        }
    }
    //Why this happens ?
    //two coincident points must have been of same resolution, but the two cells owning it may be different
    //But either case this should have been treated as usual, that is a "and" relationship, but here is a or
    //I guess the author implies that on the edge a interpolation comes better than calculation.

    for (i = 0; i < nsearch; i++)
    {
        iptr = interp2donor[i];
        if (iptr > -1) {
            if (interpList[iptr].cancel == 1) 
            {
                interpList[iptr].cancel = ireset[xtag[i]];//Either one of two coincident search point is not cancelled, leave them all
            }
        }
    }
    TIOGA_FREE(ireset);
}
void MeshBlock::getInterpData(int* nrecords, int** intData)
{
    int i, k;
    //
    *nrecords = 0;
    for (i = 0; i < ninterp; i++)
        if (!interpList[i].cancel) (*nrecords)++;
    //
    (*intData) = new int[3 * (*nrecords)];
    for (i = 0, k = 0; i < ninterp; i++)
        if (!interpList[i].cancel) 
        {
            (*intData)[k++] = interpList[i].receptorInfo[0];
            (*intData)[k++] = interpList[i].receptorInfo[1];
            (*intData)[k++] = interpList[i].receptorInfo[2];
        }
}

void MeshBlock::clearIblanks(void)
{
    int i;
    for (i = 0; i < nnodes; i++)
        if (iblank[i] < 0) iblank[i] = 1;
    if (iblank_reduced) //%TODO::We will come to this later after iblank_reduce is made.
    {
        for (i = 0; i < nnodes; i++)
        {
            if (iblank_reduced[i] == 0)
            {
                iblank[i] = 0;
            }
        }
        TIOGA_FREE(iblank_reduced);
    }
}

void MeshBlock::setIblanks(int inode)
{
    iblank[inode] = -1;
}

void MeshBlock::reduce_fringes()
{
    
    int inode[8];
    INTEGERLIST* clist;
    INTEGERLIST* clist_next;

    if (iblank_reduced) TIOGA_FREE(iblank_reduced);
    iblank_reduced = new int[nnodes];
    int* ibltmp = new int[nnodes];
    for (int i = 0; i < nnodes; i++)
    {
        iblank_reduced[i] = iblank[i] > 0 ? iblank[i] : 0;//Reserve only the positive ones, which are not receptors
    }
    for (int i = 0; i < nnodes; i++)
    {
        ibltmp[i] = iblank_reduced[i];//ibltmap is the initial all-positive-one or zero
    }

    for (int iter = 0; iter < nfringe + 1; iter++)
    {
        // writeGridFile2("iblankReduced_"+std::to_string(iter), iblank_reduced);
        for(int n=0;n<ntypes;n++)
        {
            int nvert = nv[n];
            
            for (int i = 0; i < nc[n]; i++)
            {
                int ncount = 0;
                for (int m = 0; m < nvert; m++)
                {
                    inode[m] = vconn[n][nvert * i + m] - BASE;
                    ncount += (iblank_reduced[inode[m]] == 1 or iblank_reduced[inode[m]] < 0);//Not zero
                }
                if (ncount > 0 and ncount < nvert)//There is a node in holeMap
                {
                    for (int m = 0; m < nvert; m++)
                    {
                        if (ibltmp[inode[m]] == 0 and iblank[inode[m]] < 0)
                        {
                            ibltmp[inode[m]] = iblank[inode[m]];
                        }
                    }
                }
            }
        }
        for (int i = 0; i < nnodes; i++)
        {
            iblank_reduced[i] = ibltmp[i];
        }
    }
    // writeGridFile2("iblankReduced_"+std::to_string(nfringe+1), iblank_reduced);
    if (cancelList) deallocateLinkList2(cancelList);
    cancelList = NULL;
    ncancel = 0;
    clist = cancelList;
    int count = 0;
    for (int i = 0; i < nnodes; i++)
    {
        if (iblank[i] < 0 and iblank_reduced[i] == 0)
        {
            ncancel++;
        }
    }
    if (ncancel > 0)
    {
        int count = 0;
        cancelArray = new int[ncancel];
        for (int i = 0; i < nnodes; i++)
        {
            if (iblank[i] < 0 and iblank_reduced[i] == 0)
            {
                cancelArray[count++] = i;
            }
        }
    }
    
    // writeGridFile2("iblankReduced", iblank_reduced);
    // writeGridFile2("iblank_nReduced", iblank);
    TIOGA_FREE(ibltmp);
}

void MeshBlock::getCellIblanks(void)
{
    int i;
    int n, nvert, m;
    int icell;
    int inode[8];
    int ncount, flag;
    int verbose;
    int* ibl;

    if (iblank_reduced)
    {
        // std::cout<<"Is reduced\n";
        // std::cin.get();
        ibl = iblank_reduced;
    }
    else
    {
        // std::cout<<"Is NOT reduced\n";
        // std::cin.get();
        ibl = iblank;
    }
    icell = 0;
    if (iblank_cell == NULL) iblank_cell = new int[ncells];
    for (n = 0; n < ntypes; n++)
    {
        nvert = nv[n];
        for (i = 0; i < nc[n]; i++)
        {
            flag = 1;
            iblank_cell[icell] = 1;
            verbose = 0;
            ncount = 0;
            for (m = 0; m < nvert && flag; m++)
            {
                inode[m] = vconn[n][nvert * i + m] - BASE;

                if (ibl[inode[m]] == 0)
                {
                    iblank_cell[icell] = 0;
                    flag = 0;
                }
                ncount = ncount + (ibl[inode[m]] == -1);
            }
            if (flag)
            {
                if (ncount == nvert) iblank_cell[icell] = -1;
            }
            icell++;
        }
    }
}

void MeshBlock::getCellIblanks2()
{
    int i;
    int n, nvert, m;
    int icell;
    int inode[8];
    int ncount, flag;
    int verbose;

    icell = 0;
    if (iblank_cell == NULL) 
        iblank_cell = new int[ncells];

    for (n = 0; n < ntypes; n++)
    {
        nvert = nv[n];
        for (i = 0; i < nc[n]; i++)
        {
            flag = 1;
            iblank_cell[icell] = 1;
            verbose = 0;
            //if (myid==763 && icell==7975) {
            // verbose=1;
            //}
            ncount = 0;
            for (m = 0; m < nvert && flag; m++)
            {
                inode[m] = vconn[n][nvert * i + m] - BASE;
                if (iblank[inode[m]] == 0)//There is a wall, then it is wall
                {
                    iblank_cell[icell] = 0;
                    flag = 0;
                }
                ncount = ncount + (iblank[inode[m]] == -1);
            }
            if (flag)//Not a wall and all node's iblank is -1
            {
                if (ncount == nvert) iblank_cell[icell] = -1;
                //            if (ncount > 0)  iblank_cell[icell]=-1;
                //            if (ncount >= nvert/2) iblank_cell[icell]=-1;
            }
            icell++;
        }
    }
}