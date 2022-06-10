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
#include "codetypes.h"
#include "MeshBlock.h"
#include <unordered_map>
#include <iostream>
#include"ADT.h"
bool real_TriangleInclusionTest(double* vp, double xv[3][3]);
void computeNodalWeight2D(double xv[8][3], double* xp, double* frac, int nvert);
void findOBB(int dim,double* x, double xc[3], double dxc[3], double vec[3][3], int nnodes);
//void writebbox(OBB* obb, int bid);
//void writePoints(double* x, int nsearch, int bid);
void uniquenodes(double* x, int* meshtag, double* rtag, int* itag, int* nn);
void uniquenodes_octree(int dim, double* x, int* meshtag, double* rtag, int* itag, int* nn);


namespace {

    /** Determine the unique nodes by a global identifier
     *
     *  The function will create a mapping such that all duplicate nodes will point
     *  to the original node (as determined by a global identifier) in the `itag`
     *  array. It will update the nodal resolutions the shared nodes such that the
     *  resolutions upon exit will be the maximum resolution amongst all the
     *  duplicate nodes.
     *
     *  \param[in] node_ids Global IDs for the nodes across all MPI ranks
     *  \param[inout] node_res The nodal resolutions
     *  \param[out] itag The local index of the original node (duplicate to original mapping)
     *  \param[in] nnodes The size of the arrays
     */
    void uniquenode_map(uint64_t* node_ids, double* node_res, int* itag, int nnodes)
    {
        std::unordered_map<uint64_t, int> lookup;

        for (int i = 0; i < nnodes; i++) {
            auto found = lookup.find(node_ids[i]);
            if (found != lookup.end()) {
                // This is a duplicate node, store the index to the original node
                // found previously
                itag[i] = found->second;

                // Update the original node's resolution to be the max of either
                // node resolution
                node_res[found->second] = std::max(node_res[found->second], node_res[i]);
            }
            else {
                // This is the first appearance of the unique ID, stash it in the
                // lookup table
                lookup[node_ids[i]] = i;
                itag[i] = i;
            }
        }

        // The max node resolution was stored off in the original node, propagate
        // this to all the duplicates
        for (int i = 0; i < nnodes; i++)
            node_res[i] = node_res[itag[i]];
    }
}


#include<fstream>
void MeshBlock::search()
{
    
    // std::string filename = "FoundQueries" + std::to_string(myid) + std::to_string(meshtag) + ".sct";
    // std::ofstream fout;
    // fout.open(filename);
    int ndim, id;
    int iptr, isum, nvert;
    OBB* obq;
    int* icell = new int[ncells];

    int cell_count;
    int cellindex;
    double xd[3];
    double dxc[3];
    double xmin[3];
    double xmax[3];
    int* dId;
    bool cross = true;
    //
    // form the bounding box of the 
    // query points
    //
    if (nsearch == 0) {
        donorCount = 0;
        return;
    }
    obq = new OBB[1];
    findOBB(d_dim, xsearch, obq->xc, obq->dxc, obq->vec, nsearch);
    writeOBB2("OBQ", obq);
    // find all the cells that may have intersections with
    // the OBB
    for (int i = 0; i < ncells; i++) { icell[i] = -1; }
    iptr = -1;
    cell_count = 0;
    int p = 0;
    for (int n = 0; n < ntypes; n++)
    {
        int nvert = nv[n];
        for (int i = 0; i < nc[n]; i++)
        {
            //
            // find each cell that has
            // overlap with the bounding box
            //
            xmin[0] = xmin[1] = xmin[2] = BIGVALUE;
            xmax[0] = xmax[1] = xmax[2] = -BIGVALUE;
            for (int m = 0; m < nvert; m++)
            {
                id = d_dim * (vconn[n][nvert * i + m] - BASE);
                for (int j = 0; j < d_dim; j++)
                {
                    xd[j] = 0;
                    for (int k = 0; k < d_dim; k++)
                        xd[j] += (x[id + k] - obq->xc[k]) * obq->vec[j][k];
                    xmin[j] = TIOGA_MIN(xmin[j], xd[j]);
                    xmax[j] = TIOGA_MAX(xmax[j], xd[j]);
                }
                for (int j = 0; j < d_dim; j++)
                {
                    xd[j] = (xmax[j] + xmin[j]) * 0.5;
                    dxc[j] = (xmax[j] - xmin[j]) * 0.5;
                }
            }
            cross = true;
            for (int j = 0; j < d_dim; j++)
            {
                cross = cross and (fabs(xd[j]) <= dxc[j] + obq->dxc[j]);

            }
            if (cross)
            {
                //iblank_cell[p] = 1;
                icell[p] = iptr;//icell starting from 0;
                iptr = p;
                cell_count++;
            }
            p++;
        }
    }

    //
    // now find the axis aligned bounding box
    // of each cell in the LIFO stack to build the
    // ADT
    //
    
    if (elementBbox) TIOGA_FREE(elementBbox);
    if (elementList) TIOGA_FREE(elementList);
    elementBbox = new double[2 * d_dim * cell_count];//BoundingBox upper range and lower
    elementList = new int[cell_count];//cellIndex within the box.
    int k = iptr;
    int i;
    int n;
    int l = 0;
    //for(k=0;k<ncells;k++)
    while (k != -1)
    {
        cellindex = k;
        isum = 0;
        for (n = 0; n < ntypes; n++)
        {
            isum += nc[n];
            if (cellindex < isum)
            {
                i = cellindex - (isum - nc[n]);
                break;
            }
        }
        nvert = nv[n];
        xmin[0] = xmin[1] = xmin[2] = BIGVALUE;
        xmax[0] = xmax[1] = xmax[2] = -BIGVALUE;
        for (int m = 0; m < nvert; m++)
        {
            id = d_dim * (vconn[n][nvert * i + m] - BASE);
            for (int j = 0; j < d_dim; j++)
            {
                xmin[j] = TIOGA_MIN(xmin[j], x[id + j]);
                xmax[j] = TIOGA_MAX(xmax[j], x[id + j]);
            }
        }
        //
        for (int q = 0; q < d_dim; q++)
        {

            elementBbox[(2 * l) * d_dim + q] = xmin[q];
            elementBbox[(2 * l + 1) * d_dim + q] = xmax[q];
            //std::cout <<"dim "<<q+1<<" " << elementBbox[(2 * l) * d_dim + q] << " " << elementBbox[(2 * l + 1) * d_dim + q] << '\n';
        }
        elementList[l] = k;
        k = icell[k];
        l++;
    }
    
    //
    // build the ADT now
    //
    
    if (adt)
    {
        adt->clearData();
    }
    else
    {
        adt = new ADT[1];
        adt->setDim(d_dim);
    }
    
    ndim = 2 * d_dim;
    //
    std::cout << "Building ADT\n";
    adt->buildADT(ndim, cell_count, elementBbox);
    std::cout << "Allocating donorID\n";
    //
    if (donorId) TIOGA_FREE(donorId);
    donorId = new int[nsearch];
    if (xtag) TIOGA_FREE(xtag);
    xtag = new int[nsearch];
    
    //
    // create a unique hash
    //
#ifdef TIOGA_HAS_NODEGID
    uniquenode_map(gid_search.data(), res_search, xtag, nsearch);
#else
    uniquenodes_octree(d_dim, xsearch, tagsearch, res_search, xtag, &nsearch);
#endif

    dId = new int[2];

    for (int i = 0; i < nsearch; i++)
    {
        if (xtag[i] == i) {

            //adt->searchADT(this,&(donorId[i]),&(xsearch[3*i]));
            adt->searchADT(this, dId, &(xsearch[d_dim * i]));
            donorId[i] = dId[0];
        }
        else
        {
            donorId[i] = donorId[xtag[i]];
        }
    }
    
    TIOGA_FREE(dId);
    TIOGA_FREE(icell);
    TIOGA_FREE(obq);
    std::cout << "Returing from search\n";
    return;
}

void MeshBlock::search_distance(ADT* cur_adt,double* distance,int WorO)
{
    int* iflag = new int[nnodes];
    for (int i = 0; i < nnodes; i++)
    {
        iflag[i] = 1;
    }
    if (WorO == 1)//wall
    {
        for (int i = 0; i < nwbc; i++)
        {
            iflag[wbcnode[i] - 1] = -1;
            distance[wbcnode[i] - 1] = 0;
        }
        if (nwbc == 0)
        {
            double minD = BIGVALUE;
            double tempD = 0;
            iflag[0] = -1;
            //Find a random point(or points if this could be improved improve)
            for (int j = 0; j < ngWbc; j++)
            {
                tempD = 0;
                for (int k = 0; k < d_dim; k++)
                {
                    tempD += pow(x[k] - gWbcXYZ[d_dim * j + k], 2);
                }
                tempD = sqrt(tempD);
                minD = TIOGA_MIN(minD, tempD);
            }
            distance[0] = minD;
        }
    }
    else if (WorO == 0)
    {
        for (int i = 0; i < nobc; i++)
        {
            //std::cout << obcnode[i] << '\n';
            iflag[obcnode[i] - 1] = -1;
            distance[obcnode[i] - 1] = 0;
        }
        if (nobc == 0)
        {
            double minD = BIGVALUE;
            double tempD = 0;
            iflag[0] = -1;
            //Find a random point(or points if this could be improved improve)
            for (int j = 0; j < ngObc; j++)
            {
                tempD = 0;
                for (int k = 0; k < d_dim; k++)
                {
                    tempD += pow(x[k] - gObcXYZ[d_dim * j + k], 2);
                }
                tempD = sqrt(tempD);
                minD = TIOGA_MIN(minD, tempD);
            }
            distance[0] = minD;
        }
    }
    bool signal = true;
    int nvert;
    int nchecked;
    int checked[8];
    int ichecked;
    double d2;
    int i1;
    int i2;
    int tag_count = 0;
    for (int i = 0; i < nnodes; i++)
    {
        if (iflag[i] == -1)
        {
            tag_count++;
        }
    }
    int iter_count = 0;
    while (signal and iter_count<10)
    {
        iter_count++;
        std::cout << "Inside loop\n";
        std::cout << signal << '\n';
        signal = false;
        for (int n = 0; n < ntypes; n++)
        {
            nvert = nv[n];
            for (int i = 0; i < nc[n]; i++)
            {
                nchecked = 0;
                for (int m = 0; m < nvert; m++)
                {
                    if (iflag[vconn[n][i*nvert+m]-BASE] == -1)//If there is a vertex already found 
                    {
                        checked[m] = 1;
                        ++nchecked;
                        ichecked = m;
                    }
                    else//There is a vertex not found yet
                    {
                        checked[m] = 0;
                        signal = true;
                    }
                }
                //std::cout << nchecked << " checked out of " << nvert << '\n';
                if (nchecked > 0 and nchecked < nvert)//There is a node to search
                {
                    for (int m = 0; m < nvert; m++)
                    {
                        if (checked[m] == 0)
                        {
                            i1 = vconn[n][i * nvert + ichecked] - BASE;
                            i2 = vconn[n][i * nvert + m] - BASE;
                            d2 = 0;
                            //search with the already existing node ichecked
                            for (int j = 0; j < d_dim; j++)
                            {
                                
                                d2 += pow(x[d_dim * i1 + j] - x[d_dim * i2 + j], 2);
                            }
                            d2 = sqrt(d2) + distance[i1];
                            cur_adt->searchNearest(&d2, &(x[d_dim * i2]),myid);
                            distance[i2] = d2;
                            iflag[i2] = -1;
                        }
                        //Now call ADT to search , adt need to know coordinates of x2 and searchBox, which will be updated as search goes on.
                        
                    }
                }
            }
        }
        tag_count = 0;
        for (int i = 0; i < nnodes; i++)
        {
            if (iflag[i] == -1)
            {
                tag_count++;
            }
        }
        std::cout << "Tag count is " << tag_count << "out of " << nnodes << " at proc " << myid << " meshtag " << meshtag << '\n';
        std::cout << ngWbc << " " << ngObc << '\n';
    }
    delete[] iflag;
}

void MeshBlock::set_global_bc(int inwbc, int inobc, double* iwbcXYZ, double* iobcXYZ)
{
    ngWbc = inwbc;
    ngObc = inobc;
    gWbcXYZ = iwbcXYZ;
    gObcXYZ = iobcXYZ;

}

void MeshBlock::get_distances()//We might only want this to work when there are two bc at the same time
{
    if (wallDistance) { TIOGA_FREE(wallDistance); }
    if (overDistance) { TIOGA_FREE(overDistance); }

    int ndim = 2 * d_dim;
    double* wallElementBox = NULL;
    double* overElementBox = NULL;
    ADT* wallADT = new ADT[1];
    ADT* overADT = new ADT[1];
    if (ngWbc != 0)
    {
        wallDistance = new double[nnodes];
        wallElementBox = new double[2 * d_dim * ngWbc];
    }
    if (ngObc != 0)
    {
        overDistance = new double[nnodes];
        overElementBox = new double[2 * d_dim * ngObc];
    }
    for (int i = 0; i < ngWbc; i++)
    {
        for (int j = 0; j < d_dim; j++)
        {
            wallElementBox[(2 * i) * d_dim + j] = wallElementBox[(2 * i + 1) * d_dim + j] = gWbcXYZ[d_dim * i + j];
        }
    }
    for (int i = 0; i < ngObc; i++)
    {
        for (int j = 0; j < d_dim; j++)
        {
            overElementBox[(2 * i) * d_dim + j] = overElementBox[(2 * i + 1) * d_dim + j] = gObcXYZ[d_dim * i + j];
        }
    }
    if (ngWbc != 0)
    {
        
        wallADT->buildADT(ndim, ngWbc, wallElementBox);
        search_distance(wallADT, wallDistance, 1);
    }
    if (ngObc != 0)
    {
        overADT->buildADT(ndim, ngObc, overElementBox);
        //This is time-consuming
        //search_distance(overADT, overDistance, 0);
    }
    if (wallElementBox) { TIOGA_FREE(wallElementBox); }
    if (overElementBox) { TIOGA_FREE(overElementBox); }
    delete[] wallADT;
    delete[] overADT;
}

void MeshBlock::reResolution()
{
    if (wallDistance == NULL or overDistance == NULL)
    {
        return;
    }
    double* nodeEff = new double[nnodes];
    double* cellEff = new double[ncells];
    for (int i = 0; i < nnodes; i++)
    {
        nodeEff[i] = wallDistance[i] / (wallDistance[i] + overDistance[i]);
        if (nodeRes[i] == BIGVALUE)
        {

        }
        else
        {
            nodeRes[i] = nodeRes[i] * nodeEff[i];
        }
    }
    int temp_int = 0;
    int nvert;
    for (int n = 0; n < ntypes; n++)
    {
        nvert = nv[n];
        for (int i = 0; i < nc[n]; i++)
        {
            cellEff[temp_int] = 0;
            for (int m = 0; m < nvert; m++)
            {
                cellEff[temp_int] += nodeEff[vconn[n][i * nvert + m] - BASE];
            }
            cellEff[temp_int] = cellEff[temp_int] / nvert;
            temp_int++;
        }
    }
    for (int i = 0; i < ncells; i++)
    {
        if (cellRes[i] == BIGVALUE)
        {

        }
        else
        {
            cellRes[i] = cellRes[i] * cellEff[i];
        }
    }
    writeGridFile2("nodeweight", nodeEff);
    writeCellFile2("cellweight", cellEff);
    delete[] nodeEff;
    delete[] cellEff;
}