#include"MeshBlock.h"
#include<iostream>
#include<fstream>

//double computeCellVolume(double xv[8][3], int nvert);

double computeVolume(int dim, double xv[8][3], int nvert, bool counter_clock = false);
void computeNodalWeights(int dim, double xv[8][3], double* xp, double* frac, int nvert);
void deallocateLinkList(DONORLIST* temp);
void deallocateLinkList2(INTEGERLIST* temp);

//Me abandoned cuz it doesn't seem like we could have a uniform hex
//double tdot_product(double a[3], double b[3], double c[3]);
void findOBB(int dim, double* x, double xc[3], double dxc[3], double vec[3][3], int nnodes);
//void findOBB(double* x, double xc[3], double dxc[3], double vec[3][3], int nnodes);
void getobbcoords(int dim, double* xc, double* dxc, double vec[3][3], double xv[8][3]);
void transform2OBB(int dim, double* xv, double* xc, double vec[3][3], double* xd);
//The author abandoned
//void writebbox(OBB* obb, int bid);
//void writebboxdiv(OBB* obb, int bid);




MeshBlock::~MeshBlock()
{
	int i;
    if (cellRes) TIOGA_FREE(cellRes);
    if (nodeRes) TIOGA_FREE(nodeRes);
    if (elementBbox) TIOGA_FREE(elementBbox);
    if (elementList) TIOGA_FREE(elementList);
    //if (adt) delete[] adt;
    if (donorList) {
        for (i = 0; i < nnodes; i++) deallocateLinkList(donorList[i]);
        TIOGA_FREE(donorList);
    }
    if (interpList) {
        for (i = 0; i < interpListSize; i++)
        {
            if (interpList[i].inode) TIOGA_FREE(interpList[i].inode);
            if (interpList[i].weights) TIOGA_FREE(interpList[i].weights);
        }
        TIOGA_FREE(interpList);
    }
    if (interpList2) {
        for (i = 0; i < interp2ListSize; i++)
        {
            if (interpList2[i].inode) TIOGA_FREE(interpList2[i].inode);
            if (interpList2[i].weights) TIOGA_FREE(interpList2[i].weights);
        }
        TIOGA_FREE(interpList2);
    }
    if (interpListCart) {
        for (i = 0; i < interpListCartSize; i++)
        {
            if (interpListCart[i].inode) TIOGA_FREE(interpListCart[i].inode);
            if (interpListCart[i].weights) TIOGA_FREE(interpListCart[i].weights);
        }
        TIOGA_FREE(interpListCart);
    }
    // For nalu-wind API the iblank_cell array is managed on the nalu side
    // if (!ihigh) {
    //  if (iblank_cell) TIOGA_FREE(iblank_cell);
    // }
    if (obb) TIOGA_FREE(obb);
    if (obh) TIOGA_FREE(obh);
    if (isearch) TIOGA_FREE(isearch);
    if (xsearch) TIOGA_FREE(xsearch);
    if (res_search) TIOGA_FREE(res_search);
    if (xtag) TIOGA_FREE(xtag);
    if (rst) TIOGA_FREE(rst);
    if (interp2donor) TIOGA_FREE(interp2donor);
    if (cancelList) deallocateLinkList2(cancelList);
    if (ctag) TIOGA_FREE(ctag);
    if (pointsPerCell) TIOGA_FREE(pointsPerCell);
    if (rxyz) TIOGA_FREE(rxyz);
    if (picked) TIOGA_FREE(picked);
    if (rxyzCart) TIOGA_FREE(rxyzCart);
    if (donorIdCart) TIOGA_FREE(donorIdCart);
    if (pickedCart) TIOGA_FREE(pickedCart);
    if (ctag_cart) TIOGA_FREE(ctag_cart);

    if (tagsearch) TIOGA_FREE(tagsearch);
    if (donorId) TIOGA_FREE(donorId);
    if (receptorIdCart) TIOGA_FREE(receptorIdCart);
    if (icft) TIOGA_FREE(icft);
    if (mapmask) TIOGA_FREE(mapmask);
    if (uindx) TIOGA_FREE(uindx);
    if (invmap) TIOGA_FREE(invmap);
    // need to add code here for other objects as and
    // when they become part of MeshBlock object  
}


void MeshBlock::setData(int dim, int btag, int nnodesi, double* xyzi, int* ibli, int nwbci, int nobci,
    int* wbcnodei, int* obcnodei,
    int ntypesi, int* nvi, int* nci, int** vconni,
    uint64_t* cell_gid, uint64_t* node_gid)
{
    d_dim = dim;
    meshtag = btag;
    nnodes = nnodesi;
    x = xyzi;
    iblank = ibli;
    nwbc = nwbci;
    nobc = nobci;
    wbcnode = wbcnodei;
    obcnode = obcnodei;

    ntypes = ntypesi;

    nv = nvi;
    nc = nci;
    vconn = vconni;
    cellGID = cell_gid;
    nodeGID = node_gid;

    ncells = 0;
    for (int i = 0; i < ntypes; i++)
    {
        ncells += nc[i];
    }
#ifdef  TIOGA_HAS_NODEGID
    if (nodeGID == NULL)
    {
        throw std::runtime_error("tioga: global IDs for nodes not provides\n");
    }
#endif //  TIOGA_HAS_NODEGID



}


void MeshBlock::preprocess()
{
    writeWbc("raw_WBC");
    writeObc("raw_OBC");
    for (int i = 0; i < nnodes; i++)
    {
        iblank[i] = 1;
    }
    //%TODO: There is a uniform hex checking
    /*if (check_uniform_hex_flag)
    {
        check_for_uniform_hex();
        if (uniform_hex) create_hex_cell_map();
    }*/
    if (obb) TIOGA_FREE(obb);
    obb = new OBB[1];
    //std::cout << d_dim << " CALLING findOBB\n";
    findOBB(d_dim,x, obb->xc, obb->dxc, obb->vec, nnodes);//%TODO::Here vector is 2 2DArray not a pointer's pointer
    writeOBB("rawOBB");
    tagBoundary();
}


void MeshBlock::tagBoundary()
{
    //Local arraies to be used
    int inode[8];
    double xv[8][3];

    double vol;
    int* iflag = new int[nnodes];
    int* iextmp = new int[nnodes];
    int* iextmp1 = new int[nnodes];
    int* iptr = new int[nnodes];

    double xd[3], xmin[3], xmax[3];
    int idx[3];
    int kc = 0;

    //%TODO::This could be combined into my implemention of OBB
    int dim_vol = 1;
    for (int i = 0; i < d_dim; i++)
    {
        mapdims[i] = 10;
        dim_vol = dim_vol * mapdims[i];
        mapdx[i] = 2 * obb->dxc[i] / mapdims[i];
    }
    int frac[3] = { 1 ,1,1 };
    for (int i = 1; i < d_dim; i++)
    {
        frac[i] = frac[i - 1] * mapdims[i - 1];
    }

    if (cellRes) TIOGA_FREE(cellRes);
    cellRes = new double[ncells];
    if (nodeRes) TIOGA_FREE(nodeRes);
    nodeRes = new double[nnodes];

    if (icft) TIOGA_FREE(icft);
    icft = new int[dim_vol + 1];
    if (invmap) TIOGA_FREE(invmap);
    invmap = new int[nnodes];
    if (mapmask) TIOGA_FREE(mapmask);
    mapmask = new int[dim_vol];
    for (int i = 0; i < dim_vol + 1; i++)
    {
        icft[i] = -1;
    }
    icft[0] = 0;
    
    for (int i = 0; i < nnodes; i++)
    {
        iflag[i] = 0;
    }
    if(userSpecifiedNodeRes!=NULL and userSpecifiedCellRes!=NULL)
    {
        int k = 0;
        for (int n = 0; n < ntypes; n++)
        {
            for (int i = 0; i < nc[n]; i++)
            {
                cellRes[k] = userSpecifiedCellRes[k];
                k++;
            }
        }
        for (k = 0; k < nnodes; k++)
        {
            nodeRes[k] = userSpecifiedNodeRes[k];
        }
    }
    else
    {
        for (int i = 0; i < nnodes; i++)
        {
            nodeRes[i] = 0.0;
        }
        int k = 0;
        int nvert;
        int ni;
        // std::ofstream fout;
        // fout.open("FoundNegative_"+std::to_string(meshtag)+"_"+std::to_string(myid));
        for (int n = 0; n < ntypes; n++)
        {
            nvert = nv[n];
            for(int i=0;i<nc[n];i++)
            {
                ni = nvert * i;
                for (int m = 0; m < nvert; m++)
                {
                    inode[m] = vconn[n][ni + m] - BASE;
                    for (int j = 0; j < d_dim; j++)
                    {
                        xv[m][j] = x[d_dim * inode[m] + j];
                    }
                }
                vol = computeVolume(d_dim, xv, nvert);
                // if(vol <=0 )
                // {
                //     fout<<"Cell: type "<<n<<" with "<<nvert<<" verts with volume being "<<vol<<"\n";
                //     for(int m = 0;m<nvert;m++)
                //     {
                //         fout<<"\tnode "<<m<<": "<<xv[m][0]<<","<<xv[m][1]<<","<<xv[m][2]<<'\n';
                //     }
                // }
                //%TODO::Resolution scale is always set to 1
                cellRes[k] = (vol * resolutionScale);
                for (int m = 0; m < nvert; m++)
                {
                    ++iflag[inode[m]];
                    nodeRes[inode[m]] += cellRes[k];
                }
                ++k;
            }
            
            
        }
        // fout.close();

    }
    writeCellFile2("cellResImmediate", cellRes);

    //Make inverse map
    for (int i = 0; i < nnodes; i++)
    {
        if (iflag[i] != 0) nodeRes[i] = nodeRes[i] / iflag[i];
        iflag[i] = 0;
        iextmp[i] = iextmp1[i] = 0;
        for (int j = 0; j < d_dim; j++)
        {
            xd[j] = obb->dxc[j];
            for (int k = 0; k < d_dim; k++)
            {
                xd[j] += (x[d_dim * i + k] - obb->xc[k]) * obb->vec[j][k];
            }
            idx[j] = static_cast<int>(xd[j] / mapdx[j]);
        }
        int indx = 0;
        for (int i = 0; i < d_dim; i++)
        {
            indx += idx[i] * frac[i];
        }
        iptr[i] = icft[indx + 1];
        icft[indx + 1] = i;
    }
    for (int i = 0; i < dim_vol; i++)
    {
        int ip = icft[i + 1];
        int m = 0;
        while (ip != -1)
        {
            invmap[kc++] = ip;
            ip = iptr[ip];
            m++;
        }
        icft[i + 1] = icft[i] + m;
    }
    TIOGA_FREE(iptr);
    //Tag overset boundary nodes
    for (int i = 0; i < nobc; i++)
    {
        iflag[(obcnode[i] - BASE)] = 1;
    }
    //Make mapmask while tagging over boundary cells
    int itag;
    for (int n = 0; n < ntypes; n++)
    {
        int nvert = nv[n];
        for (int i = 0; i < nc[n]; i++)
        {
            itag = 0;
            for (int j = 0; j < d_dim; j++) { xmin[j] = BIGVALUE; xmax[j] = -BIGVALUE; }
            for (int m = 0; m < nvert; m++)
            {

                inode[m] = vconn[n][nvert * i + m] - BASE;
                if (iflag[inode[m]]) itag = 1;//If is  obc ,set itag to 1
                for (int j = 0; j < d_dim; j++)
                {
                    xd[j] = obb->dxc[j];
                    for (int k = 0; k < d_dim; k++)
                    {
                        xd[j] += (x[d_dim * inode[m] + k] - obb->xc[k]) * obb->vec[j][k];
                    }
                    xmin[j] = TIOGA_MIN(xd[j], xmin[j]);
                    xmax[j] = TIOGA_MAX(xd[j], xmax[j]);
                }
            }
            for (int j = 0; j < d_dim; j++)
            {
                xmax[j] += TOL;
                xmin[j] -= TOL;
            }
            //%TODO::This might be modified after testing
            if (d_dim == 3)
            {
                for (int j = xmin[0] / mapdx[0]; j <= xmax[0] / mapdx[0]; j++)
                {
                    for (int k = xmin[1] / mapdx[1]; k <= xmax[1] / mapdx[1]; k++)
                    {
                        for (int l = xmin[2] / mapdx[2]; l <= xmax[2] / mapdx[2]; l++)
                        {
                            idx[0] = TIOGA_MAX(TIOGA_MIN(j, mapdims[0] - 1), 0);
                            idx[1] = TIOGA_MAX(TIOGA_MIN(k, mapdims[1] - 1), 0);
                            idx[2] = TIOGA_MAX(TIOGA_MIN(l, mapdims[2] - 1), 0);
                            mapmask[idx[2] * mapdims[1] * mapdims[0] + idx[1] * mapdims[0] + idx[0]] = 1;
                        }
                    }
                }
            }
            else if(d_dim==2)
            {
                for (int j = xmin[0] / mapdx[0]; j <= xmax[0] / mapdx[0]; j++)
                {
                    for (int k = xmin[1] / mapdx[1]; k <= xmax[1] / mapdx[1]; k++)
                    {
                        idx[0] = TIOGA_MAX(TIOGA_MIN(j, mapdims[0] - 1), 0);
                        idx[1] = TIOGA_MAX(TIOGA_MIN(k, mapdims[1] - 1), 0);
                        mapmask[ idx[1] * mapdims[0] + idx[0]] = 1;
                    }
                }

            }
            else
            {
                throw std::runtime_error("Invalid dimension value\n");
            }
            if (itag)
            {
                for (int m = 0; m < nvert; m++)
                {
                    nodeRes[inode[m]] = BIGVALUE;//Set the outter layer of cells to UNACCEPTABLE DONORS
                    iextmp[inode[m]] = iextmp1[inode[m]] = 1;//Outtest cells' node's itxtmpl all set to 1
                }
            }
        }
    }
    //Pushing forward from cells with mandatory receptors as nodes and make mexclude layers of unacceptable donor cells
    for (int iex = 0; iex < mexclude; iex++)
    {
        int k = 0;
        for (int n = 0; n < ntypes; n++)
        {
            int nvert = nv[n];
            for (int i = 0; i < nc[n]; i++)
            {
                for (int m = 0; m < nvert; m++)
                {
                    inode[m] = vconn[n][nvert * i + m] - BASE;
                    if (iextmp[inode[m]] == 1)
                    {
                        cellRes[k] = BIGVALUE;
                        break;
                    }
                }
                if (cellRes[k] == BIGVALUE)
                {
                    for (int m = 0; m < nvert; m++)
                    {
                        inode[m] = vconn[n][nvert * i + m] - BASE;
                        if (iextmp[inode[m]] != 1) iextmp1[inode[m]] = 1;//%TODO::I still think this if statement is unnecessary
                    }
                }
                k++;
            }
        }
        for (int i = 0; i < nnodes; i++)
        {
            iextmp[i] = iextmp1[i];
        }
    }
    TIOGA_FREE(iflag);
    TIOGA_FREE(iextmp);
    TIOGA_FREE(iextmp1);
}


void MeshBlock::writeGridFile(const std::string& filename)//This is 2d only for test
{
    std::string localFilename = filename + std::to_string(meshtag) + std::to_string(myid) + ".dat";
    std::ofstream fout;

    fout.open(localFilename);
    if (!fout.is_open())
    {
        throw std::runtime_error("WriteGridFile failure\n");
    }
    //Assume that the vconn is organized in a Tecplot friendly style
    fout << "TITLE =\"Tioga output\"\n";
    fout << "VARIABLES=\"X\",\"Y\",";
    if (d_dim == 3)
    {
        fout << "\"Z\",";
    }
    fout << "\"IBLANK\"\n";
    fout << "ZONE T=\"VOL_MIXED\",N=" << nnodes << " E=" << ncells;
    if (d_dim == 2)
    {
        fout<< " ET = QUADRILATERAL, F = FEPOINT\n";
    }
    else if (d_dim == 3)
    {
        fout << " ET = BRICK, F = FEPOINT\n";
    }
    for (int i = 0; i < nnodes; i++)
    {
        for (int j = 0; j < d_dim; j++)
        {
            fout << x[d_dim * i + j] << " ";
        }
        fout << iblank[i] << '\n';
    }
    int ba = 1 - BASE;
    if (d_dim == 2)
    {
        for (int n = 0; n < ntypes; n++)
        {
            int nvert = nv[n];
            switch (nvert)
            {
            case 3:
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 2] + ba << "\n";
                }
                break;
            case 4:
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 3] + ba << "\n";
                }
                break;
            default:
                break;
            }
        }
    }
    else if (d_dim == 3)
    {
        for (int n = 0; n < ntypes; n++)
        {
            int nvert = nv[n];
            switch (nvert)
            {
            case 4://Tetrahedron [1233,4444]
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 2] + ba << " ";
                    fout << vconn[n][nvert * i + 3] + ba << " " << vconn[n][nvert * i + 3] + ba << " " << vconn[n][nvert * i + 3] + ba << " " << vconn[n][nvert * i + 3] + ba << "\n";
                }
                break;
            case 5://[1234 5555]
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 3] + ba << " ";
                    fout << vconn[n][nvert * i + 4] + ba << " " << vconn[n][nvert * i + 4] + ba << " " << vconn[n][nvert * i + 4] + ba << " " << vconn[n][nvert * i + 4] + ba << "\n";
                }
                break;
            case 6://[1233 4566]
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 2] + ba << " ";
                    fout << vconn[n][nvert * i + 3] + ba << " " << vconn[n][nvert * i + 4] + ba << " " << vconn[n][nvert * i + 5] + ba << " " << vconn[n][nvert * i + 5] + ba << "\n";
                }
                break;
            case 8:
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 3] + ba << " ";
                    fout << vconn[n][nvert * i + 4] + ba << " " << vconn[n][nvert * i + 5] + ba << " " << vconn[n][nvert * i + 6] + ba << " " << vconn[n][nvert * i + 7] + ba << "\n";
                }
                break;
            default:
                break;
            }
        }
    }
    fout.close();
    return;

}

void MeshBlock::writeCellFile(const std::string& filename)
{
    std::string localFilename = filename + std::to_string(meshtag) + std::to_string(myid) + ".dat";
    std::ofstream fout;

    fout.open(localFilename);
    if (!fout.is_open())
    {
        throw std::runtime_error("WriteCellFile failure\n");
    }
    //Assume that the vconn is organized in a Tecplot friendly style
    fout << "TITLE =\"Tioga output\"\n";
    fout << "VARIABLES=\"X\",\"Y\",";
    if (d_dim == 3)
    {
        fout << "\"Z\",";
    }
    fout << "\"IBLANK\",\"IBLANK_CELL\"\n";
    fout << "ZONE T=\"VOL_MIXED\",N=" << nnodes << " E =" << ncells;
    if (d_dim == 2)
    {
        fout << " ET = QUADRILATERAL, F = FEBLOCK\n";
    }
    else if (d_dim == 3)
    {
        fout << " ET = BRICK, F = FEBLOCK\n";
    } ;
    fout << "VARLOCATION =  (";
    for (int i = 0; i < d_dim; i++)
    {
        fout << i + 1 << "=NODAL,";
    }
    fout << d_dim + 1 << "=NODAL," << d_dim + 2 << "=CELLCENTERED)\n";
    for (int j = 0; j < d_dim; j++)
    {
        for (int i = 0; i < nnodes; i++)
        {
            fout << x[d_dim * i + j]<<'\n';
        }
    }
    for (int i = 0; i < nnodes; i++)
    {
        fout << iblank[i] << '\n';
    }
    for (int i = 0; i < ncells; i++)
    {
        fout << ' ' << iblank_cell[i] << ' ' << '\n';
    }
    int ba = 1 - BASE;
    if (d_dim == 2)
    {
        for (int n = 0; n < ntypes; n++)
        {
            int nvert = nv[n];
            switch (nvert)
            {
            case 3:
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 2] + ba << "\n";
                }
                break;
            case 4:
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 3] + ba << "\n";
                }
                break;
            default:
                break;
            }
        }
    }
    else if (d_dim == 3)
    {
        for (int n = 0; n < ntypes; n++)
        {
            int nvert = nv[n];
            switch (nvert)
            {
            case 4://Tetrahedron [1233,4444]
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 2] + ba << " ";
                    fout << vconn[n][nvert * i + 3] + ba << " " << vconn[n][nvert * i + 3] + ba << " " << vconn[n][nvert * i + 3] + ba << " " << vconn[n][nvert * i + 3] + ba << "\n";
                }
                break;
            case 5://[1234 5555]
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 3] + ba << " ";
                    fout << vconn[n][nvert * i + 4] + ba << " " << vconn[n][nvert * i + 4] + ba << " " << vconn[n][nvert * i + 4] + ba << " " << vconn[n][nvert * i + 4] + ba << "\n";
                }
                break;
            case 6://[1233 4566]
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 2] + ba << " ";
                    fout << vconn[n][nvert * i + 3] + ba << " " << vconn[n][nvert * i + 4] + ba << " " << vconn[n][nvert * i + 5] + ba << " " << vconn[n][nvert * i + 5] + ba << "\n";
                }
                break;
            case 8:
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 3] + ba << " ";
                    fout << vconn[n][nvert * i + 4] + ba << " " << vconn[n][nvert * i + 5] + ba << " " << vconn[n][nvert * i + 6] + ba << " " << vconn[n][nvert * i + 7] + ba << "\n";
                }
                break;
            default:
                break;
            }
        }
    }
    fout.close();
    return;
}

void MeshBlock::getWallBounds(int* mtag, int* existWall, double wbox[6])
{
    int id;
    int inode;

    *mtag = meshtag + (1 - BASE);
    if (nwbc <= 0) 
    {
        *existWall = 0;
        for (int i = 0; i < 2 * d_dim; i++) wbox[i] = 0;
        return;
    }

    *existWall = 1;
    for (int i = 0; i < d_dim; i++)
    {
        wbox[i] = BIGVALUE;
        wbox[d_dim + i] = -BIGVALUE;
    }


    for (int i = 0; i < nwbc; i++)
    {
        inode = wbcnode[i] - BASE;
        id = d_dim * inode;
        for (int j = 0; j < d_dim; j++)
        {
            wbox[j] = TIOGA_MIN(wbox[j], x[id + j]);
            wbox[d_dim + j] = TIOGA_MAX(wbox[j + d_dim], x[id + j]);
        }
    }

}

void MeshBlock::markWallBoundary(int* sam, int nx[3], double extents[6])
{
    int nvert;
    int ii, jj, kk, mm;
    int id, iv;
    int* iflag = new int[ncells];
    int* inode = new int[nnodes];
   
    double ds[3];
    double xv;
    int imin[3];
    int imax[3];

    //
    for (int i = 0; i < ncells; i++) iflag[i] = 0;
    for (int i = 0; i < nnodes; i++) inode[i] = 0;
    //
    for (int i = 0; i < nwbc; i++)
    {
        ii = wbcnode[i] - BASE;
        inode[ii] = 1;
    }

    // mark wall boundary cells
    //
    int m = 0;
    for (int n = 0; n < ntypes; n++)
    {
        nvert = nv[n];
        for (int i = 0; i < nc[n]; i++)
        {
            for (int j = 0; j < nvert; j++)
            {
                ii = vconn[n][nvert * i + j] - BASE;
                if (inode[ii] == 1)
                {
                    iflag[m] = 1;
                    break;
                }
            }
            m++;
        }
    }
    //
    // find delta's in each directions
    //
    for (int k = 0; k < d_dim; k++) ds[k] = (extents[k + d_dim] - extents[k]) / nx[k];
    //
    // mark sam cells with wall boundary cells now
    //
    m = 0;
    for (int n = 0; n < ntypes; n++)
    {
        nvert = nv[n];
        for (int i = 0; i < nc[n]; i++)
        {
            if (iflag[m] == 1)
            {
                //
                // find the index bounds of each wall boundary cell
                // bounding box
                //
                for (int j = 0; j < d_dim; j++)
                {
                    imin[j] = BIGINT;
                    imax[j] = -BIGINT;
                }
                for (int j = 0; j < nvert; j++)
                {
                    id = d_dim * (vconn[n][nvert * i + j] - BASE);
                    for (int k = 0; k < d_dim; k++)
                    {
                        xv = x[id + k];
                        iv = floor((xv - extents[k]) / ds[k]);
                        imin[k] = TIOGA_MIN(imin[k], iv);
                        imax[k] = TIOGA_MAX(imax[k], iv);
                    }
                }
                for (int j = 0; j < d_dim; j++)
                {
                    imin[j] = TIOGA_MAX(imin[j], 0);
                    imax[j] = TIOGA_MIN(imax[j], nx[j] - 1);
                }
                //
                // mark sam to 1
                //
                if (d_dim == 3)
                {
                    for (kk = imin[2]; kk < imax[2] + 1; kk++)
                        for (jj = imin[1]; jj < imax[1] + 1; jj++)
                            for (ii = imin[0]; ii < imax[0] + 1; ii++)
                            {
                                mm = kk * nx[1] * nx[0] + jj * nx[0] + ii;
                                sam[mm] = 2;
                            }
                }
                else if (d_dim == 2)
                {
                    for (jj = imin[1]; jj < imax[1] + 1; jj++)
                        for (ii = imin[0]; ii < imax[0] + 1; ii++)
                        {
                            mm =  jj * nx[0] + ii;
                            sam[mm] = 2;
                        }
                }
            }
            m++;
        }
    }
    TIOGA_FREE(iflag);
    TIOGA_FREE(inode);
}

void MeshBlock::getReducedOBB(OBB* obc, double* realData)
{
    /*std::cout << "obc's xc :" << obc->xc[0] << " " << obc->xc[1] << " " << obc->xc[2] << '\n';
    std::cout << "obc's dxc :" << obc->dxc[0] << " " << obc->dxc[1] << " " << obc->dxc[2] << '\n';
    std::cout << "obc's vec :" << obc->vec[0][0] << " " << obc->vec[0][1] << " " << obc->vec[0][2] << '\n';
    std::cout << "obc's vec :" << obc->vec[1][0] << " " << obc->vec[1][1] << " " << obc->vec[1][2] << '\n';
    std::cout << "obc's vec :" << obc->vec[2][0] << " " << obc->vec[2][1] << " " << obc->vec[2][2] << '\n';*/

    double bbox[6], xd[3];
    bool iflag;
    for (int i = 0; i < d_dim; i++)
    {
        realData[i] = BIGVALUE;
        realData[d_dim + i] = -BIGVALUE;
    }
    int count_of_cell_involved = 0;
    int id;
    for (int n = 0; n < ntypes; n++)
    {
        int nvert = nv[n];
        for (int i = 0; i < nc[n]; i++)
        {
            for (int j = 0; j < d_dim; j++)
            {
                bbox[j] = BIGVALUE;
                bbox[j + d_dim] = -BIGVALUE;
            }
            for (int m = 0; m < nvert; m++)
            {
                id = d_dim * (vconn[n][nvert * i + m] - BASE);
                for (int j = 0; j < d_dim; j++)
                {
                    xd[j] = 0;
                }
                for (int j = 0; j < d_dim; j++)
                {
                    for (int k = 0; k < d_dim; k++)
                    {
                        xd[j] += (x[id + k] - obc->xc[k]) * obc->vec[j][k];
                    }
                }
                for (int j = 0; j < d_dim; j++)
                {
                    bbox[j] = TIOGA_MIN(bbox[j], xd[j]);
                    bbox[d_dim + j] = TIOGA_MAX(bbox[d_dim + j], xd[j]);
                }
            }
            iflag = false;
            for (int j = 0; j < d_dim; j++)
            {
                iflag = iflag or (bbox[j] > obc->dxc[j]);//if the bbox's minimal value is bigger than the maximum
            }
            if (iflag) continue;
            iflag = false;
            
            for (int j = 0; j < d_dim; j++)
            {
                iflag = iflag or (bbox[d_dim + j] < -obc->dxc[j]);
            }
            if (iflag) continue;
            count_of_cell_involved++;
            for (int m = 0; m < nvert; m++)
            {
                id = d_dim * (vconn[n][nvert * i + m] - BASE);
                for (int j = 0; j < d_dim; j++)
                {
                    xd[j] = 0;
                }
                for (int j = 0; j < d_dim; j++)
                {
                    for (int k = 0; k < d_dim; k++)
                    {
                        xd[j] += (x[id + k] - obb->xc[k]) * obb->vec[j][k];
                    }
                }
                for (int j = 0; j < d_dim; j++)
                {
                    realData[j] = TIOGA_MIN(realData[j], xd[j]);
                    realData[d_dim + j] = TIOGA_MAX(realData[d_dim + j], xd[j]);
                }
                for (int j = 0; j < d_dim; j++)
                {
                    realData[j] = TIOGA_MIN(realData[j], xd[j]);
                    realData[d_dim + j] = TIOGA_MAX(realData[d_dim + j], xd[j]);
                }
            }
        }
    }
   // std::cout << "reduced data raw is "<<count_of_cell_involved<<" cell's result:" << realData[0] << " " << realData[1] << " " << realData[2] << '\n';
    //std::cout << realData[3] << " " << realData[4] << " " << realData[5] << '\n';
    for (int i = 0; i < 2 * d_dim; i++)
    {
        bbox[i] = realData[i];
    }
    for (int i = 0; i < d_dim; i++)
    {
        realData[i] = obb->xc[i];
        for (int j = 0; j < d_dim; j++)
        {
            realData[i] += ((bbox[j] + bbox[d_dim + j]) * 0.5) * obb->vec[j][i];
        }
        realData[d_dim + i] = (bbox[d_dim + i] - bbox[i]) * 0.51;
    }
    if (count_of_cell_involved == 0)
    {
        for (int i = 0; i < 2 * d_dim; i++)
        {
            realData[i] = 0;
        }
    }
    std::cout << "reduced data on meshtag" <<meshtag<<" proc "<<myid<<" is " << count_of_cell_involved << " cell's result:" << realData[0] << " " << realData[1] << " " << realData[2] << '\n';
    std::cout << realData[3] << " " << realData[4] << " " << realData[5] << '\n';
    return;
}

void MeshBlock::getQueryPoints(OBB* obc,
    int* nints, int** intData,
    int* nreals, double** realData)
{
    int id;
    double xd[3];
    bool isIn = true;
    int* inode = new int[nnodes];
    *nints = *nreals = 0;
    for (int i = 0; i < nnodes; i++)
    {
        id = d_dim * i;
        for (int j = 0; j < d_dim; j++) xd[j] = 0;
        for (int j = 0; j < d_dim; j++)
            for (int k = 0; k < d_dim; k++)
                xd[j] += (x[id + k] - obc->xc[k]) * obc->vec[j][k];
        isIn = true;
        for (int j = 0; j < d_dim; j++)
        {
            isIn = isIn and fabs(xd[j]) <= obc->dxc[j];
        }
        if (isIn)
        {
            inode[*nints] = i;
            (*nints)++;
            (*nreals) += d_dim+1;

        }
    }
    (*intData) = new int[*nints];
    (*realData) = new double[*nreals];

    for (int i = 0; i < *nints; i++)
    {
        id = d_dim * inode[i];
        (*intData)[i] = inode[i];
        for (int j = 0; j < d_dim; j++)
        {
            (*realData)[(d_dim+1) * i + j] = x[id + j];
        }
        (*realData)[(d_dim + 1) * i + d_dim] = nodeRes[inode[i]];
    }
    //
    TIOGA_FREE(inode);
}

void MeshBlock::getQueryPoints2(OBB* obc,
    int* nints, int** intData,
    int* nreals, double** realData)
{
    double xv[8][3];
    double xd[3];
    double xc[3];
    int imin[3], imax[3];
    double mdx[3];
    double xmin[3], xmax[3];
    bool iflag;
    int* inode = new int[nnodes];
    int* signOfTheTime = new int[nnodes];
    for (int i = 0; i < nnodes; i++)
    {
        signOfTheTime[i] = 1;
    }
    *nints = *nreals = 0;
    //First decide the range of obb subblocks to be checked with obc's bbox in obb coordinates
    getobbcoords(d_dim, obc->xc, obc->dxc, obc->vec, xv);//xv of 2D obc is stored at the first 4 
    for (int i = 0; i < d_dim; i++)
    {
        xmin[i] = BIGVALUE;
        xmax[i] = -BIGVALUE;
    }
    for (int n = 0; n < pow(2,d_dim); n++)
    {
        transform2OBB(d_dim, xv[n], obb->xc, obb->vec, xd);//Transfer obc's nodes into obb coordinates
        for (int j = 0; j < d_dim; j++)
        {
            xmin[j] = TIOGA_MIN(xmin[j], xd[j] + obb->dxc[j]);//Ranging from 0 to 2*dxc
            xmax[j] = TIOGA_MAX(xmax[j], xd[j] + obb->dxc[j]);
        }
    }
    for (int i = 0; i < d_dim; i++)
    {
        double delta = 0.01 * (xmax[i] - xmin[i]);
        //std::cout << "xmax " << xmax[i] << " xmin " << xmin[i] << " mapdx " << mapdx[i] << " imax " << imax[i] << " imin " << imin[i] << '\n';
        xmin[i] -= delta;
        xmax[i] += delta;

        imin[i] = TIOGA_MAX(static_cast<int>(xmin[i] / mapdx[i]), 0);
        imax[i] = TIOGA_MIN(static_cast<int>(xmax[i] / mapdx[i]), mapdims[i] - 1);
        mdx[i] = 0.5 * mapdx[i];
    }
    getobbcoords(d_dim,obc->xc, mdx, obb->vec, xv);//Transfer a single sub-block( in obb)located at obc center  into  xv
    for (int i = 0; i < d_dim; i++)
    {
        xmin[i] = BIGVALUE;
        xmax[i] = -BIGVALUE;
    }
    for (int n = 0; n < pow(2, d_dim); n++)
    {
        transform2OBB(d_dim, xv[n], obc->xc, obc->vec, xd);//Get the min_max bouding box for a ll
        for (int i = 0; i < d_dim; i++)
        {
            xmin[i] = TIOGA_MIN(xmin[i], xd[i]);
            xmax[i] = TIOGA_MAX(xmax[i], xd[i]);
        }
    }
    //std::cout << "imin: " << imin[0] << " " << imin[1] << '\n';
    //std::cout << "imax: " << imax[0] << " " << imax[1] << '\n';
    //Count how many query points there are and store in inode
    if (d_dim == 2)
    {
        //std::cout << "Entering d_dim==2 getQueryPoint2\n";
        for (int k = imin[1]; k < imax[1] + 1; k++)
        {
            for (int j = imin[0]; j < imax[0] + 1; j++)
            {
                xd[0] = -obb->dxc[0] + j * mapdx[0] + mapdx[0] * 0.5;
                xd[1] = -obb->dxc[1] + k * mapdx[1] + mapdx[1] * 0.5;
                for (int n = 0; n < 2; n++)
                {
                    xc[n] = obb->xc[n];
                    for (int ij = 0; ij < 2; ij++)
                        xc[n] += (xd[ij] * obb->vec[ij][n]);//Transfer to cartesian coordinates
                }
                transform2OBB(d_dim, xc, obc->xc, obc->vec, xd);//Transfer the center of the subblock to the obc coordinates, then check if there is overlap
                // check if this sub-block overlaps OBC
                //
                iflag = false;
                for (int ij = 0; ij < 2 && !iflag; ij++)
                    iflag = (iflag or (xmin[ij] + xd[ij] > obc->dxc[ij]));
                if (iflag) continue;
                iflag = 0;
                for (int ij = 0; ij < 2 && !iflag; ij++)
                    iflag = (iflag or (xmax[ij] + xd[ij] < -obc->dxc[ij]));
                if (iflag) continue;
                //std::cout << "There is overlap!!!!!!!!!!!!!!!\n";
                // if there overlap
                // go through points within the sub-block
                // to figure out what needs to be send
                //
                int indx = k * mapdims[0] + j;
                for (int m = icft[indx]; m < icft[indx + 1]; m++)
                {
                    int i2 = 2 * invmap[m];
                    for (int ik = 0; ik < 2; ik++)
                        xc[ik] = x[i2 + ik];
                    transform2OBB(d_dim, xc, obc->xc, obc->vec, xd);
                    if (fabs(xd[0]) <= obc->dxc[0] &&
                        fabs(xd[1]) <= obc->dxc[1])
                    {
                        inode[*nints] = invmap[m];
                        (*nints)++;
                        (*nreals) += 3;
                    }
                }

            }
        }
    }
    else if (d_dim == 3)
    {
        for (int l = imin[2]; l < imax[2] + 1; l++)
        {
            for (int k = imin[1]; k < imax[1] + 1; k++)
            {
                for (int j = imin[0]; j < imax[0] + 1; j++)
                {
                    xd[0] = -obb->dxc[0] + j * mapdx[0] + mapdx[0] * 0.5;
                    xd[1] = -obb->dxc[1] + k * mapdx[1] + mapdx[1] * 0.5;
                    xd[2] = -obb->dxc[2] + l * mapdx[2] + mapdx[2] * 0.5;
                    for (int n = 0; n < 3; n++)
                    {
                        xc[n] = obb->xc[n];
                        for (int ij = 0; ij < 3; ij++)
                            xc[n] += (xd[ij] * obb->vec[ij][n]);
                    }

                    transform2OBB(d_dim, xc, obc->xc, obc->vec, xd);
                    // check if this sub-block overlaps OBC
                    //
                    iflag = false;
                    for (int ij = 0; ij < 3 && !iflag; ij++)
                        iflag = (iflag or (xmin[ij] + xd[ij] > obc->dxc[ij]));
                    if (iflag) continue;
                    iflag = 0;
                    for (int ij = 0; ij < 3 && !iflag; ij++)
                        iflag = (iflag or (xmax[ij] + xd[ij] < -obc->dxc[ij]));
                    if (iflag) continue;
                    
                // if there overlap
                // go through points within the sub-block
                // to figure out what needs to be send
                //
                    int indx = l * mapdims[1] * mapdims[0] + k * mapdims[0] + j;
                    for (int m = icft[indx]; m < icft[indx + 1]; m++)
                    {
                        int i3 = 3 * invmap[m];
                        for (int ik = 0; ik < 3; ik++) 
                            xc[ik] = x[i3 + ik];
                        transform2OBB(d_dim, xc, obc->xc, obc->vec, xd);
                        if (fabs(xd[0]) <= obc->dxc[0] &&
                            fabs(xd[1]) <= obc->dxc[1] &&
                            fabs(xd[2]) <= obc->dxc[2])
                        {
                            inode[*nints] = invmap[m];
                            (*nints)++;
                            (*nreals) += 4;
                        }
                    }

                }

            }
        }
    }
#ifdef TIOGA_HAS_NODEGID
    int nintsPerNode = 3;
#else
    int nintsPerNode = 1;
#endif
    (*intData) = new int[(*nints) * nintsPerNode];
    (*realData) = new double[(*nreals)];
    //
    int m = 0;
    int iidx = 0;
    int id;
    for (int i = 0; i < *nints; i++) 
    {
        signOfTheTime[inode[i]] = -1;
        id = d_dim * inode[i];//inode starting from 0;
        (*intData)[iidx++] = inode[i];
#ifdef TIOGA_HAS_NODEGID
        std::memcpy(&(*intData)[iidx], &nodeGID[inode[i]], sizeof(uint64_t));
        iidx += 2;
#endif
        for (int j = 0; j < d_dim; j++)
        {
            (*realData)[m++] = x[id+j];
        }
        (*realData)[m++] = nodeRes[inode[i]];
    }
    // Adjust nints to the proper array size
    *nints *= nintsPerNode;
    writeGridFile2("searchResult", signOfTheTime);
    TIOGA_FREE(inode);
    TIOGA_FREE(signOfTheTime)
}

void MeshBlock::getReducedOBB2(OBB* obc, double* realData)
{
//Not implemented yet

}

void MeshBlock::writeOBB(const std::string& filename)
{
    std::string localFilename = filename + std::to_string(meshtag) + std::to_string(myid) + ".dat";
    std::ofstream fout;

    fout.open(localFilename);
    if (!fout.is_open())
    {
        throw std::runtime_error("writeOBB failure\n");
    }
    //Assume that the vconn is organized in a Tecplot friendly style
    fout << "TITLE =\"Box file\"\n";
    if (d_dim == 2)
    {
        fout << "VARIABLES=\"X\",\"Y\"\n";
        fout << "ZONE T=\"VOL_MIXED\",N=4 E=1 ET=QUADRILATERAL, F=FEPOINT\n";
        double xx[2];
        for (int l = 0; l < 2; l++)
        {
            int il = 2 * (l % 2) - 1;
            for (int k = 0; k < 2; k++)
            {
                int ik = 2 * (k % 2) - 1;
                xx[0] = xx[1] = 0;
                for (int m = 0; m < 2; m++)
                {
                    xx[m] = obb->xc[m] + il * obb->vec[0][m] * obb->dxc[0] + ik * obb->vec[1][m] * obb->dxc[1];
                }

                fout << xx[0] << " " << xx[1] << '\n';
            }
        }
        fout << "1 2 4 3\n";
        fout << obb->xc[0] << " " << obb->xc[1] << '\n';
        fout << obb->dxc[0] << " " << obb->dxc[1] << '\n';
    }
    else if (d_dim == 3)
    {
        fout << "VARIABLES=\"X\",\"Y\",\"Z\"\n";
        fout << "ZONE T=\"VOL_MIXED\",N=8 E=1 ET=BRICK, F=FEPOINT\n";
        double xx[3];
        for (int l = 0; l < 2; l++)
        {
            int il = 2 * (l % 2) - 1;
            for (int k = 0; k < 2; k++)
            {
                int ik = 2 * (k % 2) - 1;
                for (int j = 0; j < 2; j++)
                {
                    int ij = 2 * (j % 2) - 1;
                    xx[0] = xx[1] = xx[2] = 0;
                    for (int m = 0; m < 3; m++)
                    {
                        xx[m] = obb->xc[m] + il * obb->vec[0][m] * obb->dxc[0] + ik * obb->vec[1][m] * obb->dxc[1] + ij * obb->vec[2][m] * obb->dxc[2];
                    }
                    fout << xx[0] << " " << xx[1] << " " << xx[2] << '\n';
                }
            }
        }
        fout << "1 2 4 3 5 6 8 7\n";

        fout << obb->xc[0] << " " << obb->xc[1] << " " << obb->xc[2] << '\n';
        fout << obb->dxc[0] << " " << obb->dxc[1] << " " << obb->dxc[2] << '\n';
    }
    fout.close();
}

void MeshBlock::writeOBB2(const std::string& filename,OBB* obc)
{
    std::string localFilename = filename + std::to_string(meshtag) + std::to_string(myid) + ".dat";
    std::ofstream fout;
    fout.open(localFilename);
    if (!fout.is_open())
    {
        throw std::runtime_error("writeOBB failure\n");
    }
    //Assume that the vconn is organized in a Tecplot friendly style
    fout << "TITLE =\"Box file\"\n";
    if (d_dim == 2)
    {
        fout << "VARIABLES=\"X\",\"Y\"\n";
        fout << "ZONE T=\"VOL_MIXED\",N=4 E=1 ET=QUADRILATERAL, F=FEPOINT\n";
        double xx[2];
        for (int l = 0; l < 2; l++)
        {
            int il = 2 * (l % 2) - 1;
            for (int k = 0; k < 2; k++)
            {
                int ik = 2 * (k % 2) - 1;
                xx[0] = xx[1] = 0;
                for (int m = 0; m < 2; m++)
                {
                    xx[m] = obc->xc[m] + il * obc->vec[0][m] * obc->dxc[0] + ik * obc->vec[1][m] * obc->dxc[1];
                }

                fout << xx[0] << " " << xx[1] << '\n';
            }
        }
        fout << "1 2 4 3\n";
        fout << obc->xc[0] << " " << obc->xc[1] << '\n';
        fout << obc->dxc[0] << " " << obc->dxc[1] << '\n';
    }
    else if (d_dim == 3)
    {
        fout << "VARIABLES=\"X\",\"Y\",\"Z\"\n";
        fout << "ZONE T=\"VOL_MIXED\",N=8 E=1 ET=BRICK, F=FEPOINT\n";
        double xx[3];
        for (int l = 0; l < 2; l++)
        {
            int il = 2 * (l % 2) - 1;
            for (int k = 0; k < 2; k++)
            {
                int ik = 2 * (k % 2) - 1;
                for (int j = 0; j < 2; j++)
                {
                    int ij = 2 * (j % 2) - 1;
                    xx[0] = xx[1] = xx[2] = 0;
                    for (int m = 0; m < 3; m++)
                    {
                        xx[m] = obc->xc[m] + il * obc->vec[0][m] * obc->dxc[0] + ik * obc->vec[1][m] * obc->dxc[1] + ij * obc->vec[2][m] * obc->dxc[2];
                    }
                    fout << xx[0] << " " << xx[1] << " " << xx[2] << '\n';
                }
            }
        }
        fout << "1 2 4 3 5 6 8 7\n";
        
        fout << obc->xc[0] << " " << obc->xc[1] << " " << obc->xc[2] << '\n';
        fout << obc->dxc[0] << " " << obc->dxc[1] << " " << obc->dxc[2] << '\n';
    }
    fout.close();
}

void MeshBlock::writeOBBPair(const std::string& filename, OBB* obc)
{
    std::string localFilename = filename + std::to_string(meshtag) + std::to_string(myid) + ".dat";
    std::ofstream fout;
    fout.open(localFilename);
    if (!fout.is_open())
    {
        throw std::runtime_error("writeOBB failure\n");
    }
    //Assume that the vconn is organized in a Tecplot friendly style
    fout << "TITLE =\"Box file\"\n";
    if (d_dim == 2)
    {
        fout << "VARIABLES=\"X\",\"Y\"\n";
        fout << "ZONE T=\"VOL_MIXED\",N=8 E=2 ET=QUADRILATERAL, F=FEPOINT\n";
        double xx[2];
        for (int l = 0; l < 2; l++)
        {
            int il = 2 * (l % 2) - 1;
            for (int k = 0; k < 2; k++)
            {
                int ik = 2 * (k % 2) - 1;
                xx[0] = xx[1] = 0;
                for (int m = 0; m < 2; m++)
                {
                    xx[m] = obb->xc[m] + il * obb->vec[0][m] * obb->dxc[0] + ik * obb->vec[1][m] * obb->dxc[1];
                }

                fout << xx[0] << " " << xx[1] << '\n';
            }
        }
        for (int l = 0; l < 2; l++)
        {
            int il = 2 * (l % 2) - 1;
            for (int k = 0; k < 2; k++)
            {
                int ik = 2 * (k % 2) - 1;
                xx[0] = xx[1] = 0;
                for (int m = 0; m < 2; m++)
                {
                    xx[m] = obc->xc[m] + il * obc->vec[0][m] * obc->dxc[0] + ik * obc->vec[1][m] * obc->dxc[1];
                }

                fout << xx[0] << " " << xx[1] << '\n';
            }
        }
        fout << "1 2 4 3\n";
        fout << "5 6 8 7\n";
        fout << obc->xc[0] << " " << obc->xc[1] << '\n';
        fout << obc->dxc[0] << " " << obc->dxc[1] << '\n';
    }
    else if (d_dim == 3)
    {
        fout << "VARIABLES=\"X\",\"Y\",\"Z\"\n";
        fout << "ZONE T=\"VOL_MIXED\",N=16 E=2 ET=BRICK, F=FEPOINT\n";
        double xx[3];
        for (int l = 0; l < 2; l++)
        {
            int il = 2 * (l % 2) - 1;
            for (int k = 0; k < 2; k++)
            {
                int ik = 2 * (k % 2) - 1;
                for (int j = 0; j < 2; j++)
                {
                    int ij = 2 * (j % 2) - 1;
                    xx[0] = xx[1] = xx[2] = 0;
                    for (int m = 0; m < 3; m++)
                    {
                        xx[m] = obb->xc[m] + il * obb->vec[0][m] * obb->dxc[0] + ik * obb->vec[1][m] * obb->dxc[1] + ij * obb->vec[2][m] * obb->dxc[2];
                    }
                    fout << xx[0] << " " << xx[1] << " " << xx[2] << '\n';
                }
                

                
            }
        }
        for (int l = 0; l < 2; l++)
        {
            int il = 2 * (l % 2) - 1;
            for (int k = 0; k < 2; k++)
            {
                int ik = 2 * (k % 2) - 1;
                for (int j = 0; j < 2; j++)
                {
                    int ij = 2 * (j % 2) - 1;
                    xx[0] = xx[1] = xx[2] = 0;
                    for (int m = 0; m < 3; m++)
                    {
                        xx[m] = obc->xc[m] + il * obc->vec[0][m] * obc->dxc[0] + ik * obc->vec[1][m] * obc->dxc[1] + ij * obc->vec[2][m] * obc->dxc[2];
                    }
                    fout << xx[0] << " " << xx[1] << " " << xx[2] << '\n';
                }



            }
        }
        fout << "1 2 4 3 5 6 8 7\n";
        fout << "9 10 12 11 13 14 16 15\n";
        fout << obc->xc[0] << " " << obc->xc[1] << " " << obc->xc[2] << '\n';
        fout << obc->dxc[0] << " " << obc->dxc[1] << " " << obc->dxc[2] << '\n';

    }
    
    
    fout.close();
}
//Do it later

void MeshBlock::outPutSearch(const std::string& filename)
{
    if (iblank_cell) TIOGA_FREE(iblank_cell);
    if (iblank) TIOGA_FREE(iblank);
    iblank_cell = new int[ncells];
    iblank = new int[nnodes];

    for (int i = 0; i < nnodes; i++)
    {
        iblank[i] = -1;
    }
    for (int i = 0; i < ncells; i++)
    {
        iblank_cell[i] = 1;
    }
    for (int i = 0; i < nsearch; i++)
    {
        if (donorId[i] >= 0)
        {
            iblank_cell[donorId[i]] = -1;
        }
    }
    writeCellFile("searchCell");
    return;
    
}

void MeshBlock::writeMandatoryReceptor(const std::string& filename)
{
    std::string localFilename = filename + std::to_string(meshtag) + std::to_string(myid) + ".dat";
    std::ofstream fout;

    fout.open(localFilename);
    if (!fout.is_open())
    {
        throw std::runtime_error("WriteGridFile failure\n");
    }
    //Assume that the vconn is organized in a Tecplot friendly style
    fout << "TITLE =\"MandatoryReceptors\"\n";
    fout << "VARIABLES=\"X\",\"Y\",";
    if (d_dim == 3)
    {
        fout << "\"Z\",";
    }
    fout<<"\"IBLANK\"\n";
    fout << "ZONE T=\"VOL_MIXED\",N=" << nnodes << " E=" << ncells;
    if (d_dim == 2)
    {
        fout << " ET = QUADRILATERAL, F = FEPOINT\n";
    }
    else if (d_dim == 3)
    {
        fout << " ET = BRICK, F = FEPOINT\n";
    }
    
    for (int i = 0; i < nnodes; i++)
    {
        for (int j = 0; j < d_dim; j++)
        {
            fout << x[d_dim * i + j] << " ";
        }
        fout << (nodeRes[i] == BIGVALUE ? -1 : 1) << '\n';
    }
    int ba = 1 - BASE;
    if (d_dim == 2)
    {
        for (int n = 0; n < ntypes; n++)
        {
            int nvert = nv[n];
            switch (nvert)
            {
            case 3:
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 2] + ba << "\n";
                }
                break;
            case 4:
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 3] + ba << "\n";
                }
                break;
            default:
                break;
            }
        }
    }
    else if (d_dim == 3)
    {
        for (int n = 0; n < ntypes; n++)
        {
            int nvert = nv[n];
            switch (nvert)
            {
            case 4://Tetrahedron [1233,4444]
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 2] + ba << " ";
                    fout << vconn[n][nvert * i + 3] + ba << " " << vconn[n][nvert * i + 3] + ba << " " << vconn[n][nvert * i + 3] + ba << " " << vconn[n][nvert * i + 3] + ba << "\n";
                }
                break;
            case 5://[1234 5555]
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 3] + ba << " ";
                    fout << vconn[n][nvert * i + 4] + ba << " " << vconn[n][nvert * i + 4] + ba << " " << vconn[n][nvert * i + 4] + ba << " " << vconn[n][nvert * i + 4] + ba << "\n";
                }
                break;
            case 6://[1233 4566]
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 2] + ba << " ";
                    fout << vconn[n][nvert * i + 3] + ba << " " << vconn[n][nvert * i + 4] + ba << " " << vconn[n][nvert * i + 5] + ba << " " << vconn[n][nvert * i + 5] + ba << "\n";
                }
                break;
            case 8:
                for (int i = 0; i < nc[n]; i++)
                {
                    fout << vconn[n][nvert * i] + ba << " " << vconn[n][nvert * i + 1] + ba << " " << vconn[n][nvert * i + 2] + ba << " " << vconn[n][nvert * i + 3] + ba << " ";
                    fout << vconn[n][nvert * i + 4] + ba << " " << vconn[n][nvert * i + 5] + ba << " " << vconn[n][nvert * i + 6] + ba << " " << vconn[n][nvert * i + 7] + ba << "\n";
                }
                break;
            default:
                break;
            }
        }
    }
    fout.close();
    
    return;

}

void MeshBlock::writeObc(const std::string& filename)
{
    std::string localFilename = filename + std::to_string(meshtag) + std::to_string(myid) + ".sct";
    std::ofstream fout;

    fout.open(localFilename);
    if (!fout.is_open())
    {
        throw std::runtime_error("WriteGridFile failure\n");
    }
    for (int i = 0; i < nobc; i++)
    {
        int inode = obcnode[i] - 1;
        for (int j = 0; j < d_dim; j++)
        {
            fout << x[d_dim * inode + j] << '\t';
        }
        fout << '\n';
    }
    fout.close();
}
void MeshBlock::writeWbc(const std::string& filename)
{
    std::string localFilename = filename + std::to_string(meshtag) + std::to_string(myid) + ".sct";
    std::ofstream fout;

    fout.open(localFilename);
    if (!fout.is_open())
    {
        throw std::runtime_error("WriteGridFile failure\n");
    }
    for (int i = 0; i < nwbc; i++)
    {
        int inode = wbcnode[i] - 1;
        for (int j = 0; j < d_dim; j++)
        {
            fout << x[d_dim * inode + j] << '\t';
        }
        fout << '\n';
    }
    fout.close();
}

void writeFlowFile(int bid, double* q, int nvar, int type);
void setResolutions(double* nres, double* cres);
void search();
void search_uniform_hex();


void writeOBB2(OBB* obc, int bid);

#include "hdf5.h"

void MeshBlock::writeGridHDF(const std::string& filename)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id,comm,info);
    hid_t file_id = H5Fcreate(filename.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,plist_id);
    
    std::string curGroupName = "Mesh"+std::to_string(meshtag);
    hid_t MeshGroup_id = H5Gcreate(file_id,curGroupName.c_str(),H5P_DEFAULT,H5P_DEFAULT,plist_id);
    //return;
    hid_t ProcMeshGroup_id = H5Gcreate(MeshGroup_id,("Proc"+std::to_string(myid)).c_str(),H5P_DEFAULT,H5P_DEFAULT,plist_id);
    hsize_t dims[2] = {nnodes,d_dim};
    //return;
    hid_t dataspace_id = H5Screate_simple(2,dims,NULL);
    //return;
    hid_t dataset_id = H5Dcreate(ProcMeshGroup_id,"Coordinates",H5T_NATIVE_DOUBLE,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,plist_id);
    
    hid_t status = H5Dwrite(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,rxyz);
    
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
    //return;

    
    int buffer[8];
    if(d_dim==2)
    {
        hsize_t count[2] = {1,4};
        hsize_t offset[2] = {0,0};
        hsize_t stride[2] = {1,1};
        hsize_t block[2] = {1,1};
        dims[0] = ncells;
        dims[1] = 4;
        dataspace_id = H5Screate_simple(2,dims,NULL);
        dataset_id = H5Dcreate(ProcMeshGroup_id,"Connectivity",H5T_NATIVE_INT,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,plist_id);
        hid_t memspace_id = H5Screate_simple(2,dims,NULL);
        int temp_int = 0;
        for (int n = 0; n < ntypes; n++)
        {
            int nvert = nv[n];
            switch (nvert)
            {
            case 3:
                for (int i = 0; i < nc[n]; i++)
                {
                    offset[0] = temp_int++;
                    hid_t status = H5Sselect_hyperslab(dataspace_id,H5S_SELECT_SET,offset,stride,count,block);
                    for(int j = 0;j<3;j++)
                    {
                        buffer[j] = vconn[n][nvert*i+j];
                    }
                    buffer[3] = vconn[n][nvert*i+2];
                    //status = H5Dwrite(dataset_id,H5T_NATIVE_INT,memspace_id,dataspace_id,H5P_DEFAULT,buffer);
                }
                break;
            case 4:
                for (int i = 0; i < nc[n]; i++)
                {
                    offset[0] = temp_int++;
                    hid_t status = H5Sselect_hyperslab(dataspace_id,H5S_SELECT_SET,offset,stride,count,block);
                    for(int j = 0;j<4;j++)
                    {
                        buffer[j] = vconn[n][nvert*i+j];
                    }
                    //status = H5Dwrite(dataset_id,H5T_NATIVE_INT,memspace_id,dataspace_id,H5P_DEFAULT,buffer);
                }
                break;
            default:
                break;
            }
        }
        H5Sclose(memspace_id);
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);


    }
    else if(d_dim==3)
    {
        throw std::runtime_error("Not implemented yet\n");
    }

    H5Gclose(MeshGroup_id);
    H5Gclose(ProcMeshGroup_id);
    H5Fclose(file_id);
    H5Pclose(plist_id);
}