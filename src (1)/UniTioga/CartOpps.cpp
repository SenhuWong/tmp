#include "codetypes.h"
#include "MeshBlock.h"
#include "parallelComm.h"
#include "CartGrid.h"
#include "assert.h"
int obbIntersectCheck(int dim, double vA[3][3], double xA[3], double dxA[3],
                      double vB[3][3], double xB[3], double dxB[3]);
void computeNodalWeights(int dim, double xv[8][3], double *xp, double *frac, int nvert);
void deallocateLinkList3(INTEGERLIST2 *temp);
void get_amr_index_xyz(int dim, int nq, int i, int j, int k,
                       int pBasis,
                       int nX, int nY, int nZ,
                       int nf,
                       double *xlo,
                       double *dx,
                       double qnodes[],
                       int *index,
                       double *xyz);
//Set iblank or iblank_cell depending on the final result of valid donor
void MeshBlock::setCartIblanks()
{
    if (ihigh)
    {
        int m = 0;
        for (int i = 0; i < nreceptorCellsCart; i++)
        {

            int icount = 0;
            for (int j = 0; j < pointsPerCell[i]; j++)
            {
                if (donorIdCart[m++] != -1)
                    icount++;
            }
            //If all points inside this cell is valid
            if (icount == pointsPerCell[i])
                iblank_cell[ctag_cart[i]] = -1;
        }
    }
    else
    {
        int m = 0;
        for (int i = 0; i < nnodes; i++)
        {
            if (pickedCart[i] == 1)
            {
                //If the donor from Cart is valid and not cancelled, then tag iblank here to be -1;
                if (donorIdCart[m++] != -1)
                    iblank[i] = -1;
            }
        }
    }
}

void MeshBlock::findInterpListCart()
{
    if (interpListCart)
    {
        for (int i = 0; i < interpListCartSize; i++)
        {
            if (interpListCart[i].inode)
                TIOGA_FREE(interpListCart[i].inode);
            if (interpListCart[i].weights)
                TIOGA_FREE(interpListCart[i].weights);
        }
        TIOGA_FREE(interpListCart);
        interpListCartSize = 0;
    }
    //Here comes the result from this mb to cartBlocks
    
    for (int i = 0; i < nsearch; i++)
    {
        if (donorId[i] != -1)
            interpListCartSize++;
    }
    interpListCart = new INTERPLIST[interpListCartSize];
    int interpcount = 0;
    double xp[3];
    int inode[8];
    double xv[8][3];
    double frac[8];
    
    for (int i = 0; i < nsearch; i++)
    {
        if (donorId[i] != -1)
        {
            int i3 = 3 * i;
            int procid = isearch[i3];
            int local_id = isearch[i3 + 1];
            int point_id = isearch[i3 + 2];
            for (int j = 0; j < d_dim; j++)
            {
                xp[j] = xsearch[d_dim * i + j];
            }
            int isum = 0;
            //Locate the donor
            int n = 0;
            int temp_count = 0;
            for (; n < ntypes; n++)
            {
                isum += nc[n];
                if (donorId[i] < isum)
                {
                    temp_count = donorId[i] - (isum - nc[n]);
                    break;
                }
            }
            int nvert = nv[n];
            for (int m = 0; m < nvert; m++)
            {
                inode[m] = vconn[n][nvert * temp_count + m] - BASE;
                for (int j = 0; j < d_dim; j++)
                {
                    xv[m][j] = x[inode[m] * d_dim + j];
                }
            }
            computeNodalWeights(d_dim, xv, xp, frac, nvert);

            interpListCart[interpcount].receptorInfo[0] = procid;
            interpListCart[interpcount].receptorInfo[1] = point_id;
            interpListCart[interpcount].receptorInfo[1] = local_id;
            interpListCart[interpcount].nweights = nvert;
            interpListCart[interpcount].cancel = 0;

            interpListCart[interpcount].inode = new int[nvert];
            interpListCart[interpcount].weights = new double[nvert];
            for (int m = 0; m < nvert; m++)
            {
                interpListCart[interpcount].inode[m] = inode[m];
                interpListCart[interpcount].weights[m] = frac[m];
            }
            interpcount++;
        }
    }
    ninterpCart = interpcount;
    assert(ninterpCart == interpcount);
}
//Get receptors on  Cartesian Grids
void MeshBlock::getCartReceptors(CartGrid *cg, parallelComm *pc)
{
    int iflag;
    OBB *obCart = new OBB[1];
    double xd[3];
    int *pmap = new int[pc->numprocs]; //Showing if a proc has intersecting CartBlocks
    for (int i = 0; i < pc->numprocs; i++)
    {
        pmap[i] = 0;
    }
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            obCart->vec[i][j] = 0;
        }
    }
    obCart->vec[0][0] = obCart->vec[1][1] = obCart->vec[2][2] = 1.0;

    INTEGERLIST2 *head = new INTEGERLIST2[1];
    INTEGERLIST2 *curPtr;
    head->intData = NULL;
    head->realData = NULL;
    curPtr = head;
    curPtr->next = NULL;
    nsearch = 0;
    //Make a basic OBB to check for each CartGrid
    for (int c = 0; c < cg->ngrids; c++)
    {
        for (int n = 0; n < d_dim; n++)
        {
            obCart->dxc[n] = (cg->dx[d_dim * c + n]) * (cg->ncells[d_dim * c + n]) * 0.5;
            obCart->xc[n] = cg->xlo[d_dim * c + n] + obCart->dxc[n];
        }
        if (obbIntersectCheck(d_dim, obb->vec, obb->xc, obb->dxc, obCart->vec, obCart->xc, obCart->dxc) or obbIntersectCheck(d_dim, obCart->vec, obCart->xc, obCart->dxc, obb->vec, obb->xc, obb->dxc))
        {
            int ntm = (cg->porder[c] + 1) * (cg->porder[c] + 1) * (d_dim == 2 ? 1 : cg->porder[c] + 1);
            double *xtm = new double[d_dim * ntm];
            int *itm = new int[ntm];
            int ploc = (cg->porder[c]) * (cg->porder[c] + 1) / 2; //Used to refer to the p-th order's x loc
            if (d_dim == 2)
            {
                for (int j = 0; j < cg->ncells[2 * c]; j++)
                {
                    for (int k = 0; k < cg->ncells[2 * c + 1]; k++)
                    {
                        get_amr_index_xyz(d_dim, cg->qstride,
                                          j, k, 0, cg->porder[c],
                                          cg->ncells[2 * c], cg->ncells[2 * c + 1], 0, cg->nf, &cg->xlo[2 * c], &cg->dx[2 * c], &cg->qnode[ploc], itm, xtm);
                        iflag = 0;
                        for (int n = 0; n < ntm; n++)
                        {
                            int i2 = 2 * n;
                            for (int jj = 0; jj < 2; jj++)
                            {
                                xd[jj] = 0;
                            }
                            for (int jj = 0; jj < 2; jj++)
                            {
                                for (int kk = 0; kk < 2; kk++)
                                {
                                    xd[jj] += (xtm[i2 + kk] - obb->xc[kk]) * obb->vec[jj][kk];
                                }
                            }
                            if (fabs(xd[0]) <= obb->dxc[0] and fabs(xd[1]) <= obb->dxc[1])
                            {
                                iflag++;
                            }
                        }
                        if (iflag > 0)
                        {
                            pmap[cg->proc_id[c]] = 1;
                            curPtr->next = new INTEGERLIST2[1];
                            curPtr = curPtr->next;
                            curPtr->realData = new double[3 * ntm];
                            curPtr->intData = new int[ntm + 2];
                            curPtr->intData[0] = cg->proc_id[c];
                            curPtr->intData[1] = cg->local_id[c];
                            curPtr->intDataSize = ntm + 2;
                            curPtr->realDataSize = ntm;
                            nsearch += ntm;
                            for (int n = 0; n < ntm; n++)
                            {
                                for (int kk = 0; kk < 2; kk++)
                                {
                                    curPtr->realData[2 * n + kk] = xtm[2 * n + kk];
                                }
                                curPtr->intData[2 + n] = itm[n];
                            }
                            curPtr->next = NULL;
                        }
                    }
                }
            }
            else if (d_dim == 3)
            {
                for (int j = 0; j < cg->ncells[3 * c]; j++)
                {
                    for (int k = 0; k < cg->ncells[3 * c + 1]; k++)
                    {
                        for (int l = 0; l < cg->ncells[3 * c + 2]; l++)
                        {
                            //For this cells get all amr index
                            get_amr_index_xyz(d_dim, cg->qstride, j, k, l, cg->porder[c],
                                              cg->ncells[3 * c], cg->ncells[3 * c + 1], cg->ncells[3 * c + 2],
                                              cg->nf, &cg->xlo[3 * c], &cg->dx[3 * c], &cg->qnode[ploc], itm, xtm);
                            iflag = 0;
                            for (int n = 0; n < ntm; n++)
                            {
                                int i3 = 3 * n;
                                for (int jj = 0; jj < 3; jj++)
                                {
                                    xd[jj] = 0;
                                }
                                for (int jj = 0; jj < 3; jj++)
                                {
                                    for (int kk = 0; kk < 3; kk++)
                                    {
                                        xd[jj] += (xtm[i3 + kk] - obb->xc[kk]) * obb->vec[jj][kk];
                                    }
                                }
                                if (fabs(xd[0]) <= obb->dxc[0] and fabs(xd[1]) <= obb->dxc[1] and fabs(xd[2]) <= obb->dxc[2])
                                {
                                    iflag++;
                                }
                            }
                            if (iflag > 0)
                            {
                                pmap[cg->proc_id[c]] = 1;
                                curPtr->next = new INTEGERLIST2[1];
                                curPtr = curPtr->next;
                                curPtr->realData = new double[3 * ntm];
                                curPtr->intData = new int[ntm + 2];
                                curPtr->intData[0] = cg->proc_id[c];
                                curPtr->intData[1] = cg->local_id[c];
                                curPtr->intDataSize = ntm + 2;
                                curPtr->realDataSize = ntm;
                                nsearch += ntm;
                                for (int n = 0; n < ntm; n++)
                                {
                                    for (int kk = 0; kk < 3; kk++)
                                    {
                                        curPtr->realData[3 * n + kk] = xtm[3 * n + kk];
                                    }
                                    curPtr->intData[n + 2] = itm[n];
                                }
                                curPtr->next = NULL;
                            }
                        }
                    }
                }
            }
            TIOGA_FREE(xtm);
            TIOGA_FREE(itm);
        }
    }
    //Create the communication map
    int nsend = 0;
    for (int i = 0; i < pc->numprocs; i++)
    {
        if (pmap[i] == 1)
        {
            nsend++;
        }
    }
    //Assume there is always a overlap because Cart is background
    int nrecv = nsend;
    int *sndMap = new int[nsend];
    int *rcvMap = new int[nrecv];
    int m = 0;
    int n = 0;
    for (int i = 0; i < pc->numprocs; i++)
    {
        if (pmap[i] == 1)
        {
            sndMap[m] = rcvMap[m] = i;
            m++;
        }
    }
    pc->setMap(nsend, nrecv, sndMap, rcvMap);
    if (xsearch)
        TIOGA_FREE(xsearch);
    if (isearch)
        TIOGA_FREE(isearch);
    if (donorId)
        TIOGA_FREE(donorId);
    if (rst)
        TIOGA_FREE(rst);
    xsearch = new double[d_dim * nsearch];
    isearch = new int[3 * nsearch];
    donorId = new int[nsearch];
    rst = new double[d_dim * nsearch];

    curPtr = head->next;
    m = 0;
    n = 0;
    while (curPtr != NULL)
    {
        for (int j = 2; j < curPtr->intDataSize; j++)
        {
            isearch[m++] = curPtr->intData[0];
            isearch[m++] = curPtr->intData[1];
            isearch[m++] = curPtr->intData[j];
        }
        for (int j = 0; j < d_dim * curPtr->realDataSize; j++)
        {
            xsearch[n++] = curPtr->realData[j];
        }
        curPtr = curPtr->next;
    }
    deallocateLinkList3(head);
    TIOGA_FREE(obCart);
    TIOGA_FREE(pmap);
    TIOGA_FREE(sndMap);
    TIOGA_FREE(rcvMap);
}

void MeshBlock::getUnresolvedMandatoryReceptors()
{

    int *iflag = new int[ncells];
    if (pickedCart)
        TIOGA_FREE(pickedCart);
    pickedCart = new int[nnodes];
    for (int i = 0; i < ncells; i++)
        iflag[i] = 0;
    for (int i = 0; i < nnodes; i++)
        pickedCart[i] = 0;
    int temp_count = 0;
    int temp_int = 0;
    int inode[8];

    for (int n = 0; n < ntypes; n++)
    {
        int nvert = nv[n];
        for (int i = 0; i < nc[n]; i++)
        {
            temp_count = 0;
            for (int m = 0; m < nvert; m++)
            {
                inode[m] = vconn[n][nvert * i + m] - BASE;
                if (nodeRes[inode[m]] == BIGVALUE)
                    temp_count++;
            }
            if (temp_count == nvert and iblank_cell[temp_int] == 1)
            {
                iflag[temp_int] = 1;
                for (int m = 0; m < nvert; m++)
                {
                    pickedCart[inode[m]] = 1;
                }
            }
            temp_int++;
        }
    }

    if (ctag_cart)
        TIOGA_FREE(ctag_cart);
    ctag_cart = new int[ncells];
    nreceptorCellsCart = 0;
    for (int i = 0; i < ncells; i++)
    {
        //==1 i guess
        if (iflag[i] == -1)
            ctag_cart[nreceptorCellsCart++] = i + 1; //Starting from 1
    }

    if (ihigh)
    {
        if (pointsPerCell)
            TIOGA_FREE(pointsPerCell);
        pointsPerCell = new int[nreceptorCellsCart];
        maxPointsPerCell = 0;
        ntotalPointsCart = 0;
        for (int i = 0; i < nreceptorCellsCart; i++)
        {
            get_nodes_per_cell(&(ctag_cart[i]), &pointsPerCell[i]);
            ntotalPointsCart += pointsPerCell[i];
            maxPointsPerCell = TIOGA_MAX(maxPointsPerCell, pointsPerCell[i]);
        }
        if (rxyzCart)
            TIOGA_FREE(rxyzCart);
        if (donorIdCart)
            TIOGA_FREE(donorIdCart);
        rxyzCart = new double[ntotalPointsCart * d_dim];
        donorIdCart = new int[ntotalPointsCart];
        temp_int = 0;
        for (int i = 0; i < nreceptorCellsCart; i++)
        {
            get_receptor_nodes(&(ctag_cart[i]), &(pointsPerCell[i]), &(rxyzCart[temp_int]));
            temp_int += (d_dim * pointsPerCell[i]);
        }
    }
    else
    {
        ntotalPointsCart = 0;
        for (int i = 0; i < nnodes; i++)
        {
            if (pickedCart[i] == 1)
            {
                ntotalPointsCart++;
            }
        }

        if (rxyzCart)
            TIOGA_FREE(rxyzCart);
        if (donorIdCart)
            TIOGA_FREE(donorIdCart);
        if (receptorIdCart)
            TIOGA_FREE(receptorIdCart);
        rxyzCart = new double[d_dim * ntotalPointsCart];
        donorIdCart = new int[ntotalPointsCart];
        receptorIdCart = new int[ntotalPointsCart];
        temp_int = 0;
        temp_count = 0;

        for (int i = 0; i < nnodes; i++)
        {
            if (pickedCart[i] == 1)
            {
                receptorIdCart[temp_int++] = i;

                int id = d_dim * i;
                for (int j = 0; j < d_dim; j++)
                {
                    rxyzCart[temp_count++] = x[id + j];
                }
            }
        }
    }
    TIOGA_FREE(iflag);
}