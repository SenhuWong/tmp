#include"codetypes.h"
#include"MeshBlock.h"
#include<iostream>
void find_median_wrap(int* ix, double* x, int n);

//Recursive functions used by ADT class;
void buildADTrecursion(double* coord, double* adtReals, double* adtWork, int* adtIntegers,
    int* elementsAvailable, int* adtCount, int side, int parent,
    int level, int ndim, int nelem, int nav)
{
    int nd = ndim / 2;
    int i, j;
    int dimcut;
    int nleft;
    int ii, iip, jj, jjp;
    int parentToChild;

    if (nav > 1)//Call recursion
    {
        // find the dimension to create the cut
        dimcut = (level % ndim);
        //
        // collect coordinates along the dimension dimcut
        //
        for (i = 0; i < nav; i++)
            adtWork[i] = coord[ndim * elementsAvailable[i] + dimcut];
        //
        // reorder elements with nleft elements to
        // the left of median of adtWork
        //
        find_median_wrap(elementsAvailable, adtWork, nav);
        nleft = (nav + 1) / 2;
        (*adtCount)++;
        ii = (*adtCount) * 4;
        //nleft is placed at the node;
        adtIntegers[ii] = elementsAvailable[nleft - 1];
        adtIntegers[ii + 1] = -1;//left
        adtIntegers[ii + 2] = -1;//right
        adtIntegers[ii + 3] = -1;//This is to its parent;
        //
        // find minimum and maximum bounds of the elements
        // contained in this leaf
        //
        for (i = 0; i < nd; i++)
        {
            adtReals[ndim * (*adtCount) + i] = BIGVALUE;//The first nd is min
            adtReals[ndim * (*adtCount) + i + nd] = -BIGVALUE;//The last nd is max
        }
        //
        for (i = 0; i < nav; i++)
            for (j = 0; j < nd; j++)
            {
                ii = ndim * (*adtCount) + j;
                iip = ii + nd;
                jj = ndim * elementsAvailable[i] + j;
                jjp = jj + nd;
                //
                adtReals[ii] = TIOGA_MIN(adtReals[ii], coord[jj]);
                adtReals[iip] = TIOGA_MAX(adtReals[iip], coord[jjp]);
            }
        //
        // specify that the new element is the child of parent
        // unless root
        //
        if (side > 0)
        {
            adtIntegers[4 * parent + side] = elementsAvailable[nleft - 1];
        }
        parentToChild = *adtCount;
        //
        // build the left side of the tree
        //
        if (nleft > 1)
        {
            buildADTrecursion(coord, adtReals, adtWork, adtIntegers, elementsAvailable,
                adtCount, 1, parentToChild, level + 1, ndim, nelem, nleft - 1);//side is 1 when left child
        }
        //
        // build the right side of the tree
        //
        buildADTrecursion(coord, adtReals, adtWork, adtIntegers, &(elementsAvailable[nleft]),
            adtCount, 2, parentToChild, level + 1, ndim, nelem, nav - nleft);//side is 2 when right child

    }
    else if (nav == 1)//The only one left is placed at the node;
    {
        (*adtCount)++;
        ii = 4 * (*adtCount);
        jj = ndim * (*adtCount);
        adtIntegers[ii] = elementsAvailable[0];//at index zero is the current element Bbox indice inside elementBox
        adtIntegers[ii + 1] = -1;
        adtIntegers[ii + 2] = -1;
        adtIntegers[ii + 3] = -1;
        for (j = 0; j < ndim; j++)
        {
            adtReals[jj + j] = coord[ndim * elementsAvailable[0] + j];//get the coordinate according to the currrent element Bbox indice
        }
        if (side > 0)
        {
            //side decide which side the current tree node lies at its parent side;
            adtIntegers[4 * parent + side] = elementsAvailable[0];
        }
    }

}

void searchIntersections(MeshBlock* mb, int* cellIndex, int* adtIntegers, double* adtReals,
    double* coord, int level, int node, double* xsearch, int nelem, int ndim)
{
    int i;
    int d, nodeChild;
    double* element = new double[ndim];
    bool flag;
    //
    for (i = 0; i < ndim; i++)
        element[i] = coord[ndim * (adtIntegers[4 * node]) + i];
    //
    flag = 1;
    for (i = 0; i < ndim / 2; i++)
        flag = (flag && (xsearch[i] >= element[i] - TOL));
    for (i = ndim / 2; i < ndim; i++)
        flag = (flag && (xsearch[i - ndim / 2] <= element[i] + TOL));
    //
    if (flag)
    {
        
        mb->checkContainment(cellIndex, adtIntegers[4 * node], xsearch);
        if (cellIndex[0] > -1 && cellIndex[1] == 0) return;
    }
    //
    // check the left and right children
    // now
    //
    for (d = 1; d < 3; d++)
    {
        nodeChild = adtIntegers[4 * node + d];
        if (nodeChild > -1)
        {
            nodeChild = adtIntegers[4 * nodeChild + 3];
            for (i = 0; i < ndim; i++)
            {
                element[i] = adtReals[ndim * nodeChild + i];
            }
            flag = 1;
            for (i = 0; i < ndim / 2; i++)
                flag = (flag && (xsearch[i] >= element[i] - TOL));
            for (i = ndim / 2; i < ndim; i++)
                flag = (flag && (xsearch[i - ndim / 2] <= element[i] + TOL));
            if (flag)
            {
                searchIntersections(mb, cellIndex, adtIntegers, adtReals, coord, level + 1,
                    nodeChild, xsearch, nelem, ndim);
                if (cellIndex[0] > -1 && cellIndex[1] == 0) return;
            }
        }
    }
    delete[] element;
    return;
}


void searchNearIntersections(double* minD, double* x_search, int* adtIntegers, double* adtReals,
    double* coord, int level, int node, int nelem, int ndim, int cur_proc)
{
    double element[6];//This is a point ,not a cell
    //First check if the current node is inside the search_box
    for (int i = 0; i < ndim / 2; i++)
    {
        element[i] = coord[ndim * adtIntegers[4 * node] + i];
    }
    bool flag = true;
    for (int i = 0; i < ndim / 2; i++)
    {
        flag = flag and (abs(element[i] - x_search[i]) <= *minD);
    }
    //If inside ,update the search_box when a new minimal is found
    //Else nothing changes
    if (flag)
    {
        double new_minD = 0;
        for (int i = 0; i < ndim / 2; i++)
        {
            new_minD += pow(x_search[i] - element[i], 2);
        }
        new_minD = sqrt(new_minD);
        if (new_minD < *minD)
        {
            *minD = new_minD;
        }

    }

    
    //Check the left and right subrigion for overlap with the search_box
    //If overlap exists, call searchNearIntersections at the next level
    for (int d = 1; d < 3; d++)
    {
        int nodeChild = adtIntegers[4 * node + d];
        if (nodeChild > -1)
        {
            nodeChild = adtIntegers[4 * nodeChild + 3];
            for (int i = 0; i < ndim; i++)
            {
                element[i] = adtReals[ndim * nodeChild + i];
            }
            flag = true;
            for (int i = 0; i < ndim / 2; i++)
            {
                flag = flag and (element[i] <= x_search[i] + *minD);
                flag = flag and (element[ndim / 2 + i] >= x_search[i] - *minD);
            }
            if (flag)
            {
                //std::cout << level << '\t' << level + 1 << '\n';
                searchNearIntersections(minD, x_search, adtIntegers, adtReals,
                    coord, level + 1, nodeChild, nelem, ndim, cur_proc);
            }
        }
    }

}