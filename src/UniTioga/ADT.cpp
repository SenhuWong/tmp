#include"ADT.h"
#include"MeshBlock.h"

void buildADTrecursion(double* coord, double* adtReals, double* adtWork, int* adtIntegers,
    int* elementsAvailable, int* adtCount, int side, int parent,
    int level, int ndim, int nelem, int nav);

void searchIntersections(MeshBlock* mb, int* cellIndex, int* adtIntegers, double* adtReals,
    double* coord, int level, int node, double* xsearch, int nelem, int ndim);

void searchNearIntersections(double* minD, double* x_search, int* adtIntegers, double* adtReals,
    double* coord, int level, int node, int nelem , int ndim, int cur_proc);
//ADT class methods;

void ADT::buildADT(int d, int nelements, double* elementBbox)
{
    int* elementsAvailable;
    double* adtWork;
    int adtCount, parent, level, nav;
    int side;
    double tolerance, delta;

    ndim = d;
    nelem = nelements;
    coord = elementBbox;//First 4 is min ,last 4 is max

    //Allocate work arrays
    elementsAvailable = new int[nelem];
    adtWork = new double[nelem];

    //Allocate arrays in the class

    if (adtExtents) TIOGA_FREE(adtExtents);
    if (adtIntegers) TIOGA_FREE(adtIntegers);
    if (adtReals) TIOGA_FREE(adtReals);
    adtExtents = new double[ndim];
    adtIntegers = new int[4 * nelem];
    adtReals = new double[nelem * ndim];

    //Find extent of the elements;
    for (int i = 0; i < ndim / 2; i++)
    {
        int i2 = 2 * i;
        adtExtents[i2] = BIGVALUE;
        adtExtents[i2 + 1] = -BIGVALUE;
    }
    for (int j = 0; j < nelem; j++)
    {
        int jd = ndim * j;
        for (int i = 0; i < ndim / 2; i++)
        {
            int i2 = 2 * i;
            adtExtents[i2] = TIOGA_MIN(adtExtents[i2], coord[jd + i]);            
        }
        for (int i = 0; i < ndim / 2; i++)
        {
            int i2 = 2 * i + 1;
            adtExtents[i2] = TIOGA_MAX(adtExtents[i2], coord[jd + i + ndim / 2]);
        }
    }
    tolerance = 0.01;
    for (int i = 0; i < ndim / 2; i++)
    {
        int i2 = 2 * i;
        delta = tolerance * (adtExtents[i2 + 1] - adtExtents[i2]);
        adtExtents[i2] -= delta;
        adtExtents[i2 + 1] += delta;
    }

    //Build ADT with recursive method
    for (int i = 0; i < nelem; i++)
    {
        elementsAvailable[i] = i;
    }
    adtCount = -1;
    side = 0;
    parent = 0;
    level = 0;
    nav = nelem;
    //make the ptr from parent to its children
    buildADTrecursion(coord, adtReals, adtWork, adtIntegers, elementsAvailable,
        &adtCount, side, parent, level, ndim, nelem, nav);
    //make the ptr from children to its parent
    for (int i = 0; i < nelem; i++)
    {
        int i4 = 4 * adtIntegers[4 * i];
        adtIntegers[i4 + 3] = i;
    }
    TIOGA_FREE(elementsAvailable);
    TIOGA_FREE(adtWork);
    //finally, the ADT structure is stored at adtInteger.
    //the adtReals reaveals the extent of each level
    //the adtWorks has every element Bbox's bounding box.
}

void ADT::searchADT(MeshBlock* mb, int* cellIndex, double* xsearch)
{
    int i;
    int flag;
    int rootNode;
    //
    // check if the given point is in the bounds of
    // the ADT
    //
    rootNode = 0;
    cellIndex[0] = -1;//If it is not inside even the greatest BBOX, it returns with [-1,0];
    cellIndex[1] = 0;
    //
    flag = 1;
    for (i = 0; i < ndim / 2; i++)
        flag = (flag && (xsearch[i] >= adtExtents[2 * i] - TOL));
    for (i = 0; i < ndim / 2; i++)
        flag = (flag && (xsearch[i] <= adtExtents[2 * i + 1] + TOL));
    //
    // call recursive routine to check intersections with 
    // ADT nodes
    //
    if (flag) searchIntersections(mb, cellIndex, adtIntegers, adtReals,
        coord, 0, rootNode, xsearch, nelem, ndim);
}



void ADT::searchNearest(double* minD, double* x_search, int cur_proc)
{
    int rootNode = 0;
    int flag = 1;

    for (int i = 0; i < ndim / 2; i++)
    {
        flag = (flag and (x_search[i] - *minD <= adtExtents[2 * i + 1] + TOL));//box_search's minimal <= region's maximal
    }
    for (int i = 0; i < ndim / 2; i++)
    {
        flag = (flag and (x_search[i] + *minD >= adtExtents[2 * i] - TOL));//box_search's maximal >= regions's minimal
    }
    if (flag)
        searchNearIntersections(minD, x_search, adtIntegers, adtReals, coord, 0, rootNode, nelem, ndim, cur_proc);
}

