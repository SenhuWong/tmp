#include <memory>
#pragma once
// amrFeeder should serve as an intermidiator which is maintained by SAMRAI, and used by TIOGA and SAMRAI,
// Every time Regrid is performed in SAMRAI, AMRFeeder's mesh info should be updated.
// 
class AMRFeeder
{
private:
    /* data */
    int d_ncartGrids = 0; //Number of global Cart Grids
    int d_ncartBlocks = 0;//Number of local Cart Block
    int *idata = NULL;
    double *rdata = NULL;
    // Iblank is the flag telling if a cell is tagged as donor or receiver
    // It should be provided by SAMRAI(Updated every time regrid is performed)
    int nf = 0;
    int qstride = 0;

    int qnodesize = 0;
    double* qnodein = NULL;

    int** iblks = NULL;
    // It should be provided by SAMRAI with its PatchDatas all compressed in one. 
    double** qblks = NULL;
    // 
public:
    AMRFeeder(/* args */);
    ~AMRFeeder();
};

AMRFeeder::AMRFeeder(/* args */)
{

}

AMRFeeder::~AMRFeeder()
{

}
