#pragma once
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SamraiFeeder.h"
#include <memory>
//This class manages the behaviour like
//1.Ask the Integrator for a overset related Patch and get the output.
//
//2.Zip all the CellData on that Patch to fit into TIOGA's interface
//
//3.Procedure 2 should be done after scratch workspace is filled in PureGeoStrategy.
//
//4.We might want to use PureGeoStrategy's Variable and Context,should we use friend or ask Integrator for it?
//
//5.What does this class need to know to performwell?
//qnode
//qnodesize
//nf(not here)
//
class PureGeometricIntegrator;
class InterpStrategy
{
    struct DataPack
    {
        int idata[11];   //9 in 2D,11 in 3D
        double rdata[6]; //4 in 2D,6 in 3D
        int *ibl;        //This should be get from PureGeoMetric
        double *q;       //This could be easily get from here
    };
    std::string d_object_name;

    SAMRAI::tbox::Dimension d_dim;
    std::shared_ptr<SAMRAI::geom::CartesianGridGeometry> d_grid_geometry;

    SAMRAI::tbox::ResourceAllocator d_allocator;
    //The valid(1.related to body 2.over some minimal refinement time)
    //The first is done through OBB Check
    //The later is easily done by selecting level.
    int cur_proc = -1;
    int num_proc = -1;

    //Obtained from input_db
    int qnodesize = -1;
    int least_refine = -1;
    double *qnode = NULL;

    //Variable is dynamically assigned a depth.
    std::shared_ptr<SAMRAI::pdat::CellVariable<double>> zip_variable;
    int ZipDepth = -1;
    int ngridsin = -1;
    int *ZipIndexes = NULL;
    std::list<DataPack> related_patches;

    //Set by Integrator's Constructor
    std::shared_ptr<SAMRAI::hier::VariableContext> workspace;

public:
    InterpStrategy();
    InterpStrategy(const std::string &object_name,
                   const SAMRAI::tbox::Dimension &dim,
                   std::shared_ptr<SAMRAI::tbox::Database> input_db,
                   std::shared_ptr<SAMRAI::geom::CartesianGridGeometry> grid_geometry);
    
    void getFromInput(const std::shared_ptr<SAMRAI::tbox::Database> &input_db, bool is_from_restart);
    inline void set_Context(std::shared_ptr<SAMRAI::hier::VariableContext> current_context)
    {
        workspace = current_context;
    }
    inline int geo_fine_level()
    {
        return least_refine;
    }
    //Get all the context and Variables on PureGeoStrategy.
    void registerModelVariables(PureGeometricIntegrator *integrator);
    void registerFeeder(SamraiFeeder *fder)
    {
        //Find a list of Bounding Box from reader for later use.
    }
    //Counting could totally be done on integrator and oversetStrategy side.
    //Now all integrator has to do is pass a patch in and get its variable pressed into TIOGA format

    //For now ,just make a cheap "depthed" version of all these variables, later we might need to implement a PatchData
    void make_amr(SAMRAI::hier::Patch &patch);
    double *make_zip(SAMRAI::hier::Patch &patch, double zip_data_time);
    int *make_ibl(SAMRAI::hier::Patch &patch, PureGeometricIntegrator *integrator, double make_block_time);
    void formCartBlock(SAMRAI::hier::Patch &patch, PureGeometricIntegrator *integrator, double make_block_time);
    void formCartGrid(int *nf, int *qstride, double **qnodein,
     int *ngridsin, int *qnodeSize,
      int **idata, double **rdata, int** block_global_id);
    void takeCartBlock(double ***q, int ***ibl, int *nblocks);
};