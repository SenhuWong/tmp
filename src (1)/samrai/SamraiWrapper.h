#pragma once
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
//#include "SAMRAI/geom/CartesianGridGeometry.h"
//This is a wrapper class of SAMRAI.
#include "PureGeoIntegrator.h"
class WrapSAMRAI
{
    friend class Orchestra;
    //All the components that need to keep track of.
    std::shared_ptr<SAMRAI::geom::CartesianGridGeometry> grid_geometry;
    std::shared_ptr<SAMRAI::hier::PatchHierarchy> patch_hierarchy;
    PureGeoStrategy *geo_model = NULL;
    OversetStrategy *geo_overset = NULL;
    InterpStrategy *geo_interp = NULL;

    std::shared_ptr<PureGeometricIntegrator> geo_integrator;
    std::shared_ptr<SAMRAI::mesh::StandardTagAndInitialize> error_detector;
    std::shared_ptr<SAMRAI::mesh::BergerRigoutsos> box_generator;
    std::shared_ptr<SAMRAI::mesh::TreeLoadBalancer> load_balancer;
    std::shared_ptr<SAMRAI::mesh::GriddingAlgorithm> gridding_algorithm;
    std::shared_ptr<SAMRAI::appu::VisItDataWriter> visit_data_writer;
    std::shared_ptr<TecplotWriter> tec_writer;
    std::vector<int> tag_buffer_array;
    std::vector<double> regrid_start_time;
    //For output
    bool viz_dump_data = true;

    //For input and restart
    std::string input_filename;
    std::string restart_read_dirname;
    int restore_num;
    bool is_from_restart;

public:
    WrapSAMRAI();
    WrapSAMRAI(int argc, char *argv[]);
    void getFromInput(const std::shared_ptr<SAMRAI::tbox::Database> &input_db, bool is_from_restart);
    void init(SamraiFeeder *feeder);

    void make_geometric_refine();

    void pack_patch_for_overset(
        int *nf, int *qstride, double **qnodein,
        int **idata, double **rdata,
        int *ngridsin, int *qnodesize,
        int ***ibl, double ***q,int** block_global_id, int *nblockin);

    ~WrapSAMRAI();

    void extractGrid(
        int *nf, int *qstride, double **qnodein,
        int **idata, double **rdata,
        int *ngridsout, int *qnodesize);

    void fillCartGridData(int *idata, double *rdata, int ngrids);
    void dump_data(double time);
};