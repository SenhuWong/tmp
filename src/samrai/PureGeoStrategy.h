#pragma once
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/xfer/CoarsenPatchStrategy.h"
#include "SAMRAI/tbox/Dimension.h"

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/appu/VisItDataWriter.h"
#include <string>
#include "tecplotWriter.h"
class PureGeometricIntegrator;

class PureGeoStrategy : public SAMRAI::xfer::RefinePatchStrategy,
                        public SAMRAI::xfer::CoarsenPatchStrategy
{

public:
    
    //These are problem specific stuff need to be done
    void initializeDataOnPatch(SAMRAI::hier::Patch &patch,
                               const double time,
                               const bool initial_time) const;

    void tagGradientDetectorCells(
        SAMRAI::hier::Patch &patch,
        const double regrid_time,
        const bool initial_error,
        const int tag_index,
        const bool uses_richardson_extrapolation_too);
    
    //These are Interface to be implemented for R/C PatchStrategy
    void setPhysicalBoundaryConditions(
        SAMRAI::hier::Patch &patch,
        const double fill_time,
        const SAMRAI::hier::IntVector &ghost_width_tofill);

    SAMRAI::hier::IntVector getRefineOpStencilWidth(
        const SAMRAI::tbox::Dimension &dim) const;

    void preprocessRefine(
        SAMRAI::hier::Patch &fine,
        const SAMRAI::hier::Patch &coarse,
        const SAMRAI::hier::Box &fine_box,
        const SAMRAI::hier::IntVector &ratio);

    void postprocessRefine(
        SAMRAI::hier::Patch &fine,
        const SAMRAI::hier::Patch &coarse,
        const SAMRAI::hier::Box &fine_box,
        const SAMRAI::hier::IntVector &ratio);

    SAMRAI::hier::IntVector getCoarsenOpStencilWidth(
        const SAMRAI::tbox::Dimension &dim) const;

    void preprocessCoarsen(
        SAMRAI::hier::Patch &coarse,
        const SAMRAI::hier::Patch &fine,
        const SAMRAI::hier::Box &coarse_box,
        const SAMRAI::hier::IntVector &ratio);

    void postprocessCoarsen(
        SAMRAI::hier::Patch &coarse,
        const SAMRAI::hier::Patch &fine,
        const SAMRAI::hier::Box &coarse_box,
        const SAMRAI::hier::IntVector &ratio);
    //Trivial Method so Strategy could refer to Integrator's context
    void set_Context(std::shared_ptr<SAMRAI::hier::VariableContext> current,
                     std::shared_ptr<SAMRAI::hier::VariableContext> scratch)
    {
        d_current_one = current;
        d_scratch_one = scratch;
    }

public:
    PureGeoStrategy(const std::string &object_name,
                    const SAMRAI::tbox::Dimension &dim,
                    std::shared_ptr<SAMRAI::tbox::Database> input_db,
                    std::shared_ptr<SAMRAI::geom::CartesianGridGeometry> grid_geometry);
    ~PureGeoStrategy();
    void registerModelVariables(PureGeometricIntegrator *integrator);

    std::string d_object_name;
    SAMRAI::tbox::Dimension d_dim;
    //Computating Variable
    std::shared_ptr<SAMRAI::pdat::CellVariable<double>> d_XSquare;
    std::shared_ptr<SAMRAI::pdat::CellVariable<double>> d_YSquare;
    
    std::shared_ptr<SAMRAI::pdat::CellVariable<int>> d_iblank;
    
    std::shared_ptr<SAMRAI::geom::CartesianGridGeometry> d_grid_geometry;
#ifdef HAVE_HDF5
    void registerVisItDataWriter(std::shared_ptr<SAMRAI::appu::VisItDataWriter> viz_writer);
    std::shared_ptr<SAMRAI::appu::VisItDataWriter> d_visit_writer;
#endif
    void registerTecplotWriter(std::shared_ptr<TecplotWriter> tec_writer);
    std::shared_ptr<TecplotWriter> d_tecplot_writer;

    // //FLOW REFLECT DIRICHLET NEUMANN
    // std::vector<int> d_scalar_bdry_edge_conds;
    // std::vector<int> d_scalar_bdry_node_conds;
    // //Post processed for scalar or vector, passes into boundary routines
    // std::vector<int> d_node_bdry_edge;
    std::shared_ptr<SAMRAI::hier::VariableContext> d_current_one;
    std::shared_ptr<SAMRAI::hier::VariableContext> d_scratch_one;

    SAMRAI::tbox::ResourceAllocator d_allocator;

    //Boundary values for dirichlet case
    std::vector<double> d_bdry_edge_val;
private:
    
};