#include "PureGeoIntegrator.h"
#include "SAMRAI/tbox/AllocatorDatabase.h"
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/tbox/Utilities.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

PureGeoStrategy::PureGeoStrategy(const std::string &object_name,
                                 const SAMRAI::tbox::Dimension &dim,
                                 std::shared_ptr<SAMRAI::tbox::Database> input_db,
                                 std::shared_ptr<SAMRAI::geom::CartesianGridGeometry> grid_geometry) : SAMRAI::xfer::RefinePatchStrategy(),
                                                                                                       SAMRAI::xfer::CoarsenPatchStrategy(),
                                                                                                       d_object_name(object_name),
                                                                                                       d_dim(dim),
                                                                                                       d_grid_geometry(grid_geometry),
                                                                                                       d_allocator(SAMRAI::tbox::AllocatorDatabase::getDatabase()->getDefaultAllocator()),
                                                                                                       d_XSquare(new SAMRAI::pdat::CellVariable<double>(dim, "XSquareVar", d_allocator)),
                                                                                                       d_YSquare(new SAMRAI::pdat::CellVariable<double>(dim, "YSquareVar", d_allocator)),
                                                                                                       d_iblank(new SAMRAI::pdat::CellVariable<int>(dim,"Iblank",d_allocator))
                                                                            
{
}
PureGeoStrategy::~PureGeoStrategy()
{
}

void simple_initialize_2d(const double *dx, const double *xlo, const double *xhi,
                          int ifirst0, int ilast0, int ifirst1, int ilast1, double *x_val_ptr, double *y_val_ptr)
{
    //std::cout<<"Simple Initialize 2D being called\n";
    int x_stride = ilast0 - ifirst0 + 1;
    int y_stride = ilast1 - ifirst1 + 1;
    for (int i = 0; i <= ilast1 - ifirst1; i++) //Looping through index 1
    {
        for (int j = 0; j <= ilast0 - ifirst0; j++) //Looping through index 0
        {
            double x_loc = xlo[0] + dx[0] * j;
            double y_loc = xlo[1] + dx[1] * i;
            //Wrong way
            //x_val_ptr[y_stride * j + i] = x_loc;
            //y_val_ptr[y_stride * j + i] = y_loc;
            //Right way
            x_val_ptr[x_stride * i + j] = x_loc;
            y_val_ptr[x_stride * i + j] = y_loc;
        }
    }
}

void PureGeoStrategy::initializeDataOnPatch(SAMRAI::hier::Patch &patch,
                                            const double time,
                                            const bool initial_time) const
{
    //std::cout<<"PureGeoStrategy::initializeDataOnPatch being called\n";
    NULL_USE(time);
    if (initial_time)
    {
        const std::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> patch_geom(
            SAMRAI_SHARED_PTR_CAST<SAMRAI::geom::CartesianPatchGeometry, SAMRAI::hier::PatchGeometry>(patch.getPatchGeometry()));
        TBOX_ASSERT(patch_geom);
        const double *dx = patch_geom->getDx();
        const double *xlo = patch_geom->getXLower();
        const double *xhi = patch_geom->getXUpper();

        std::shared_ptr<SAMRAI::pdat::CellData<double>> XSquares(
            SAMRAI_SHARED_PTR_CAST<SAMRAI::pdat::CellData<double>, SAMRAI::hier::PatchData>(patch.getPatchData(d_XSquare, d_current_one)));
        TBOX_ASSERT(XSquares);
        std::shared_ptr<SAMRAI::pdat::CellData<double>> YSquares(
            SAMRAI_SHARED_PTR_CAST<SAMRAI::pdat::CellData<double>, SAMRAI::hier::PatchData>(patch.getPatchData(d_YSquare, d_current_one)));
        TBOX_ASSERT(YSquares);

        const SAMRAI::hier::Box &pbox = patch.getBox();
        const SAMRAI::hier::Index ifirst = pbox.lower();
        const SAMRAI::hier::Index ilast = pbox.upper();

        if (d_dim == SAMRAI::tbox::Dimension(2))
        {
            simple_initialize_2d(dx, xlo, xhi, ifirst(0), ilast(0), ifirst(1), ilast(1), XSquares->getPointer(), YSquares->getPointer());
        }
        else
        {
            TBOX_ASSERT(false);
        }
    }
}
#ifdef HAVE_HDF5
void PureGeoStrategy::registerVisItDataWriter(std::shared_ptr<SAMRAI::appu::VisItDataWriter> viz_writer)
{
    TBOX_ASSERT(viz_writer);
    d_visit_writer = viz_writer;
}
#endif

void PureGeoStrategy::registerTecplotWriter(std::shared_ptr<TecplotWriter> tec_writer)
{
    d_tecplot_writer = tec_writer;
}
void PureGeoStrategy::registerModelVariables(PureGeometricIntegrator *integrator)
{
    std::cout << "registerModelVariable being called\n";
    // //This flag variable telling if it's inside is used only in current context
    // //This is used only in this strategy therefore introduces a Selector
    
    std::cout << "FIrst flag being registered\n";
    //I'd better find a common ghost width first;
    integrator->registerVariable(
        d_XSquare,
        SAMRAI::hier::IntVector::getZero(d_dim),
        d_grid_geometry,
        "CONSERVATIVE_COARSEN",
        "LINEAR_REFINE",
        PureGeometricIntegrator::VariableContextRegistry::CURRENT_SCRATCH);
    integrator->registerVariable(
        d_YSquare,
        SAMRAI::hier::IntVector::getZero(d_dim),
        d_grid_geometry,
        "CONSERVATIVE_COARSEN",
        "LINEAR_REFINE",
        PureGeometricIntegrator::VariableContextRegistry::CURRENT_SCRATCH);
    integrator->registerVariable(
        d_iblank,
        SAMRAI::hier::IntVector::getZero(d_dim),
        d_grid_geometry,
        "NO_COARSEN",
        "NO_REFINE",
        PureGeometricIntegrator::VariableContextRegistry::CURRENT_ONLY);
//     //We might need some HDF5 here

    SAMRAI::hier::VariableDatabase *vardb = SAMRAI::hier::VariableDatabase::getDatabase();

    int xSquare_id = vardb->mapVariableAndContextToIndex(d_XSquare, d_current_one);
    int ySquare_id = vardb->mapVariableAndContextToIndex(d_YSquare, d_current_one);
    int iblank_id = vardb->mapVariableAndContextToIndex(d_iblank,d_current_one);
    
    std::string dump_nameX = "XSquare Var #";
    std::string dump_nameY = "YSquare Var #";
    std::string dump_nameIblank = "IBLANK Var #";
    d_tecplot_writer->register_var(iblank_id);


#ifdef HAVE_HDF5
    if (d_visit_writer)
    {
        d_visit_writer->registerPlotQuantity(dump_nameX, "SCALAR", xSquare_id, 0);
        d_visit_writer->registerPlotQuantity(dump_nameY, "SCALAR", ySquare_id, 0);
        d_visit_writer->registerPlotQuantity(dump_nameIblank, "SCALAR", iblank_id, 0);
    }
    else
    {
        TBOX_WARNING(d_object_name << ":registerModelVariables()\n"
                                   << "VisIt data writer was not registered.\n"
                                   << "Consequently, no plot data will\n"
                                   << "be written" << std::endl);
    }
#endif
}


void PureGeoStrategy::tagGradientDetectorCells(SAMRAI::hier::Patch &patch,
                                               const double regrid_time,
                                               const bool initial_error,
                                               const int tag_index,
                                               const bool uses_richardson_extrapolation_too)
{
    std::cout << "tagGradientDetectorCells being called\n";

    NULL_USE(regrid_time);
    NULL_USE(initial_error);
    NULL_USE(uses_richardson_extrapolation_too);

    const std::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> patch_geom(
        SAMRAI_SHARED_PTR_CAST<SAMRAI::geom::CartesianPatchGeometry, SAMRAI::hier::PatchGeometry>(patch.getPatchGeometry()));
    TBOX_ASSERT(patch_geom);
    const double *dx = patch_geom->getDx();
    const double *xlo = patch_geom->getXLower();
    const double *xhi = patch_geom->getXUpper();

    std::shared_ptr<SAMRAI::pdat::CellData<int>> tags(
        SAMRAI_SHARED_PTR_CAST<SAMRAI::pdat::CellData<int>, SAMRAI::hier::PatchData>(patch.getPatchData(tag_index)));
        
    std::shared_ptr<SAMRAI::pdat::CellData<double>> XSquares(
        SAMRAI_SHARED_PTR_CAST<SAMRAI::pdat::CellData<double>, SAMRAI::hier::PatchData>(patch.getPatchData(d_XSquare, d_current_one)));
    TBOX_ASSERT(XSquares);
    std::shared_ptr<SAMRAI::pdat::CellData<double>> YSquares(
        SAMRAI_SHARED_PTR_CAST<SAMRAI::pdat::CellData<double>, SAMRAI::hier::PatchData>(patch.getPatchData(d_YSquare, d_current_one)));
    TBOX_ASSERT(YSquares);

    const SAMRAI::hier::Box &pbox = patch.getBox();
    const SAMRAI::hier::Index ifirst = pbox.lower();
    const SAMRAI::hier::Index ilast = pbox.upper();

    //int *inside_ptr = inside->getPointer();
    int *tag_ptr = tags->getPointer();
    int vol = (ilast(0) - ifirst(0) + 1) * (ilast(1) - ifirst(1) + 1);
    for (int i = 0; i < vol; i++)
    {
        //Do nothing with Gradient Detector now.
    }
}

void PureGeoStrategy::setPhysicalBoundaryConditions(SAMRAI::hier::Patch &patch,
                                                    const double fill_time,
                                                    const SAMRAI::hier::IntVector &ghost_width_tofill)
{
    //It has something to do with boudary utility, so I let it do nothing for now.
}

SAMRAI::hier::IntVector PureGeoStrategy::getRefineOpStencilWidth(const SAMRAI::tbox::Dimension &dim) const
{
    return SAMRAI::hier::IntVector::getZero(dim);
}

void PureGeoStrategy::preprocessRefine(
    SAMRAI::hier::Patch &fine,
    const SAMRAI::hier::Patch &coarse,
    const SAMRAI::hier::Box &fine_box,
    const SAMRAI::hier::IntVector &ratio)
{
    NULL_USE(fine);
    NULL_USE(coarse);
    NULL_USE(fine_box);
    NULL_USE(ratio);
}

void PureGeoStrategy::postprocessRefine(
    SAMRAI::hier::Patch &fine,
    const SAMRAI::hier::Patch &coarse,
    const SAMRAI::hier::Box &fine_box,
    const SAMRAI::hier::IntVector &ratio)
{
    NULL_USE(fine);
    NULL_USE(coarse);
    NULL_USE(fine_box);
    NULL_USE(ratio);
}

SAMRAI::hier::IntVector PureGeoStrategy::getCoarsenOpStencilWidth(
    const SAMRAI::tbox::Dimension &dim) const
{
    return SAMRAI::hier::IntVector::getZero(dim);
}

void PureGeoStrategy::preprocessCoarsen(
    SAMRAI::hier::Patch &coarse,
    const SAMRAI::hier::Patch &fine,
    const SAMRAI::hier::Box &coarse_box,
    const SAMRAI::hier::IntVector &ratio)
{
    NULL_USE(fine);
    NULL_USE(coarse);
    NULL_USE(coarse_box);
    NULL_USE(ratio);
}

void PureGeoStrategy::postprocessCoarsen(
    SAMRAI::hier::Patch &coarse,
    const SAMRAI::hier::Patch &fine,
    const SAMRAI::hier::Box &coarse_box,
    const SAMRAI::hier::IntVector &ratio)
{
    NULL_USE(fine);
    NULL_USE(coarse);
    NULL_USE(coarse_box);
    NULL_USE(ratio);
}
