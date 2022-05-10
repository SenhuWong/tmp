#pragma once
//This class now handles all that is needed to refine Boxes crossed by overset grid
//Should this class handle the interpolation of patch data with tioga?
//The output of all patch to tioga is done at integrator level(or this one if I like)
//The interpolation could also be done at this one.
#include <string>
#include "SAMRAI/tbox/Dimension.h"
#include "SamraiFeeder.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/appu/VisItDataWriter.h"
#include "tecplotWriter.h"
class PureGeometricIntegrator;
//OversetStrategy for now is just for postprocessing
class OversetStrategy
{
    enum SideTAG
    {
        insideTag = 10,crossedTag,outsideTag
    };
private:
    //This strategy holds a collection of lines in 2D here
    int nmeshes = 0;
    //These are maintained by SamraiFeeder Class and DON'T NEED TO FREE IN THIS CLASS;
    int *novers = NULL;
    int *nwalls = NULL;
    double **meshes_xyz = NULL;
    int **meshes_over_ind = NULL;
    int **meshes_over_ptr = NULL;
    int **meshes_wall_ind = NULL;
    int **meshes_wall_ptr = NULL;
    std::string d_object_name;
    SAMRAI::tbox::Dimension d_dim;
    //Flag Variable
    std::shared_ptr<SAMRAI::pdat::CellVariable<int>> d_over_flag;
    std::shared_ptr<SAMRAI::pdat::CellVariable<int>> d_wall_flag;

    //std::shared_ptr<SAMRAI::pdat::CellVariable<int>> d_interp_flag;
    //Flag Context
    std::shared_ptr<SAMRAI::hier::VariableContext> d_flag_contex;
    //Cartesian Geometry
    std::shared_ptr<SAMRAI::geom::CartesianGridGeometry> d_grid_geometry;
    SAMRAI::tbox::ResourceAllocator d_allocator;

public:
    OversetStrategy(const std::string &object_name, const SAMRAI::tbox::Dimension &dim, std::shared_ptr<SAMRAI::geom::CartesianGridGeometry> grid_geometry, SamraiFeeder *fder);
    OversetStrategy(const std::string &object_name, const SAMRAI::tbox::Dimension &dim, std::shared_ptr<SAMRAI::geom::CartesianGridGeometry> grid_geometry);
    void registerFeeder(SamraiFeeder *fder);
    ~OversetStrategy();
    void set_Context(std::shared_ptr<SAMRAI::hier::VariableContext> flag_contxt)
    {
        d_flag_contex = flag_contxt;
    }
    void set_geometry(int i_nmesh, int *i_meshes_nover, double **i_meshes_xyzs, int **i_meshes_over_ind, int **i_meshes_over_ptr, int **i_meshes_wall_ind, int ** imeshes_wall_ptr);
    void registerModelVariables(PureGeometricIntegrator *integrator);
    void applyOversetBoundaryDetector(SAMRAI::hier::Patch &patch, const double init_data_time, const bool initial_time, const int tag_index);
    void tagOverBoundaryCrossCells(SAMRAI::hier::Patch &patch, const double init_data_time);
    //This is for tagging 2D oversetBoundary
    void tagBoundary2D(int WallorOver,const double *dx, const double *xlo, const double *xhi, int ifirst0, int ilast0, int ifirst1, int ilast1, int *tag_inside_ptr);
    void FillBoundary2D(int WallorOver,const double *dx, const double *xlo, const double *xhi, int ifirst0, int ilast0, int ifirst1, int ilast1, int *tag_inside_ptr);
    //This is for tagging 3D oversetBoundary
    void makeTriangles();
    void tagOversetBoundary3D(const double *dx, const double *xlo, const double *xhi, int ifirst0, int ilast0, int ifirst1, int ilast1, int ifirst2, int ilast2, int *tag_inside_ptr);
    void eleminateOversetBoundary3D(const double *dx, const double *xlo, const double *xhi, int ifirst0, int ilast0, int ifirst1, int ilast1, int ifirst2, int ilast2, int *tag_inside_ptr);

    void checking(
        SAMRAI::hier::Patch &patch,
        int level,
        int max_level)
    {
        SAMRAI::hier::VariableDatabase *var_db = SAMRAI::hier::VariableDatabase::getDatabase();
        //%TODO::This is needed to be reused so that the output is good
        int index_flag = var_db->mapVariableAndContextToIndex(d_over_flag, d_flag_contex);

        const SAMRAI::hier::Box &pbox = patch.getBox();
        const SAMRAI::hier::Index ifirst = pbox.lower();
        const SAMRAI::hier::Index ilast = pbox.upper();

        int strides[2] = {ilast(0) - ifirst(0) + 1, ilast(1) - ifirst(1) + 1};

        std::shared_ptr<SAMRAI::pdat::CellData<int>> flag_data =
            SAMRAI_SHARED_PTR_CAST<SAMRAI::pdat::CellData<int>, SAMRAI::hier::PatchData>(patch.getPatchData(index_flag));
        int *ptr = flag_data->getPointer();

        for (int j = 0; j < strides[1]; j++)
        {
            for (int i = 0; i < strides[0]; i++)
            {
                int index = strides[0] * j + i;
                if (ptr[index] == -1)
                {
                    ptr[index] = 1;
                    std::cout << "Some data on level " << level << "(" << i << "," << j << ") out of " << max_level << " is -1\n";
                }
            }
        }
    }
#ifdef HAVE_HDF5
    void registerVisItDataWriter(std::shared_ptr<SAMRAI::appu::VisItDataWriter> viz_writer);
    std::shared_ptr<SAMRAI::appu::VisItDataWriter> d_visit_writer;
#endif
    void registerTecplotWriter(std::shared_ptr<TecplotWriter> tec_writer);
    std::shared_ptr<TecplotWriter> d_tecplot_writer;
};