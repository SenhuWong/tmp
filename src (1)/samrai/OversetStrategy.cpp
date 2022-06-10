#include "OversetStrategy.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "../mathTools/vector3d.h"
#include "PureGeoIntegrator.h"
#define SMALL_DISTURBTION 0.000000001

OversetStrategy::OversetStrategy(const std::string &object_name, const SAMRAI::tbox::Dimension &dim, std::shared_ptr<SAMRAI::geom::CartesianGridGeometry> grid_geometry)
    : d_object_name(object_name),
      d_dim(dim),
      d_grid_geometry(grid_geometry),
      d_allocator(SAMRAI::tbox::AllocatorDatabase::getDatabase()->getDefaultAllocator()),
      d_over_flag(new SAMRAI::pdat::CellVariable<int>(d_dim, "overFlag", d_allocator)),
      d_wall_flag(new SAMRAI::pdat::CellVariable<int>(d_dim, "wallFlag", d_allocator))
{
    std::cout << "ENtering getOverBoundary\n";
}

void OversetStrategy::registerFeeder(SamraiFeeder* fder)
{
    fder->getOverBoundary(nmeshes, &novers, &nwalls, &meshes_xyz, &meshes_over_ind, &meshes_over_ptr, &meshes_wall_ind, &meshes_wall_ptr);
}

OversetStrategy::OversetStrategy(const std::string &object_name, const SAMRAI::tbox::Dimension &dim, std::shared_ptr<SAMRAI::geom::CartesianGridGeometry> grid_geometry, SamraiFeeder *fder)
    : d_object_name(object_name),
      d_dim(dim),
      d_grid_geometry(grid_geometry),
      d_allocator(SAMRAI::tbox::AllocatorDatabase::getDatabase()->getDefaultAllocator()),
      d_over_flag(new SAMRAI::pdat::CellVariable<int>(d_dim, "overFlag", d_allocator)),
      d_wall_flag(new SAMRAI::pdat::CellVariable<int>(d_dim, "wallFlag", d_allocator))
{
    std::cout << "ENtering getOverBoundary\n";
    fder->getOverBoundary(nmeshes, &novers, &nwalls, &meshes_xyz, &meshes_over_ind, &meshes_over_ptr, &meshes_wall_ind, &meshes_wall_ptr);
}

OversetStrategy::~OversetStrategy()
{
}

void OversetStrategy::set_geometry(int i_nmesh, int *i_meshes_nover, double **i_meshes_xyzs, int **i_meshes_over_ind, int **i_meshes_over_ptr, int **i_meshes_wall_ind, int **i_meshes_wall_ptr)
{
    nmeshes = i_nmesh;
    novers = i_meshes_nover;
    meshes_xyz = i_meshes_xyzs;
    meshes_over_ind = i_meshes_over_ind;
    meshes_over_ptr = i_meshes_over_ptr;
    meshes_wall_ind = i_meshes_wall_ind;
    meshes_wall_ptr = i_meshes_wall_ptr;
}

#ifdef HAVE_HDF5
void OversetStrategy::registerVisItDataWriter(std::shared_ptr<SAMRAI::appu::VisItDataWriter> viz_writer)
{
    TBOX_ASSERT(viz_writer);
    d_visit_writer = viz_writer;
}
#endif

void OversetStrategy::registerTecplotWriter(std::shared_ptr<TecplotWriter> tec_writer)
{
    d_tecplot_writer = tec_writer;
}

void OversetStrategy::registerModelVariables(PureGeometricIntegrator *integrator)
{
    integrator->registerVariable(
        d_over_flag,
        SAMRAI::hier::IntVector::getZero(d_dim),
        d_grid_geometry,
        "These two doesn't",
        "Matter",
        PureGeometricIntegrator::VariableContextRegistry::FLAG_ONLY);
    integrator->registerVariable(
        d_wall_flag,
        SAMRAI::hier::IntVector::getZero(d_dim),
        d_grid_geometry,
        "Nothing really matters",
        "TO me",
        PureGeometricIntegrator::VariableContextRegistry::FLAG_ONLY);

    std::cout << "FIrst flag being registered\n";

    SAMRAI::hier::VariableDatabase *vardb = SAMRAI::hier::VariableDatabase::getDatabase();

    int overflag_id = vardb->mapVariableAndContextToIndex(d_over_flag, d_flag_contex);
    int wallflag_id = vardb->mapVariableAndContextToIndex(d_wall_flag, d_flag_contex);
    
    std::string dump_nameOver = "OVER #";
    std::string dump_nameWall = "WALL #";
    //Should I make the tecplot writer a Singleton?
    d_tecplot_writer->register_var(overflag_id);
    d_tecplot_writer->register_var(wallflag_id);

#ifdef HAVE_HDF5
    if (d_visit_writer)
    {
        d_visit_writer->registerPlotQuantity(dump_nameOver, "SCALAR", overflag_id, 0);
        d_visit_writer->registerPlotQuantity(dump_nameWall, "SCALAR", wallflag_id, 0);
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
void OversetStrategy::applyOversetBoundaryDetector(SAMRAI::hier::Patch &patch, const double init_data_time, const bool initial_time, const int tag_index)
{
    //%TODO::THis is done at last
    const std::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> patch_geom(
        SAMRAI_SHARED_PTR_CAST<SAMRAI::geom::CartesianPatchGeometry, SAMRAI::hier::PatchGeometry>(patch.getPatchGeometry()));
    TBOX_ASSERT(patch_geom);
    const double *dx = patch_geom->getDx();
    const double *xlo = patch_geom->getXLower();
    const double *xhi = patch_geom->getXUpper();

    std::shared_ptr<SAMRAI::pdat::CellData<int>> tags(
        SAMRAI_SHARED_PTR_CAST<SAMRAI::pdat::CellData<int>, SAMRAI::hier::PatchData>(patch.getPatchData(tag_index)));

    int over_index = SAMRAI::hier::VariableDatabase::getDatabase()->mapVariableAndContextToIndex(d_over_flag, d_flag_contex);
    std::shared_ptr<SAMRAI::pdat::CellData<int>> inside_over(
        SAMRAI_SHARED_PTR_CAST<SAMRAI::pdat::CellData<int>, SAMRAI::hier::PatchData>(patch.getPatchData(over_index)));
    int wall_index = SAMRAI::hier::VariableDatabase::getDatabase()->mapVariableAndContextToIndex(d_wall_flag, d_flag_contex);
    std::shared_ptr<SAMRAI::pdat::CellData<int>> inside_wall(
        SAMRAI_SHARED_PTR_CAST<SAMRAI::pdat::CellData<int>, SAMRAI::hier::PatchData>(patch.getPatchData(wall_index)));

    const SAMRAI::hier::Box &pbox = patch.getBox();
    const SAMRAI::hier::Index ifirst = pbox.lower();
    const SAMRAI::hier::Index ilast = pbox.upper();

    int *inside_over_ptr = inside_over->getPointer();
    int *inside_wall_ptr = inside_wall->getPointer();
    int *tag_ptr = tags->getPointer();
    int vol = 0;
    if (d_dim.getValue() == 2)
    {
        vol = (ilast(0) - ifirst(0) + 1) * (ilast(1) - ifirst(1) + 1);
    }
    else if (d_dim.getValue() == 3)
    {
        vol = (ilast(0) - ifirst(0) + 1) * (ilast(1) - ifirst(1) + 1) * (ilast(2) - ifirst(2));
    }
    for (int i = 0; i < vol; i++)
    {
        if (inside_wall_ptr[i] == SideTAG::insideTag or inside_over_ptr[i] == SideTAG::outsideTag)
        {
            tag_ptr[i] = false;
        }
        else
        {
            tag_ptr[i] = true;
        }
    }
}
void OversetStrategy::tagOverBoundaryCrossCells(SAMRAI::hier::Patch &patch, const double init_data_time)
{
    std::cout << "OversetStrategy::Tag over BOundary Cells being called\n";
    //The flag context should be registered with this overset strategy
    SAMRAI::hier::VariableDatabase *var_db = SAMRAI::hier::VariableDatabase::getDatabase();
    int over_flag_ind = var_db->mapVariableAndContextToIndex(d_over_flag, d_flag_contex);
    int wall_flag_ind = var_db->mapVariableAndContextToIndex(d_wall_flag, d_flag_contex);

    const std::shared_ptr<SAMRAI::pdat::CellData<int>> patch_overFlag(
        SAMRAI_SHARED_PTR_CAST<SAMRAI::pdat::CellData<int>, SAMRAI::hier::PatchData>(patch.getPatchData(over_flag_ind)));
    const std::shared_ptr<SAMRAI::pdat::CellData<int>> patch_wallFlag(
        SAMRAI_SHARED_PTR_CAST<SAMRAI::pdat::CellData<int>, SAMRAI::hier::PatchData>(patch.getPatchData(wall_flag_ind)));

    std::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> patch_geom(
        SAMRAI_SHARED_PTR_CAST<SAMRAI::geom::CartesianPatchGeometry, SAMRAI::hier::PatchGeometry>(patch.getPatchGeometry()));
    const double *xlo = patch_geom->getXLower();
    const double *dx = patch_geom->getDx();
    const double *xhi = patch_geom->getXUpper();

    SAMRAI::hier::Box pBox = patch.getBox();
    SAMRAI::hier::Index ifirst = pBox.lower();
    SAMRAI::hier::Index ilast = pBox.upper();
    tagBoundary2D(0, dx, xlo, xhi, ifirst(0), ilast(0), ifirst(1), ilast(1), patch_overFlag->getPointer());
    tagBoundary2D(1, dx, xlo, xhi, ifirst(0), ilast(0), ifirst(1), ilast(1), patch_wallFlag->getPointer());
    std::cout << "End of tagBoundary\n";
    FillBoundary2D(0, dx, xlo, xhi, ifirst(0), ilast(0), ifirst(1), ilast(1), patch_overFlag->getPointer());
    FillBoundary2D(1, dx, xlo, xhi, ifirst(0), ilast(0), ifirst(1), ilast(1), patch_wallFlag->getPointer());
    std::cout << "End of fillBoundary\n";
}

double getYatX(double *end_node, double *xloc)
{
    //get the value at xloc given 2 end point;
    double xA = end_node[0];
    double yA = end_node[1];
    double xB = end_node[2];
    double yB = end_node[3];
    double x = *xloc;
    double y = yA + (x - xA) * (yB - yA) / (xB - xA);
    return y;
}

void OversetStrategy::FillBoundary2D(int WallorOver, const double *dx, const double *xlo, const double *xhi, int ifirst0, int ilast0, int ifirst1, int ilast1, int *tag_inside_ptr)
{
    int *ncurs;
    int **cur_ind;
    int **cur_ptr;
    if (WallorOver == 1)
    {
        ncurs = nwalls;
        cur_ind = meshes_wall_ind;
        cur_ptr = meshes_wall_ptr;
    }
    else if (WallorOver == 0)
    {
        ncurs = novers;
        cur_ind = meshes_over_ind;
        cur_ptr = meshes_over_ptr;
    }
    
    //Use the standard way
    int strides[2] = {ilast0 - ifirst0 + 1, ilast1 - ifirst1 + 1};
    int infinity = 100000;
    int vol = (strides[0] + 1) * (strides[1] + 1);
    int *nodeFlag = new int[vol];
    for (int i = 0; i < vol; i++)
    {
        nodeFlag[i] = 0; //Initial value set to 0;
    }
    bool find_inside = false;
    bool find_outside = false;

    //For every box,check the 4 vert with standard method
    for (int i = 0; i < strides[0]; i++)
    {
        for (int j = 0; j < strides[1]; j++)
        {
            find_inside = false;
            find_outside = false;
            if (tag_inside_ptr[strides[0] * j + i] == SideTAG::crossedTag)
            {
                continue;
            }
            //If find a node inside ,then set all others inside
            int inside_count;
            int outside_count;
            for (int k = 0; k < 4; k++)
            {
                int i_incre = (k & 1);
                int j_incre = (k & 2) >> 1;
                int nodei = i + i_incre;
                int nodej = j + j_incre;
                if (nodeFlag[(strides[0] + 1) * nodej + nodei] == -SideTAG::outsideTag) //-SideTAG::outsideTag means that this node has been tagged as outside
                {
                    outside_count++;
                    tag_inside_ptr[strides[0] * j + i] = SideTAG::outsideTag;
                    find_outside = true;
                }
                else if (nodeFlag[(strides[0] + 1) * nodej + nodei] == -SideTAG::insideTag) //1 means that this node has been tagged as inside
                {
                    inside_count++;
                    tag_inside_ptr[strides[0] * j + i] = SideTAG::insideTag;
                    find_inside = true;
                }
            }
            if (find_inside and find_outside)
            {
                //Two things happened at the same time, suggesting there is a node right on the line
                std::cout << "Can't be inside and outside at the same time\n";
            }
            if (!find_inside and !find_outside) //Not inside and outside,we need to calculate
            {
                //std::cout<<"Need to be specified\n";
                Vector3D this_pt = Vector3D(xlo[0] + dx[0] * i, xlo[1] + dx[1] * j, 0);
                Vector3D far_left = this_pt - Vector3D(infinity, 0, 0);
                Vector3D hor_line = far_left - this_pt;
                int cross_count = 0;
                for (int k = 0; k < nmeshes; k++)
                {
                    for (int j = 0; j < ncurs[k]; j++)
                    {
                        int inodeA = cur_ind[k][2 * j] - 1;
                        int inodeB = cur_ind[k][2 * j + 1] - 1;
                        Vector3D ptA = Vector3D(meshes_xyz[k][2 * inodeA], meshes_xyz[k][2 * inodeA + 1], 0);
                        Vector3D ptB = Vector3D(meshes_xyz[k][2 * inodeB], meshes_xyz[k][2 * inodeB + 1], 0);

                        Vector3D AB = (ptB - ptA);
                        double result1 = (AB.cross_product(far_left - ptA)).dot_product(AB.cross_product(this_pt - ptA));
                        double result2 = (hor_line.cross_product(ptA - this_pt)).dot_product(hor_line.cross_product(ptB - this_pt));
                        if (result1 <= 0 and result2 <= 0)
                        {
                            cross_count++;
                        }
                    }
                }
                if (cross_count % 2 == 0) //Outside
                {
                    find_outside = true;
                    tag_inside_ptr[strides[0] * j + i] = SideTAG::outsideTag;
                }
                else if (cross_count % 2 == 1)
                {
                    find_inside = true;
                    tag_inside_ptr[strides[0] * j + i] = SideTAG::insideTag;
                }
            }
            for (int k = 0; k < 4; k++)
            {

                int i_incre = (k & 1);
                int j_incre = (k & 2) >> 1;
                int nodei = i + i_incre;
                int nodej = j + j_incre;
                if (find_outside) //2 means that this node has been tagged as outside
                {
                    nodeFlag[(strides[0] + 1) * nodej + nodei] = -SideTAG::outsideTag;
                }
                else if (find_inside) //1 means that this node has been tagged as inside
                {
                    nodeFlag[(strides[0] + 1) * nodej + nodei] = -SideTAG::insideTag;
                }
            }
        }
    }
    delete[] nodeFlag;
}

void OversetStrategy::tagBoundary2D(int WallorOver, const double *dx, const double *xlo, const double *xhi, int ifirst0, int ilast0, int ifirst1, int ilast1, int *tag_ptr)
{
    int **cur_ind;
    int **cur_ptr;
    int *ncurs;
    if (WallorOver == 1)
    {
        ncurs = nwalls;
        cur_ind = meshes_wall_ind;
        cur_ptr = meshes_wall_ptr;
    }
    else if (WallorOver == 0)
    {
        ncurs = novers;
        cur_ind = meshes_over_ind;
        cur_ptr = meshes_over_ptr;
    }
    const int BIGVALUE = 198964;
    bool intersect = true;
    int dim = d_dim.getValue();
    double x_nodes[12]; //4*3 or 2*2
    int count = 0;
    double xmin[3];
    double xmax[3]; //For line segment
    int strides[2] = {ilast0 - ifirst0 + 1, ilast1 - ifirst1 + 1};
    int imin[3];
    int imax[3];
    //return;
    //std::cout<<nmeshes<<" nmeshes\n";
    for (int i = 0; i < nmeshes; i++)
    {
        //For this mesh
        for (int j = 0; j < ncurs[i]; j++)
        {
            //For this edge
            //First get the bounding box to see if there is intersection
            xmin[0] = xmin[1] = xmin[2] = BIGVALUE;
            xmax[0] = xmax[1] = xmax[2] = -BIGVALUE;
            count = 0;
            for (int k = cur_ptr[i][j]; k < cur_ptr[i][j + 1]; k++)
            {
                int cur_node = cur_ind[i][k] - 1;
                for (int l = 0; l < dim; l++)
                {
                    xmin[l] = std::min(meshes_xyz[i][dim * cur_node + l], xmin[l]);
                    xmax[l] = std::max(meshes_xyz[i][dim * cur_node + l], xmax[l]);
                    x_nodes[count++] = meshes_xyz[i][dim * cur_node + l];
                }
            }
            intersect = true;
            for (int k = 0; k < dim; k++)
            {
                imin[k] = (xmin[k] - SMALL_DISTURBTION - xlo[k]) / dx[k]; //This may be negative but doesn't matter
                imax[k] = (xmax[k] + SMALL_DISTURBTION - xlo[k]) / dx[k];
                //std::cout<<"imax - imin "<<imax[k]-imin[k]<<'\n';
                if (xlo[k] > xmax[k] or xhi[k] < xmin[k])
                {
                    intersect = false;
                    //Bounding box doesn't overlap
                }
            }
            if (intersect)
            {
                //In 2D a simple
                imin[0] = std::max(imin[0], 0);                //Starting from 0
                imax[0] = std::min(imax[0], ilast0 - ifirst0); //Ending at ilast - ifirst
                imin[1] = std::max(imin[1], 0);
                imax[1] = std::min(imax[1], ilast1 - ifirst1);
                double xloc = xlo[0] + imin[0] * dx[0];
                double yloc_at_xloc;
                double xloc_plus;
                double yloc_at_xloc_plus;
                int tag_imin;
                int tag_imax;
                for (int k = imin[0]; k <= imax[0]; k++) //Along x
                {

                    xloc_plus = xloc + dx[0];
                    yloc_at_xloc = getYatX(x_nodes, &xloc);
                    yloc_at_xloc_plus = getYatX(x_nodes, &xloc_plus);
                    //Do something
                    int i1_min = (yloc_at_xloc - SMALL_DISTURBTION - xlo[1]) / dx[1];
                    int i2_min = (yloc_at_xloc_plus - SMALL_DISTURBTION - xlo[1]) / dx[1];
                    int i1_max = (yloc_at_xloc + SMALL_DISTURBTION - xlo[1]) / dx[1];
                    int i2_max = (yloc_at_xloc_plus + SMALL_DISTURBTION - xlo[1]) / dx[1];
                    tag_imin = std::min(i1_min, i2_min);
                    tag_imin = std::max(tag_imin, imin[1]);
                    tag_imax = std::max(i1_max, i2_max);
                    tag_imax = std::min(tag_imax, imax[1]);
                    if (imin[0] == imax[0])
                    {
                        //std::cout << "Inside one" << imax[0] << '\t' << tag_imin << '\t' << tag_imax << '\n';
                    }

                    if (tag_imax < tag_imin)
                    {
                        //Don't do that
                    }
                    else
                    {
                        for (int l = tag_imin; l <= tag_imax; l++)
                        {
                            //std::cout << "Tagged" << k << ',' << l << "\n";
                            if (strides[0] * l + k >= (strides[0] + 1) * (strides[1] + 1) or strides[0] * l + k < 0)
                            {
                                //std::cout << strides[0] * l + k << " out of " << strides[0] << " * " << strides[1] << "=" << strides[0] * strides[1] << " " << l << " " << k << '\n';
                            }
                            tag_ptr[strides[0] * l + k] = SideTAG::crossedTag;
                        }
                        //Tag all the value in between
                    }
                    xloc = xloc_plus;
                }
            }
        }
    }
}

void tagOversetBoundary3D(const double *dx, const double *xlo, const double *xhi, int ifirst0, int ilast0, int ifirst1, int ilast1, int ifirst2, int ilast2, int *tag_inside_ptr)
{
    //For each face, find its bounding box and compare with the current Patch
    //If overlap exists, for every cell inside the patch, do a seperating axis checking
}
void eleminateOversetBoundary3D(const double *dx, const double *xlo, const double *xhi, int ifirst0, int ilast0, int ifirst1, int ilast1, int ifirst2, int ilast2, int *tag_inside_ptr)
{
    //Do a horizontal line ,For each face check overlap ,isf inside tag inside
    //If a cell has one node tagged as inside and not a crossed one, then tag as inside with all nodes tagged as inside
}
