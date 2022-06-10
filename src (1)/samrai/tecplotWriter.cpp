#include "tecplotWriter.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
TecplotWriter::TecplotWriter(const std::string &ifilename, const SAMRAI::tbox::Dimension &idim)
    : filename(ifilename), d_dim(idim)
{
    const SAMRAI::tbox::SAMRAI_MPI &mpi(SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld());
    int cur_proc = mpi.getRank();
    filename = ifilename + std::to_string(cur_proc);
}
TecplotWriter::~TecplotWriter()
{
}

void TecplotWriter::init()
{
    fout.open(filename);
}

void TecplotWriter::finalize()
{
    std::cout<<"Finalizing tecplotwriter\n";
    fout.close();
}

void TecplotWriter::register_var(int ind)
{
    //std::cout<<"register_var called\n\n\n\n\n\n";
    register_variable.push_back(ind);
}

void TecplotWriter::set_filename(const std::string &ifilename)
{
    filename = ifilename;
}
void TecplotWriter::attachHierarchy(const std::shared_ptr<SAMRAI::hier::PatchHierarchy> &hier)
{
    d_hierarchy = hier;
}
void TecplotWriter::writeAll()
{
    //First we might want to write the Header file
    //write_header()
    //Then we might want a static int to tell counting
    writeHeader();
    int num_lvl = d_hierarchy->getNumberOfLevels();
    for (int i = 0; i < num_lvl; i++)
    {
        writeLevel(i);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void TecplotWriter::writeLevel(int ln)
{
    std::shared_ptr<SAMRAI::hier::PatchLevel> level(
        d_hierarchy->getPatchLevel(ln));
    int count = 0;
    for (SAMRAI::hier::PatchLevel::iterator ip(level->begin()); ip != level->end(); ip++)
    {
        const std::shared_ptr<SAMRAI::hier::Patch> &patch = *ip;
        writePatch2D(*patch, ln, count);
        count++;
    }
}

void TecplotWriter::writeHeader()
{
    //For now let's require all variables are cell_centered
    int nregistered = register_variable.size();
    int dim = d_dim.getValue();
    fout << "Title = " << filename << '\n';
    fout << "VARIABLES = \"X\", \"Y\"";
    if (dim == 3)
    {
        fout << ", \"Z\"";
    }
    SAMRAI::hier::VariableDatabase *var_db = SAMRAI::hier::VariableDatabase::getDatabase();
    for (auto iter = register_variable.begin(); iter != register_variable.end(); iter++)
    {
        std::shared_ptr<SAMRAI::hier::Variable> cur_variable;
        std::shared_ptr<SAMRAI::hier::VariableContext> cur_context;
        var_db->mapIndexToVariableAndContext(*iter, cur_variable, cur_context);
        std::string var_name = cur_variable->getName() + "_" + cur_context->getName();
        std::cout<<var_name<<'\n';
        //std::cout<<var_name<<'\n';
        fout << ", \"" << var_name << "\"";
    }
    fout << '\n';
}

void TecplotWriter::writePatch2D(const SAMRAI::hier::Patch &patch, int ln, int pn)
{
    const std::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> patch_geom(
        SAMRAI_SHARED_PTR_CAST<SAMRAI::geom::CartesianPatchGeometry, SAMRAI::hier::PatchGeometry>(patch.getPatchGeometry()));
    const double *dx = patch_geom->getDx();
    const double *xlo = patch_geom->getXLower();
    const SAMRAI::hier::Box &pbox = patch.getBox();
    const SAMRAI::hier::Index ifirst = pbox.lower();
    const SAMRAI::hier::Index ilast = pbox.upper();

    int strides[2] = {ilast(0) - ifirst(0) +1, ilast(1) - ifirst(1)+1};
    int nnodes = (strides[0] + 1) * (strides[1] + 1);
    int nelems = strides[0] * strides[1];
    fout << "ZONE T = \"" << ln << "," << pn << "\", N=" << nnodes << ", E=" << nelems << ", DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL\n";
    fout << "VARLOCATION=(1=NODAL, 2=NODAL";
    for (int i = 0; i < register_variable.size(); i++)
    {
        fout << ", " << 3 + i << "=CELLCENTERED";
    }
    fout << ")\n";

    int count = 0;
    for (int j = 0; j < strides[1] + 1; j++)
    {
        for (int i = 0; i < strides[0] + 1; i++)
        {
            fout << xlo[0] + dx[0] * i << '\n';
        }
    }
    for (int j = 0; j < strides[1] + 1; j++)
    {
        for (int i = 0; i < strides[0] + 1; i++)
        {
            fout << xlo[1] + dx[1] * j << '\n';
        }
    }
    

    for (auto iter = register_variable.begin(); iter != register_variable.end(); iter++)
    {
        int ind = *iter;
        patch.getPatchData(ind);
        //TODO::How does this be solved?
        std::shared_ptr<SAMRAI::pdat::CellData<int>> tags(
            SAMRAI_SHARED_PTR_CAST<SAMRAI::pdat::CellData<int>, SAMRAI::hier::PatchData>(patch.getPatchData(ind)));
        int *cur_data = tags->getPointer();
        for (int j = 0; j < strides[1]; j++)
        {
            for (int i = 0; i < strides[0]; i++)
            {
                fout << cur_data[strides[0] * j + i] << '\n';
            }
        }
    }
    //Now perform the connectivity
    for (int j = 0; j < strides[1]; j++)
    {
        for (int i = 0; i < strides[0]; i++)
        {
            fout << j * (strides[0] + 1) + i + 1 << " " << j * (strides[0] + 1) + i + 2 << " " << (j + 1) * (strides[0] + 1) + i + 2 << " " << (j + 1) * (strides[0] + 1) + i + 1 << '\n';
        }
    }
    return;
}

void TecplotWriter::writePatch(const SAMRAI::hier::Patch &patch, int ln, int pn)
{
    std::cout << "Writing " << ln << ',' << pn << '\n';
    fout << "ZONE T=\"" << ln << "_" << pn << "\",";
    const std::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> patch_geom(
        SAMRAI_SHARED_PTR_CAST<SAMRAI::geom::CartesianPatchGeometry, SAMRAI::hier::PatchGeometry>(patch.getPatchGeometry()));
    const double *dx = patch_geom->getDx();
    const double *xlo = patch_geom->getXLower();
    const SAMRAI::hier::Box &pbox = patch.getBox();
    const SAMRAI::hier::Index ifirst = pbox.lower();
    const SAMRAI::hier::Index ilast = pbox.upper();
    int strides[3] = {
        ilast(0) - ifirst(0) + 1,
        ilast(1) - ifirst(1) + 1,
        d_dim.getValue() == 3 ? ilast(2) - ifirst(2) + 1 : 0};
    if (d_dim.getValue() == 2)
    {
        fout << "I=" << strides[0] + 1 << ",J=" << strides[1] + 1 << ",DATAPACKING=BLOCK";
        fout << ",VARLOCATION=([" << 3 << "," << 2 + register_variable.size() << "]=CELLCENTERED)\n";
    }
    else if (d_dim.getValue() == 3)
    {
        fout << "I=" << strides[0] + 1 << ",J=" << strides[1] + 1 << ",K=" << strides[2] + 1 << ",DATAPACKING=BLOCK";
        fout << ",VARLOCATION=([" << 4 << "," << 3 + register_variable.size() << "]=CELLCENTERED)\n";
    }
    int count = 0;
    //X alongf [x,y,z]
    if (d_dim.getValue() == 2)
    {
        count = 0;
        for (int i = 0; i < strides[0] + 1; i++)
        {
            for (int j = 0; j < strides[1] + 1; j++)
            {
                fout << xlo[0] + i * dx[0] << " ";
                count++;
                if (count % 10 == 0)
                {
                    fout << '\n';
                }
            }
        }
        if (count % 10)
        {
            fout << '\n';
        }
        count = 0;
        for (int i = 0; i < strides[0] + 1; i++)
        {
            for (int j = 0; j < strides[1] + 1; j++)
            {
                fout << xlo[1] + j * dx[1] << " ";
                count++;
                if (count % 10 == 0)
                {
                    fout << '\n';
                }
            }
        }
        if (count % 10)
        {
            fout << '\n';
        }
    }
    else if (d_dim.getValue() == 3)
    {
        count = 0;
        for (int i = 0; i < strides[0] + 1; i++)
        {
            for (int j = 0; j < strides[1] + 1; j++)
            {
                for (int k = 0; k < strides[2] + 1; k++)
                {
                    fout << xlo[0] + i * dx[0] << " ";
                    count++;
                    if (count % 10 == 0)
                    {
                        fout << '\n';
                    }
                }
            }
        }
        if (count % 10)
        {
            fout << '\n';
        }
        count = 0;
        for (int i = 0; i < strides[0] + 1; i++)
        {
            for (int j = 0; j < strides[1] + 1; j++)
            {
                for (int k = 0; k < strides[2] + 1; k++)
                {
                    fout << xlo[1] + j * dx[1] << " ";
                    count++;
                    if (count % 10 == 0)
                    {
                        fout << '\n';
                    }
                }
            }
        }
        if (count % 10)
        {
            fout << '\n';
        }
        count = 0;
        for (int i = 0; i < strides[0] + 1; i++)
        {
            for (int j = 0; j < strides[1] + 1; j++)
            {
                for (int k = 0; k < strides[2] + 1; k++)
                {
                    fout << xlo[2] + k * dx[2] << " ";
                    count++;
                    if (count % 10 == 0)
                    {
                        fout << '\n';
                    }
                }
            }
        }
        if (count % 10)
        {
            fout << '\n';
        }
    }
    int x_stride = 1;
    int y_stride = x_stride * strides[0];
    int z_stride = y_stride * strides[1];
    for (auto iter = register_variable.begin(); iter != register_variable.end(); iter++)
    {
        int ind = *iter;
        patch.getPatchData(ind);
        //TODO::How does this be solved?
        std::shared_ptr<SAMRAI::pdat::CellData<int>> tags(
            SAMRAI_SHARED_PTR_CAST<SAMRAI::pdat::CellData<int>, SAMRAI::hier::PatchData>(patch.getPatchData(ind)));
        if (d_dim.getValue() == 2)
        {
            count = 0;
            for (int i = 0; i < strides[0]; i++)
            {
                for (int j = 0; j < strides[1]; j++)
                {
                    fout << (tags->getPointer())[i * x_stride + j * y_stride] << " ";
                    count++;
                    if (count % 10 == 0)
                    {
                        fout << '\n';
                    }
                }
            }
            if (count % 10)
            {
                fout << '\n';
            }
        }
        else if (d_dim.getValue() == 3)
        {
            count = 0;
            for (int i = 0; i < strides[0]; i++)
            {
                for (int j = 0; j < strides[1]; j++)
                {
                    for (int k = 0; k < strides[2]; k++)
                    {
                        fout << (tags->getPointer())[i * x_stride + j * y_stride + k * z_stride] << " ";
                        count++;
                        if (count % 10 == 0)
                        {
                            fout << '\n';
                        }
                    }
                }
            }
            if (count % 10)
            {
                fout << '\n';
            }
        }
        fout << '\n';
    }
}