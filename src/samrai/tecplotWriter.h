#pragma once
#include "SAMRAI/hier/PatchHierarchy.h"
#include <string>
#include <vector>
#include <fstream>
#include "SAMRAI/tbox/Dimension.h"
//Later we might need a enum to tell the type of data that we load
//Or template maybe???
class TecplotWriter
{
    //We might need a struct to store registered variables
private:
    std::string filename;
    SAMRAI::tbox::Dimension d_dim;
    std::shared_ptr<SAMRAI::hier::PatchHierarchy> d_hierarchy;
    std::vector<int>register_variable;
    std::ofstream fout;
public:
    TecplotWriter(const std::string& ifilename, const SAMRAI::tbox::Dimension& idim);
    ~TecplotWriter();
    void init();
    void finalize();
    void set_filename(const std::string& filename);
    void register_var(int var_ind);
    void attachHierarchy(const std::shared_ptr<SAMRAI::hier::PatchHierarchy> &hier);
    void writeAll();
    void writeLevel(int ln);
    void writePatch(const SAMRAI::hier::Patch& patch,int ln, int pn);
    void writePatch2D(const SAMRAI::hier::Patch&patch,int ln,int pn);
    int size()
    {
        return register_variable.size();
    }
    void writeHeader();
    //So the regular routine is to:
    //Constructor()
    //set_filename()
    //register_variable()
    //attachHierarchy()
    //writeAll()
    
};

