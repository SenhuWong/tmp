//This is an abstract base class providing all interfaces for SAMRAI and TIOGA requests
#pragma once
#include <string>
#include <vector>
#include <mpi.h>

class Reader
{
public:
    std::vector<std::string> d_filenames;
    int d_nmesh = -1;
    int d_dim = -1;
    int num_proc = -1;
    int cur_proc = -1;
    Reader();
    virtual ~Reader();
    void addFile(const std::string &filename);
    void setComm(int icur_proc, int inum_proc);
    void setComm(MPI_Comm& scomm);
    //Not all format need a flag telling if it's wall boundary or over boundary 
    virtual void setDim(int idim) = 0;//There might be need for Brick's dimension setting
    virtual void readAll() = 0;

};
