#include "Reader.h"

Reader::Reader()
{

}
Reader::~Reader()
{
    
}

void Reader::addFile(const std::string &filename)
{
    d_filenames.push_back(filename);
}

void Reader::setComm(int icur_proc, int inum_proc)
{
    cur_proc = icur_proc;
    num_proc = inum_proc;
}
void Reader::setComm(MPI_Comm &scomm)
{
    MPI_Comm_rank(scomm, &cur_proc);
    MPI_Comm_size(scomm, &num_proc);
}