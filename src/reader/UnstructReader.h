#pragma once
//This class is a CobaltReader for UnstructSolver test

#include "../UnstructSolver/UnstructFeeder.h"
#include "../UniTioga/CellBrick.h"
class UnstructReader :public UnstructFeeder
{

    const static int cellBase = 1;
    const static int nodeBase = 1;

    int signWall = -2;
    int signOver = -1;

public:
    UnstructReader()
        :UnstructFeeder()
    {

    }

    void takeBoundary(int sign,int** nEach,int*** indEach)
    {
        *nEach = new int[d_nmesh];
        *indEach = new int*[d_nmesh];
        std::vector<int> ind;
        if(sign==signWall)
        {
            if(d_dim==2)
            {
                for(int i = 0;i<d_nmesh;i++)
                {
                    ind.clear();
                    auto& curBlk = UBs2[i];
                    for(int j = 0;j<curBlk.nEdges();j++)
                    {
                        auto& curEdge = curBlk.d_localEdges[j];
                        if(curEdge.rCInd()==GeomElements::edge3d<2>::BoundaryType::WALL)
                        {
                            ind.push_back(curEdge.lCInd());
                        }
                    }
                    (*nEach)[i] = ind.size();
                    (*indEach)[i] = new int[ind.size()];
                    int temp_int = 0;
                    for(auto iter = ind.begin();iter!=ind.end();iter++)
                    {
                        (*indEach)[i][temp_int++] = *iter;
                    }
                }

            }
            


        }
        else if(sign==signOver)
        {
            if(d_dim==2)
            {
                for(int i = 0;i<d_nmesh;i++)
                {
                    ind.clear();
                    auto& curBlk = UBs2[i];
                    for(int j = 0;j<curBlk.nEdges();j++)
                    {
                        auto& curEdge = curBlk.d_localEdges[j];
                        if(curEdge.rCInd()==GeomElements::edge3d<2>::BoundaryType::FARFIELD)
                        {
                            ind.push_back(curEdge.lCInd());
                        }
                    }
                    (*nEach)[i] = ind.size();
                    (*indEach)[i] = new int[ind.size()];
                    int temp_int = 0;
                    for(auto iter = ind.begin();iter!=ind.end();iter++)
                    {
                        (*indEach)[i][temp_int++] = *iter;
                    }
                }

            }

        }

    }

    void readAll()
    {
        initializeBlockInfo();
        readForUnstruct();
    }

    void setDim(int ndim)
    {
        d_dim= ndim;
    }

    void setOverSign(int sign)
    {
        signOver = sign;
    }

    void setWallSign(int sign)
    {
        signWall = sign;
    }

    void initializeBlockInfo()
    {
        d_nmesh = d_filenames.size();
        UnstructFeeder::initialize();
    }

    void writeFiles()
    {
        for (int i = 0; i < d_nmesh; i++)
	    {
		    UBs2[i].write_grd(d_filenames[i] + "SerialUnstruct", i + 1, cur_proc);
		    UBs2[i].write_grdSendRecv(d_filenames[i] + "SerialUnstructSRFlags", i + 1, cur_proc);
	    }
    }

    void readForUnstruct()
    {
        for(int m = 0;m<d_nmesh;m++)
        {
            std::string filename = d_filenames[m];
            std::ifstream fin;
            fin.open(filename);
            if(!fin.is_open())
            {
                throw std::runtime_error("Open file Failure\n");
            }
            UnstructBlock2D<2>* cur_ub2D=NULL;
            UnstructBlock2D<3>* cur_ub3D=NULL;
            if(d_dim==2)
            {
                cur_ub2D = &(UBs2[m]);
            }
            else if(d_dim==3)
            {
                cur_ub3D = &(UBs3[m]);
            }

            GeomElements::point3d<2>* pt2d;
            GeomElements::point3d<3>* pt3d;
            GeomElements::edge3d<2>* edg2d;
            GeomElements::edge3d<3>* edg3d;
            GeomElements::cell3d<2>* cel2d;
            GeomElements::cell3d<3>* cel3d;

            int nnode;
            int nedge;
            int ncell;
            int to_throw_away;
            fin >> to_throw_away >> to_throw_away >> to_throw_away;
			fin >> nnode >> nedge >> ncell >> to_throw_away >> to_throw_away;
            double* xyzs = new double[d_dim*nnode];
            for(int i = 0;i<nnode;i++)
            {
                for(int j = 0;j<d_dim;j++)
                {
                    fin >>xyzs[d_dim*i+j];
                }
            }

            Brick4Cobalt2D* cells2D = nullptr;
            Brick4Cobalt3D* cells3D = nullptr;
            
            int nodeCount;
            int edgeNodes[8];
            int *edge2Cell = new int[2*nedge];

            std::vector<int> edgeNodePtr(nedge + 1, -1);
            std::vector<int> edgeNodeInd;
            edgeNodeInd.reserve(2 * (d_dim - 1) * nedge);
            edgeNodePtr[0] = 0;

            if(d_dim==2)
            {
                cells2D = new Brick4Cobalt2D[ncell];
                for(int i = 0;i<nedge;i++)
                {
                    fin >> nodeCount;
                    edgeNodePtr[i+1] = edgeNodePtr[i] + nodeCount;
                    for(int j = 0;j<nodeCount;j++)
                    {
                        fin>>edgeNodes[j];
                        edgeNodeInd.push_back(edgeNodes[j]);
                    }

                    fin >>edge2Cell[2*i]>>edge2Cell[2*i+1];
                    cells2D[edge2Cell[2*i]-cellBase].insert(edgeNodes,xyzs,nodeCount);
                    if(edge2Cell[2*i+1]>0)
                    {
                        cells2D[edge2Cell[2*i+1]-cellBase].insert(edgeNodes,xyzs,nodeCount);
                    }
                }
            }
            else
            {
                throw std::runtime_error("Not implemented yet\n");
            }

            if(d_dim==2)
            {
                pt2d = new GeomElements::point3d<2>[nnode];
                for(int i = 0;i<nnode;i++)
                {
                    pt2d[i] = GeomElements::vector3d<2, double>(xyzs[2*i],xyzs[2*i+1]);
                }
            }

            if(d_dim==2)
            {
                cel2d = new GeomElements::cell3d<2>[ncell];
                for(int i = 0;i<ncell;i++)
                {
                    int size = cells2D[i].sizeIs();
                    for(int j = 0;j<size;j++)
                    {
                        cel2d[i].push_back(cells2D[i].nodeIs(j)-nodeBase);
                    }
                }
            }

            if(d_dim==2)
            {
                edg2d = new GeomElements::edge3d<2>[nedge];
                for(int i = 0;i<nedge;i++)
                {
                    for(int j = edgeNodePtr[i];j<edgeNodePtr[i+1];j++)
                    {
                        edg2d[i].push_back(edgeNodeInd[j]-nodeBase);
                    }
                    int lC = edge2Cell[2*i];
                    int rC = edge2Cell[2*i+1];

                    if(lC>0)
                    {
                        edg2d[i].setLeft(lC-cellBase);
                        cel2d[lC-cellBase].push_edge(i);
                    }
                    else
                    {
                        throw std::runtime_error("Left cell can't be negative\n");
                    }
                    if(rC>0)
                    {
                        edg2d[i].setRight(rC-cellBase);
                        cel2d[rC-cellBase].push_edge(i);
                    }
                    else
                    {
                        if(rC==signWall)
                        {
                            edg2d[i].setRight(GeomElements::edge3d<2>::BoundaryType::WALL);
                        }
                        else if(rC==signOver)
                        {
                            edg2d[i].setRight(GeomElements::edge3d<2>::BoundaryType::FARFIELD);
                        }
                        else
                        {
                            throw std::runtime_error(" Not a valid boundary type\n");
                        }
                    }
                }
            }

            if(d_dim==2)
            {
                cur_ub2D->setLocalStructure(pt2d,cel2d,edg2d,nnode,ncell,nedge,0);
                //No need to specify local communication
            }
            delete[] xyzs;
            if(d_dim==2)
            {
                delete[] cells2D;
            }
            delete[] edge2Cell;
        }
        
    }


};