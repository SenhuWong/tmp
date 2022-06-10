#pragma once
#include <vector>
#include<iostream>
#include "../reader/Reader.h"
class SamraiFeeder:virtual public Reader
{
    int insideTag;
    int outsideTag;
    int crossedTag;
protected:
    double **meshes_xyz = NULL;
    std::vector<int> *meshes_overInd = NULL;
    std::vector<int> *meshes_overPtr = NULL;

    std::vector<int> *meshes_wallInd = NULL;
    std::vector<int> *meshes_wallPtr = NULL;

public:
    void setTag(int inside,int outside,int crossed)
    {
        insideTag = inside;
        outsideTag = outside;
        crossedTag = crossed;
    }
    virtual void readForSamrai() = 0;
    void initialize()
    {
        meshes_xyz = new double *[d_nmesh];
        meshes_overPtr = new std::vector<int>[d_nmesh];
        meshes_overInd = new std::vector<int>[d_nmesh];
        meshes_wallInd = new std::vector<int>[d_nmesh];
        meshes_wallPtr = new std::vector<int>[d_nmesh];
    }
    void getOverBoundary(int &nmeshes, int **nOverEach, int **nWallEach, double ***overXyzEach, int ***overIndEach, int ***overPtrEach, int ***wallIndEach, int ***wallPtrEach)
    {
        nmeshes = d_nmesh;
        *nOverEach = new int[nmeshes];
        *nWallEach = new int[nmeshes];
        *overXyzEach = meshes_xyz;
        *overIndEach = new int *[nmeshes];
        *overPtrEach = new int *[nmeshes];
        *wallIndEach = new int *[nmeshes];
        *wallPtrEach = new int *[nmeshes];
        for (int i = 0; i < nmeshes; i++)
        {
            (*nOverEach)[i] = meshes_overPtr[i].size() - 1;
            (*nWallEach)[i] = meshes_wallPtr[i].size() - 1;
            //std::cout << "overEach for " << i << " is " << *nOverEach[i] << '\n';
            (*overIndEach)[i] = meshes_overInd[i].data();
            (*overPtrEach)[i] = meshes_overPtr[i].data();
            
            (*wallIndEach)[i] = meshes_wallInd[i].data();
            (*wallPtrEach)[i] = meshes_wallPtr[i].data();
        }
        std::cout<<"OversetBoundary get\n";
    }

    SamraiFeeder();
    ~SamraiFeeder();
};