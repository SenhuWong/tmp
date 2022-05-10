#include "SamraiFeeder.h"

SamraiFeeder::SamraiFeeder()
{
}
SamraiFeeder::~SamraiFeeder()
{
    if (d_nmesh != 0)
    {
        for (int i = 0; i < d_nmesh; i++)
        {
            delete[] meshes_xyz[i];
        }
        delete[] meshes_overInd;
        delete[] meshes_overPtr;
        delete[] meshes_wallInd;
        delete[] meshes_wallPtr;
        delete[] meshes_xyz;
    }
}
