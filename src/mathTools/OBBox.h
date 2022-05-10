#pragma once

#include "vector3d.h"
class OBBox3D
{
public:
    OBBox3D();
    OBBox3D(const double *xlo, const double *dx); //A single cube
    OBBox3D(const Vector3D &c, const Vector3D *pAxis, const double *half);
    bool isValid() const;

    Vector3D mCenter;
    Vector3D mAxis[3];
    double mHalfLength[3];
};

class OBBox2D
{
public:
    OBBox2D();
    OBBox2D(const double *xlo, const double *dx);
    OBBox2D(const Vector2D &c, const Vector2D *pAxis, const double *half);
    bool isValid() const;
    Vector2D mCenter;
    Vector2D mAxis[3];
    double mHlafLength[2];
};


//I will implement this later.
// class cartOBBox
// {
// private:
//     int d_dim;
//     double d_xlo[3];
//     double d_xhi[3];
//     cartOBBox();
// public:
//     cartOBBox(int dim);
//     cartOBBox(int dim,const double *xlo, const double *dx);
//     cartOBBox(const cartOBBox& another);
//     void findCartOBBox(double* rxyz,int nxyz);
// };

// cartOBBox::cartOBBox(int dim)
// {
//     d_dim = dim;
// }

// cartOBBox::cartOBBox(int dim, const double* xlo, const double* dx)
// {
//     d_dim = dim;
//     for(int i=0;i<dim;i++)
//     {
//         d_xlo[i] = 
//     }
// }