#include "OBBox.h"
#include "Operations.h"
OBBox3D::OBBox3D()
{
    mCenter = Vector3D(0, 0, 0);
    mAxis[0] = mAxis[1] = mAxis[2] = Vector3D(0, 0, 0);
    mHalfLength[0] = mHalfLength[1] = mHalfLength[2] = 0;
}
OBBox3D::OBBox3D(const double *xlo, const double *dx)
{
    for (int i = 0; i < 3; i++)
    {
        mHalfLength[i] = 0.5 * dx[i];
    }
    mAxis[0] = Vector3D(1, 0, 0);
    mAxis[1] = Vector3D(0, 1, 0);
    mAxis[2] = Vector3D(0, 0, 1);
    mCenter = Vector3D(xlo[0] + mHalfLength[0], xlo[1] + mHalfLength[1], xlo[2] + mHalfLength[2]);
}

OBBox3D::OBBox3D(const Vector3D &c, const Vector3D *pAxis, const double *half)
{
    mCenter = c;
    mAxis[0] = pAxis[0];
    mAxis[1] = pAxis[1];
    mAxis[2] = pAxis[2];

    mHalfLength[0] = half[0];
    mHalfLength[1] = half[1];
    mHalfLength[2] = half[2];
}
bool OBBox3D::isValid() const
{
    if (mAxis[0].dot_product(mAxis[1]) != 0 or mAxis[1].dot_product(mAxis[2]) != 0 or mAxis[2].dot_product(mAxis[0]) != 0)
    {
        return false;
    }
    for (int i = 0; i < 3; i++)
    {
        if (LEQUAL(mHalfLength[i], 0))
        {
            return false;
        }
        if (UNEQUAL(mAxis[i].dot_product(mAxis[i]), 1.0))
        {
            return false;
        }
    }
    return true;
}