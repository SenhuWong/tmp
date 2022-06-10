#include <cmath>
#include "Operations.h"
#include "triangle.h"
#include "AxisSeperating.h"
Triangle::Triangle()
{
}
Triangle::Triangle(const Triangle &another)
{
    for (int i = 0; i < 3; i++)
    {
        this->mPoints[i] = another.mPoints[i];
    }
    this->normal = another.normal;
    this->source = another.source;
}
Triangle &Triangle::operator=(const Triangle &another)
{

    for (int i = 0; i < 3; i++)
    {
        this->mPoints[i] = another.mPoints[i];
    }
    this->normal = another.normal;
    this->source = another.source;
    return *this;
}

void Triangle::findNormal()
{
    Vector3D e0 = mPoints[1] - mPoints[0];
    Vector3D e1 = mPoints[2] - mPoints[1];
    Vector3D normalUnit = e0.cross_product(e1).normalized();
    normal = normalUnit;
}

Triangle *Triangle::make_Triangles(int nvert, float *xyzs)
{
    Triangle *result;
    if (nvert == 3)
    {
        result = new Triangle[1];
        for (int i = 0; i < 3; i++)
        {
            result->mPoints[i] = Vector3D(xyzs[3 * i + 0], xyzs[3 * i + 1], xyzs[3 * i + 2]);
        }
        result[0].findNormal();
    }
    else if (nvert == 4)
    {
        result = new Triangle[2];
        //split along the shorter one
        //02 and 13
        float len02 = 0;
        float len13 = 0;
        for (int i = 0; i < 3; i++)
        {
            len02 += std::pow(xyzs[3 * 0 + i] - xyzs[3 * 2 + i], 2);
            len13 += std::pow(xyzs[3 * 1 + i] - xyzs[3 * 3 + i], 2);
        }
        if (len02 < len13) //012 and 023
        {
            for (int i = 0; i < 3; i++)
            {
                result[0].mPoints[0][i] = xyzs[3 * 0 + i];
                result[0].mPoints[1][i] = xyzs[3 * 1 + i];
                result[0].mPoints[2][i] = xyzs[3 * 2 + i];
                result[1].mPoints[0][i] = xyzs[3 * 0 + i];
                result[1].mPoints[1][i] = xyzs[3 * 2 + i];
                result[1].mPoints[2][i] = xyzs[3 * 3 + i];
            }
        }
        else //013 and 123
        {
            for (int i = 0; i < 3; i++)
            {
                result[0].mPoints[0][i] = xyzs[3 * 0 + i];
                result[0].mPoints[1][i] = xyzs[3 * 1 + i];
                result[0].mPoints[2][i] = xyzs[3 * 3 + i];
                result[1].mPoints[0][i] = xyzs[3 * 1 + i];
                result[1].mPoints[1][i] = xyzs[3 * 2 + i];
                result[1].mPoints[2][i] = xyzs[3 * 3 + i];
            }
        }
        result[0].findNormal();
        result[1].findNormal();
    }
    return result;
}

int Triangle::cubeCrossing(double *xlo, double *dx)
{
    
    OBBox3D box = OBBox3D(xlo,dx);
    //Maybe we should not depend on Patch
    bool crossed = TriangleIntersectionOBB_SeparatingAxisMethod(*this,box);
}