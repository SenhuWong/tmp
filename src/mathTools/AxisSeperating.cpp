#include "AxisSeperating.h"
#include "Operations.h"
bool IsSeparatingAxis(double p0, double p1, double r)
{
    double min = p0;
    double max = p0;
    if (min > p1)
    {
        min = p1;
    }
    if (max < p1)
    {
        max = p1;
    }
    if (LESS(max, -r) or GREATER(min, r))
    {
        return true;
    }
    return false;
}
bool TestAxisCrossEdge(const Vector3D *v, const Vector3D &e, const OBBox3D &box)
{
    double fes[3], r;
    for (int i = 0; i < 3; i++)
    {
        fes[i] = fabs(e[i]);
    }
    double p0, p2;
    //e and Axis X
    //n = (1,0,0)[cross]e = (0,-e[2],e[1])
    p0 = -e[2] * v[0][1] + e[1] * v[0][2];
    p2 = -e[2] * v[2][1] + e[1] * v[2][2];
    r = fes[2] * box.mHalfLength[1] + fes[1] * box.mHalfLength[2];
    if (IsSeparatingAxis(p0, p2, r))
        return false;
    //e and Axis Y
    //n = (0,1,0)[cross]e = (e[2],0,-e[0])
    p0 = -e[0] * v[0][2] + e[2] * v[0][0];
    p2 = -e[0] * v[2][2] + e[2] * v[2][0];
    r = fes[0] * box.mHalfLength[2] + fes[2] * box.mHalfLength[0];
    if (IsSeparatingAxis(p0, p2, r))
    {
        return false;
    }
    //e and Axis Z
    //n = (0,0,1)[cross]e = (-e[1],e[0],0)
    p0 = -e[1] * v[0][0] + e[0] * v[0][1];
    p2 = -e[1] * v[2][0] + e[0] * v[2][1];
    r = fes[1] * box.mHalfLength[0] + fes[0] * box.mHalfLength[1];
    if (IsSeparatingAxis(p0, p2, r))
    {
        return false;
    }
    return true;
}
bool TriangleIntersectionOBB_SeparatingAxisMethod(const Triangle &tri, const OBBox3D &box)
{
    Vector3D v[3];
    Vector3D p;
    for (int i = 0; i < 3; i++)
    {
        p = tri.mPoints[i] - box.mCenter;
        v[i][0] = box.mAxis[0].dot_product(p);
        v[i][1] = box.mAxis[1].dot_product(p);
        v[i][2] = box.mAxis[2].dot_product(p);
    }
    //Test the 9 first
    const Vector3D e0 = v[1] - v[0];
    Vector3D u[3];
    u[0] = v[0];
    u[1] = v[1];
    u[2] = v[2];
    if (!TestAxisCrossEdge(u, e0, box))
    {
        return false;
    }
    const Vector3D e1 = v[2] - v[1];
    u[0] = v[1];
    u[1] = v[2];
    u[2] = v[0];
    if (!TestAxisCrossEdge(u, e1, box))
    {
        return false;
    }
    const Vector3D e2 = v[0] - v[2];
    u[0] = v[2];
    u[1] = v[0];
    u[2] = v[1];
    if (!TestAxisCrossEdge(u, e2, box))
    {
        return false;
    }
    //Find minmax overlap in (x,y,z) directions
    double min, max;
    for (int i = 0; i < 3; i++) //For each axis v[][i]
    {
        min = max = v[0][i];
        for (int j = 1; j < 3; j++) //For each point along this axis v[j][]
        {
            if (v[j][i] < min)
            {
                min = v[j][i];
            }
            if (v[j][i] > max)
            {
                max = v[j][i];
            }
        }
        if (GREATER(min, box.mHalfLength[i]) or LESS(max, box.mHalfLength[i]))
        {
            return false;
        }
    }
    //Test if the box intersect the plane of the triangle.
    Vector3D normal = (v[1] - v[0]).cross_product(v[2] - v[0]);
    double d = -normal.dot_product(v[0]);
    Vector3D minPoint, maxPOint;
    for (int i = 0; i < 3; i++)
    {
        if (GREATER(normal[i], 0)) //If normal vector is positive
        {
            minPoint[i] = -box.mHalfLength[i]; //min is at negative
            maxPOint[i] = box.mHalfLength[i];
        }
        else
        {
            minPoint[i] = box.mHalfLength[i];
            maxPOint[i] = -box.mHalfLength[i];
        }
    }
    double t;
    t = normal.dot_product(minPoint) + d;
    if (GREATER(t, 0))
        return false;
    t = normal.dot_product(maxPOint) + d;
    if (LESS(t, 0))
        return false;
    return true;
}