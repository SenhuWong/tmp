#pragma once
#include "vector3d.h"

class Triangle
{
public:
    //sizeof(float)*4*3+sizeof(unsigned short int) == 50

    //Copy and Constructors
    Triangle();
    Triangle(const Triangle &another);
    Triangle &operator=(const Triangle &another);
    //Given 3 points, find its normal unit vector
    void findNormal();
    static Triangle *make_Triangles(int nvert, float *xyzs);
    int cubeCrossing(double *xlo, double *dx);
    // bool segmentHitTest(const Vector3D &point1, const Vector3D &point2) const
    // {
    //     Vector3D ret;
    //     Vector3D direction = point2 - point1;
    //     if (!linearIntersectTriangle(point1, direction, ret))
    //     {
    //         return false;
    //     }
    //     if (LESS(ret[2], 0) or GREATER(ret[2], 1))
    //     {
    //         return false;
    //     }
    //     return true;
    // }
    Vector3D normal;
    Vector3D mPoints[3];
    unsigned short int source;

//private:
    // bool linearIntersectTriangle(const Vector3D &base, const Vector3D &direction, Vector3D &result) const
    // {
    //     double u, v, tmp;
    //     Vector3D e1, e2, p, s, q;
    //     e1 = mPoints[1] - mPoints[0];
    //     e2 = mPoints[2] - mPoints[0];
    //     p = direction.cross_product(e2);
    //     tmp = p.dot_product(e1);
    //     if (EQUAL(tmp, 0)) //line is perpendicular to normal of triangle
    //     {
    //         p = e1.cross_product(e2);
    //         tmp = p.dot_product(base - mPoints[0]);
    //         if (EQUAL(tmp, 0))
    //         {
    //             return true;
    //         }
    //         return false;
    //     }
    //     s = base - mPoints[0];
    //     u = p.dot_product(s) / tmp;
    //     if (LESS(u, 0) or GREATER(u, 1))
    //     {
    //         return false;
    //     }
    //     q = s.cross_product(e1);
    //     v = q.dot_product(direction) / tmp;
    //     if (LESS(v, 0) or GREATER(v, 1))
    //     {
    //         return false;
    //     }

    //     result[0] = u;
    //     result[1] = v;
    //     result[2] = q.dot_product(e2) / tmp;
    //     return true;
    // }
};