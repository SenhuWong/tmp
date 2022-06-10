#pragma once

#include"OBBox.h"
#include"triangle.h"
bool IsSeparatingAxis(double p0, double p1, double r);
bool TestAxisCrossEdge(const Vector3D *v, const Vector3D &e, const OBBox3D &box);
bool TriangleIntersectionOBB_SeparatingAxisMethod(const Triangle &tri, const OBBox3D &box);