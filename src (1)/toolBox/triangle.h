#pragma once
#include"vector3d.h"
template<int ndim,typename T>
class Triangle
{
	GeomElements::vector3d<ndim, T> verts[3];
	Triangle<ndim, T>();
	Triangle<ndim, T>(const GeomElements::vector3d<ndim, T>& p1, const GeomElements::vector3d<ndim, T>& p2, const GeomElements::vector3d<ndim, T>& p3);
	Triangle<ndim, T>(const Triangle<ndim, T>& another);

	GeomElements::vector3d<ndim, T> find_normal();
};