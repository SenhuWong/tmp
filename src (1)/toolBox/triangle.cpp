#pragma once
#include"triangle.h"
#include <iostream>
template<int ndim,typename T>
GeomElements::vector3d<ndim, T> Triangle<ndim, T>::find_normal()
{
	if (ndim == 2)
		//Report an error and return.
		std::cout << "Something wrong\n";

	//Assume verts are arranged in certrain direction
	GeomElements::vector3d<ndim, T> result = (verts[1] - verts[0]).cross_product(verts[2] - verts[0]);
	result.normalize();
	
}