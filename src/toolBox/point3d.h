#pragma once
#include"vector3d.h"
//A point3d is just a point living in a 3d world
//It should have been used to Make a Edge in 2d(face in 3d)
//This is just a wrapper for vector3d
namespace GeomElements
{
	template<int ndim>
	class point3d
	{
	public:
		vector3d<ndim, double> d_pos;
	public:
		point3d()
		{

		}
		point3d(const vector3d<ndim, double>& another)
			:d_pos(another)
		{

		}
		point3d(const point3d<ndim>& another)
			:d_pos(another.d_pos)
		{

		}
		vector3d<ndim,double> position() const
		{
			return d_pos;
		}
		point3d<ndim>& operator=(const point3d<ndim>& another)
		{
			d_pos = another.d_pos;
			return *this;
		}
		point3d<ndim>& operator=(const vector3d<ndim, double>& another)
		{
			d_pos = another;
			return *this;
		}
		double& operator[](int ind)
		{
			return d_pos[ind];
		}
		//The add and subtraction operations have no real meaning if a point is returned;
		//The subtraction of two points leads to a vector pointing from another to this
		vector3d<ndim, double> operator-(const point3d<ndim>& another)
		{
			return (d_pos - another.d_pos);
		}
		//The add of two points also leads to a vector;
		vector3d<ndim, double>operator+(const point3d<ndim>& another)
		{
			return (d_pos + another.d_pos);
		}




		//%TODO::There might need some complex operations like inclusion test with OBB.But for now we don't need to .

	};
};

template <int ndim>
std::ofstream &operator<<(std::ofstream &sm,const GeomElements::point3d<ndim> &p)
{
	//sm << "-----neo-cell-----\n";
	for (int i = 0; i < p.size(); i++)
	{
		sm << p.position()[i]<<'\t';//c.d_boundingPoints[c.pointInd(i)];
	}
	sm << '\n';
	return sm;
}