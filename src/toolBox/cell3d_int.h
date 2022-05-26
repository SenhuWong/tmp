#pragma once
#include "point3d.h"
double computeVolume(int dim, double xv[8][3], int nvert, bool counter_clock = false);

namespace GeomElements
{
	template <int ndim>
	class edge3d;
	template <int ndim>
	class cell3d
	{
		static double Tolerance;
		static point3d<ndim> *d_boundingPoints;
		int d_vertsInd[4 * (ndim - 1)] = {-1};
		static edge3d<ndim> *d_boundingEdges;
		int d_edgesInd[2*ndim] = {-1};
		vector3d<ndim, double> d_center;
		double d_volume = 0;
		unsigned int d_size = 0;
		unsigned int d_edge_size = 0;

	public:
		static void bindPoints(point3d<ndim> *toBind)
		{
			cell3d<ndim>::d_boundingPoints = toBind;
			// d_nboundingPoints = ntoBound;
		}
		cell3d()
		{
		}
		void push_back(int vertInd)
		{
			d_vertsInd[d_size++] = vertInd;
		}

		void push_edge(int edgeInd)
		{
			d_edgesInd[d_edge_size++] = edgeInd;
		}

		void clear()
		{
			d_size = 0;
			
		}
		int pointInd(int ind) const
		{
			return d_vertsInd[ind];
		}
		int size() const
		{
			return d_size;
		}
		int edgeInd(int ind) const
		{
			return d_edgesInd[ind];
		}
		int edge_size() const
		{
			return d_edge_size;
		}
		
		double volume() const
		{
			return d_volume;
		}
		vector3d<ndim, double> center()
		{
			return d_center;
		}
		
		
		void computeMetaData()
		{
			find_volume();
			find_center();
		}
	private:
		vector3d<ndim, double> &point(int ind) const
		{
			return cell3d<ndim>::d_boundingPoints[d_vertsInd[ind]].d_pos;
		}
		void find_volume()
		{
			double xv[(4 * (ndim - 1))][3];
			for (int i = 0; i < d_size; i++)
			{
				for (int j = 0; j < ndim; j++)
				{
					xv[i][j] = cell3d<ndim>::d_boundingPoints[pointInd(i)][j];
					// std::cout<<xv[i][j]<<'\t';
				}
				// std::cout<<'\n';
			}
			double result = computeVolume(ndim, xv, d_size);
			// std::cout<<result<<'\n';
			d_volume = computeVolume(ndim, xv, d_size);
			// std::cout<<d_volume<<'\n';
		}
		void find_center()
		{
			vector3d<ndim, double> result(0, 0, 0);
			for (int i = 0; i < d_size; i++)
			{
				result = result + point(i);
			}
			d_center = result / d_size;
		}
	};
	// template<int ndim> int d_nboundingPoints = 0;
	template <int ndim>
	point3d<ndim> *cell3d<ndim>::d_boundingPoints = NULL;

};
