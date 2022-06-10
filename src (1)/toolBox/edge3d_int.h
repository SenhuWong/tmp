#pragma once
#include "cell3d_int.h"
#include <vector>
#include <memory>
#define UNSTRUCT_DEBUG
namespace GeomElements
{
	template <int ndim>
	class edge3d
	{

		static double Tolerance;
		static point3d<ndim> *d_boundingPoints;
		// static point3d<ndim>* d_nboundingPoints;

		int d_vertsInd[2 * (ndim - 1)]; // 2 in 2D, 4 in 3D
		unsigned int d_size = 0;
		int leftInd = BoundaryType::UNSET;
		int rightInd = BoundaryType::UNSET;
		vector3d<ndim, double> d_normal;
		vector3d<ndim, double> d_center;
		double d_area = 0;
#ifdef UNSTRUCT_DEBUG
		int needsToCheck = 0;
		void setNeedsToCheck(int needto)
		{
			needsToCheck = needto;
		}

#endif // DEBUG

	public:
		enum BoundaryType
		{
			UNSET = -1,
			WALL = -2,
			FARFIELD = -3,
			OVERSET = -4,
			SYMMETRY = -5
		};
		static void bindPoints(point3d<ndim> *toBind);
		edge3d()
		{
		}
		BoundaryType type()
		{
			return rightInd;
		}
		void push_back(int vertInd)
		{
			d_vertsInd[d_size++] = vertInd;
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
		int lCInd() const
		{
			return leftInd;
		}
		int rCInd() const
		{
			return rightInd;
		}
		void setLeft(int ind)
		{
			leftInd = ind;
		}
		void setRight(int ind)
		{
			rightInd = ind;
		}
		vector3d<ndim, double> center() const
		{
			return d_center;
		}
		vector3d<ndim, double> normal_vector() const
		{
			return d_normal;
		}
		double area() const
		{
			return d_area;
		}
		void computeMetaData()
		{
			find_normal();
			find_area();
			find_center();
		}
	private:
		vector3d<ndim, double> &point(int ind) const
		{
			return edge3d<ndim>::d_boundingPoints[d_vertsInd[ind]].d_pos;
		}
		
		// void make_clockwise();
		void find_center()
		{
			vector3d<ndim, double> result(0, 0, 0);
			for (int i = 0; i < d_size; i++)
			{
				result = result + point(i);
			}
			d_center = result / d_size;
			
		}
		void find_normal();
		void find_area();
		bool is_plane();
		bool is_clockwise();
	};
	template <int ndim>
	point3d<ndim> *edge3d<ndim>::d_boundingPoints = NULL;
	template <int ndim>
	// This function must be called only after static binding is done
	bool edge3d<ndim>::is_plane()
	{
		if (ndim == 2)
		{
			return true;
		}
		else if (ndim == 3)
		{
			if (d_size == 3)
				return true;
			else
			{
				GeomElements::vector3d<ndim, double> AB_cross_AC = (this->point(1) - this->point(0)).cross_product(this->point(2) - this->point(0));
				GeomElements::vector3d<ndim, double> AB_cross_AD = (this->point(1) - this->point(0)).cross_product(this->point(3) - this->point(0));
				GeomElements::vector3d<ndim, double> result = AB_cross_AC.cross_product(AB_cross_AD);
				for (int i = 0; i < ndim; i++)
				{
					if (result[i] < Tolerance)
					{
						return false;
					}
				}
				return true;
			}
		}
	}
	template <int ndim>
	double edge3d<ndim>::Tolerance = 1.0e-10;
	template <int ndim>
	void edge3d<ndim>::bindPoints(point3d<ndim> *toBind)
	{
		edge3d<ndim>::d_boundingPoints = toBind;
	}
	// This function must be called only after static binding is done
	template <int ndim>
	bool edge3d<ndim>::is_clockwise()
	{
		if (ndim == 2)
		{
			return true;
		}
		else if (ndim == 3)
		{
			if (d_size == 3)
				return true;
			else if (d_size == 4)
			{
				GeomElements::vector3d<ndim, double> AB_cross_AC = (this->point(1) - this->point(0)).cross_product(this->point(2) - this->point(0));
				GeomElements::vector3d<ndim, double> AB_cross_AD = (this->point(1) - this->point(0)).cross_product(this->point(3) - this->point(0));
				double CD_sameSide = AB_cross_AC.dot_product(AB_cross_AD);
				if (CD_sameSide > 0) // C and D are on the same side of AB
				{
					auto BC = this->point(2) - this->point(1);
					GeomElements::vector3d<ndim, double> BC_cross_BA = BC.cross_product(this->point(0) - this->point(1));
					GeomElements::vector3d<ndim, double> BC_cross_BD = BC.cross_product(this->point(3) - this->point(1));
					double AD_sameSide = BC_cross_BA.dot_product(BC_cross_BD);
					if (AD_sameSide > 0)
					{
						return true;
					}
					else
					{
						return false;
					}
				}
				else // AB is a diagonal
				{
					return false;
				}
				return true;
			}
		}
	}
	// This function must be called only after static binding is done
	template <int ndim>
	void edge3d<ndim>::find_normal()
	{
		find_area();
		if (ndim == 2)
		{
			// In 2d case, the normal vector would be the right hand side along the direction of the edge
			vector3d<ndim, double> direction = this->point(1) - this->point(0);
			d_normal[0] = direction[1] / d_area;
			d_normal[1] = -direction[0] / d_area;
		}
		else if (ndim == 3)
		{
			// In 3d case, the normal vector would be AB.cross(AC) after normalization.
			d_normal = (this->point(1) - this->point(0)).cross_product(this->point(2) - this->point(0));
			d_normal.normalize();
		}
	}
	// This function must be called only after static binding is done
	template <int ndim>
	void edge3d<ndim>::find_area()
	{
		if (ndim == 2)
		{
			vector3d<ndim, double> direction = this->point(1) - this->point(0);
			d_area = direction.normalize();
		}
		else if (ndim == 3)
		{
			d_area = 0;
			if (d_size == 3) // Triangle
			{
				d_area += ((this->point(1) - this->point(0)).cross_product(this->point(2) - this->point(0))).normalize() / 2;
			}
			else if (d_size == 4) // Quadrilateral but only in clockwise or counter clockwise
			{
				d_area += ((this->point(1) - this->point(0)).cross_product(this->point(2) - this->point(0))).normalize() / 2;
				d_area += ((this->point(2) - this->point(0)).cross_product(this->point(3) - this->point(0))).normalize() / 2;
			}
		}
	}
}
