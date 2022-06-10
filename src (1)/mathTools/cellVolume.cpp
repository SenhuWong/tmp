//This is only for 2D
#include<iostream>
#include"vector3d.h"
#include<cmath>
const double TOL = 0.0000000001;
double triangle_area(double* p1,double* p2,double*p3)
{
	double area = 0.5 * (
		(p1[0] - p2[0]) * (p1[1] + p2[1]) +
		(p2[0] - p3[0]) * (p2[1] + p3[1]) +
		(p3[0] - p1[0]) * (p3[1] + p1[1])
		);
	if(area<0)
	{
		area = -area;
	}
	return area;
}


double computeCellVolume2D(double xv[8][3], int nvert, bool counter_clock = false)
{
	double area;
	switch (nvert)
	{
	case 3:
	{
		area = triangle_area(xv[0], xv[1], xv[2]);
		return area;
		//Triangle
	}
	case 4:
	{
		area = 0;
		//find its correct orientation
		Vector3D AB = Vector3D(xv[1][0] - xv[0][0], xv[1][1] - xv[0][1], 0);
		Vector3D AC = Vector3D(xv[2][0] - xv[0][0], xv[2][1] - xv[0][1], 0);
		Vector3D AD = Vector3D(xv[3][0] - xv[0][0], xv[3][1] - xv[0][1], 0);
		//Do cross product
		Vector3D ABcAC = AB.cross_product(AC);
		//std::cout<<"ABcAC  = "<<ABcAC[0]<<" "<<ABcAC[1]<<" " <<ABcAC[2]<<'\n';
		Vector3D ABcAD = AB.cross_product(AD);
		if (ABcAC.dot_product(ABcAD) > 0)//C and D are on the same side of BA, which means BA is an edge.
		{
			//find ACcAD to see if B and D are on the same side of AD
			Vector3D ACcAD = AC.cross_product(AD);
			if (ACcAD.dot_product(ABcAD) > 0)//B and C are on the same side of AD,which means AD is an edge,and we know that AC is the diagnoal
			{
				//ACB+ADC
				area += triangle_area(xv[0], xv[2], xv[1]);
				area += triangle_area(xv[0], xv[2], xv[3]);
			}
			else//AD is the diagnoal
			{
				//ADC+ADB
				area += triangle_area(xv[0], xv[3], xv[1]);
				area += triangle_area(xv[0], xv[3], xv[2]);
			}
		}
		else//AB is the diagnoal
		{
			//ABC+ABD
			area += triangle_area(xv[0], xv[1], xv[2]);
			area += triangle_area(xv[0], xv[1], xv[3]);
			//Use the triangle formula

		}
		return area;
	}
	default:
		std::cout << "Other type of 2D cell is not implemented yet\n";
		return -1;
		break;
	}



}


double scalarProduct(double* a, double* b, double* c)
{
	double scalarProducts =
		a[0] * b[1] * c[2] - a[0] * b[2] * c[1]
		+ a[1] * b[2] * c[0] - a[1] * b[0] * c[2]
		+ a[2] * b[0] * c[1] - a[2] * b[1] * c[0];
	return scalarProducts;
}

void findCellVolume(double& vol, double coord[8][3],int numverts[4][6], int faceInfo[4][24],int itype, int* nface, int* nvert)
{
	vol = 0;
	for (int i = 0; i < *nface; i++)
	{
		if (numverts[itype][i] == 3)//This face has 3 vert
		{
			//c's index
			int inode = faceInfo[itype][4 * i + 2] -1 ;
			Vector3D c = Vector3D(coord[inode][0], coord[inode][1], coord[inode][2]);
			// std::cout << "inode is " << inode << " , " << c[0] << " " << c[1] << " " << c[2] << '\n';
			inode = faceInfo[itype][4 * i + 1] - 1;
			Vector3D b = Vector3D(coord[inode][0], coord[inode][1], coord[inode][2]);
			// std::cout << "inode is " << inode << " , " << b[0] << " " << b[1] << " " << b[2] << '\n';
			inode = faceInfo[itype][4 * i] -1;
			Vector3D a = Vector3D(coord[inode][0], coord[inode][1], coord[inode][2]);
			// std::cout << "inode is " << inode << " , " << a[0] << " " << a[1] << " " << a[2] << '\n';
			vol = vol - 0.5 * a.BoxProduct(b, c);
			// std::cout <<"face triangle's " << vol << '\n';
		}
		else if (numverts[itype][i] == 4)//This face has 4 vert
		{
			double aa[3], ba[3], ca[3],da[3];
			
			
			int inode = faceInfo[itype][4 * i + 3] -1;
			Vector3D d = Vector3D(coord[inode][0], coord[inode][1], coord[inode][2]);
			da[0] = coord[inode][0];
			da[1] = coord[inode][1];
			da[2] = coord[inode][2];
			// std::cout << "inode is " << inode << " , " << d[0] << " " << d[1] << " " << d[2] << '\n';
			inode = faceInfo[itype][4 * i + 2] -1;
			Vector3D c = Vector3D(coord[inode][0], coord[inode][1], coord[inode][2]);
			ca[0] = coord[inode][0];
			ca[1] = coord[inode][1];
			ca[2] = coord[inode][2];
			// std::cout << "inode is " << inode << " , " << c[0] << " " << c[1] << " " << c[2] << '\n';
			inode = faceInfo[itype][4 * i + 1] -1;
			ba[0] = coord[inode][0];
			ba[1] = coord[inode][1];
			ba[2] = coord[inode][2];
			Vector3D b = Vector3D(coord[inode][0], coord[inode][1], coord[inode][2]);
			// std::cout << "inode is " << inode << " , " << b[0] << " " << b[1] << " " << b[2] << '\n';
			inode = faceInfo[itype][4 * i] -1;
			Vector3D a = Vector3D(coord[inode][0], coord[inode][1], coord[inode][2]);
			aa[0] = coord[inode][0];
			aa[1] = coord[inode][1];
			aa[2] = coord[inode][2];
			// std::cout << "inode is " << inode << " , " << a[0] << " " << a[1] << " " << a[2] << '\n';
			if (false)
			{
				if (fabsf64(c.BoxProduct(a, b) - scalarProduct(aa, ba, ca)) < TOL)
				{
					std::cout << "Something Wrong!!!!!!!!!!!!!!!\n";
					std::cout << c.BoxProduct(a, b) << '\t' << scalarProduct(aa, ba, ca) << '\t' << c.BoxProduct(a, b) - scalarProduct(aa, ba, ca) << '\n';
				}
				if (fabsf64(d.BoxProduct(a, c) - scalarProduct(aa, ca, da)) < TOL)
				{
					std::cout << "Something Wrong!!!!!!!!!!!!!!!\n";
					std::cout << d.BoxProduct(a, c) << '\t' << scalarProduct(aa, ca, da) << '\t' << d.BoxProduct(a, c) - scalarProduct(aa, ca, da) << '\n';
				}
				if (fabsf64(d.BoxProduct(a, b) - scalarProduct(aa, ba, da)) < TOL)
				{
					std::cout << "Something Wrong!!!!!!!!!!!!!!!\n";
					std::cout << d.BoxProduct(a, b) << '\t' << scalarProduct(aa, ba, da) << '\t' << d.BoxProduct(a, b) - scalarProduct(aa, ba, da) << '\n';
				}
				if (fabsf64(d.BoxProduct(b, c) - scalarProduct(ba, ca, da)) < TOL)
				{
					std::cout << "Something Wrong!!!!!!!!!!!!!!!\n";
					std::cout << d.BoxProduct(b, c) << '\t' << scalarProduct(ba, ca, da) << '\t' << d.BoxProduct(b, c) - scalarProduct(ba, ca, da) << '\n';
				}
			}
			
			// vol = vol - 0.25 * c.BoxProduct(a, b);
			// vol = vol - 0.25 * d.BoxProduct(a, c);
			// vol = vol - 0.25 * d.BoxProduct(a, b);
			// vol = vol - 0.25 * d.BoxProduct(b, c);
			
			vol = vol - 0.25 * scalarProduct(aa, ba, ca);
			vol = vol - 0.25 * scalarProduct(aa, ca, da);
			vol = vol - 0.25 * scalarProduct(aa, ba, da);
			vol = vol - 0.25 * scalarProduct(ba, ca, da);
			// std::cout << "face quadrilateral's " << vol << '\n';
		}
	}
	vol = vol / 3.0;
	if(vol<=0)
	{
		std::cout<<*nface<<":"<<vol<<'\n';

		std::cin.get();
	}
	return;
}

double computeCellVolume3D(double xv[8][3], int nvert)
{
	double vol;
	double vol1;
	int itype;
	int nfaces;
	int numverts[4][6] = {
	3,3,3,3,0,0,
	4,3,3,3,3,0,
	3,4,4,4,3,0,
	4,4,4,4,4,4 };
	int faceInfo[4][24] = { 1,2,3,0,1,4,2,0,2,4,3,0,1,3,4,0,0,0,0,0,0,0,0,0,
						  1,2,3,4,1,5,2,0,2,5,3,0,4,3,5,0,1,4,5,0,0,0,0,0,
						  1,2,3,0,1,4,5,2,2,5,6,3,1,3,6,4,4,6,5,0,0,0,0,0,
						  1,2,3,4,1,5,6,2,2,6,7,3,3,7,8,4,1,4,8,5,5,8,7,6 };
	switch (nvert)
	{
	case 4:
		itype = 0;
		nfaces = 4;
		break;
	case 5:
		itype = 1;
		nfaces = 5;
		break;
	case 6:
		itype = 2;
		nfaces = 5;
		break;
	case 8:
		itype = 3;
		nfaces = 6;
		break;
	}

	findCellVolume(vol, xv, numverts, faceInfo,itype, &nfaces, &nvert);

	return vol;
}

double computeVolume(int dim,double xv[8][3], int nvert,bool counter_clock = false)
{
	if (dim==2)
	{
		return computeCellVolume2D(xv, nvert, counter_clock);
	}
	else if(dim==3)
	{
		return computeCellVolume3D(xv, nvert);
	}

}