#include<iostream>
#include<math.h>
#include"vector3d.h"
void solvec(double** a, double* b, int* iflag, int n)
{
	int i, j, k, l, flag;//, templ;
	double fact;
	double temp;
	double sum;
	double eps = 1e-8;
	for (i = 0; i < n; i++)
	{
		if (fabs(a[i][i]) < eps)//��Ԫ̫С
		{
			flag = 1;
			for (k = i + 1; k < n and flag; k++)
			{
				//std::cout << "a[k][i]  and a[i][i] is " << a[k][i] <<'\t'<<a[i][i] << '\n';
				if (!(fabs(a[k][i]) < eps))//�ҵ���С����Ԫ
				{
					flag = 0;
					for (l = 0; l < n; l++)
					{
						temp = a[k][l];
						a[k][l] = a[i][l];
						a[i][l] = temp;
					}
					temp = b[k];
					b[k] = b[i];
					b[i] = temp;
				}
			}
			if (flag) { std::cout << "Didnt find\n"; *iflag = 0; return; }//�Ҳ����Ͳ��ܼ���
		}
		for (k = i + 1; k < n; k++)
		{
			if (i != k)//I don't know why this is here
			{
				fact = -a[k][i] / a[i][i];
				for (j = 0; j < n; j++)
				{
					a[k][j] += fact * a[i][j];
				}
				b[k] += fact * b[i];

			}
		}
	}
	for (i = n - 1; i >= 0; i--)
	{
		sum = 0;
		for (j = i + 1; j < n; j++)
		{
			sum += a[i][j] * b[j];
		}
		b[i] = (b[i] - sum) / a[i][i];
	}
	*iflag = 1;
	return;
}

void newtonSolve3D(double f[7][3], double* u1, double* v1, double* w1)
{
	int i, j, k;
	int iter, itmax, isolflag;
	double u, v, w;
	double uv, wu, vw, uvw, norm, convergenceLimit;
	double* rhs;
	double** lhs;
	double alph;
	//
	lhs = new double* [3];
	for (i = 0; i < 3; i++)
		lhs[i] = new double[3];
	rhs = new double[3];
	//
	itmax = 500;
	convergenceLimit = 1e-14;
	alph = 1.0;
	isolflag = 1.0;
	//
	u = v = w = 0.5;
	//
	for (iter = 0; iter < itmax; iter++)
	{
		uv = u * v;
		vw = v * w;
		wu = w * u;
		uvw = u * v * w;

		for (j = 0; j < 3; j++)
			rhs[j] = f[0][j] + f[1][j] * u + f[2][j] * v + f[3][j] * w +
			f[4][j] * uv + f[5][j] * vw + f[6][j] * wu +
			f[7][j] * uvw;

		norm = rhs[0] * rhs[0] + rhs[1] * rhs[1] + rhs[2] * rhs[2];
		if (sqrt(norm) <= convergenceLimit) break;

		for (j = 0; j < 3; j++)
		{
			lhs[j][0] = f[1][j] + f[4][j] * v + f[6][j] * w + f[7][j] * vw;
			lhs[j][1] = f[2][j] + f[5][j] * w + f[4][j] * u + f[7][j] * wu;
			lhs[j][2] = f[3][j] + f[6][j] * u + f[5][j] * v + f[7][j] * uv;
		}

		solvec(lhs, rhs, &isolflag, 3);
		if (isolflag == 0) break;

		u -= (rhs[0] * alph);
		v -= (rhs[1] * alph);
		w -= (rhs[2] * alph);
	}
	if (iter == itmax) { u = 2.0; v = w = 0.; }
	if (isolflag == 0) {
		u = 2.0;
		v = w = 0.;
	}
	*u1 = u;
	*v1 = v;
	*w1 = w;
	for (i = 0; i < 3; i++) delete[]lhs[i];
	delete[]lhs;
	delete[]rhs;
	return;
}
void newtonSolve2D(double** f, double* u1, double* v1)//for quadraliteral only
{
	double** lhs = new double* [2];
	for (int i = 0; i < 2; i++)
	{
		lhs[i] = new double[2];
	}
	double* rhs = new double[2];
	int itmax = 500;
	double convergenceLimit = 1e-14;
	double alph = 1.0;
	int isoflag = 1;

	double u = 0.5;
	double v = 0.5;
	double uv;
	double norm;
	int iter;
	for (iter = 0; iter < itmax; iter++)
	{
		uv = u * v;
		for (int i = 0; i < 2; i++)
		{
			rhs[i] = f[0][i] + f[1][i] * u + f[2][i] * v + f[3][i] * uv;//getting f
		}
		norm = rhs[0] * rhs[0] + rhs[1] * rhs[1];
		if (sqrt(norm) <= convergenceLimit)
		{
			break;
		}
		for (int i = 0; i < 2; i++)
		{
			lhs[i][0] = f[1][i] + f[3][i] * v;//Each derivative occupies a colume, and each colume need to be specified.
			lhs[i][1] = f[2][i] + f[3][i] * u;
		}
		solvec(lhs, rhs, &isoflag, 2);
		if (isoflag == 0) break;//Not a good solution

		u -= rhs[0] * alph;
		v -= rhs[1] * alph;
	}
	if (iter == itmax)
	{
		u = 2.0;
		v = 0.0;
	}
	if (isoflag == 0)
	{
		u = 2.0;
		v = 0.0;
	}

	*u1 = u;
	*v1 = v;
	for (int i = 0; i < 2; i++)
	{
		delete[] lhs[i];
	}
	delete[] lhs;
	delete[] rhs;
	return;
}

void findClockWise(int* mapping,double xv[8][3])
{
	//find its correct orientation
	Vector3D AB = Vector3D(xv[1][0] - xv[0][0], xv[1][1] - xv[0][1], 0);
	Vector3D AC = Vector3D(xv[2][0] - xv[0][0], xv[2][1] - xv[0][1], 0);
	Vector3D AD = Vector3D(xv[3][0] - xv[0][0], xv[3][1] - xv[0][1], 0);

	//Do cross product

	Vector3D ABcAC = AB.cross_product(AC);
	Vector3D ABcAD = AB.cross_product(AD);
	if (ABcAC.dot_product(ABcAD) > 0)//C and D are on the same side of BA, which means BA is an edge.
	{
		//find ACcAD to see if B and D are on the same side of AD
		Vector3D ACcAD = AC.cross_product(AD);
		if (ACcAD.dot_product(ABcAD) > 0)//B and C are on the same side of AD,which means AD is an edge,and we know that AC is the diagnoal
		{
			return;
			//ACB+ADC
		}
		else//AD is the diagnoal
		{
			mapping[2] = 3;
			mapping[3] = 2;
			//0132
			//ADC+ADB
		}
	}
	else//AB is the diagnoal
	{
		//ABC+ABD
		mapping[1] = 2;
		mapping[2] = 1;
		//0213
		
		//Use the triangle formula

	}

}
void computeNodalWeight2D(double xv[8][3],double* xp,double* frac,int nvert)
{
	int isoflag = 0;
	switch (nvert)
	{
	case 3://Triangle
	{
		double** lhs;
		lhs = new double* [2];
		for (int i = 0; i < 2; i++)
		{
			lhs[i] = new double[2];
		}
		double* rhs = new double[2];
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				lhs[i][j] = xv[j][i] - xv[2][i];//
			}
			rhs[i] = xp[i] - xv[2][i];
		}
		//Solving the equation
		solvec(lhs, rhs, &isoflag, 2);
		if (isoflag)
		{
			for (int k = 0; k < 2; k++)
			{
				frac[k] = rhs[k];
			}
			frac[2] = 1. - frac[0] - frac[1];
		}
		else//Matrix is formed with vectors between triangle points, it shouldn't degenerate, which means two vector is parallel
		{
			throw std::runtime_error("Error matrix degenerate\n");
			frac[0] = 1.0;
			frac[1] = frac[2] = 0;
		}

		for (int i = 0; i < 2; i++)
		{
			delete[] lhs[i];
		}
		delete[] lhs;
		delete[] rhs;
	}
	break;
	case 4://Quadraliteral
	{
		double** f = new double* [4];
		for (int i = 0; i < 4; i++)
		{
			f[i] = new double[2];
		}
		double u, v;
		//First we need to find the counter clock wise direction mapping 
		int mapping[4] = { 0,1,2,3 };
		findClockWise(mapping, xv);
		for (int i = 0; i < 2; i++)
		{
			//v1-p
			f[0][i] = xv[mapping[0]][i] - xp[i];
			//v2-v1
			f[1][i] = xv[mapping[1]][i] - xv[mapping[0]][i];
			//v4-v1
			f[2][i] = xv[mapping[3]][i] - xv[mapping[0]][i];
			//v1-v2+v3-v4
			f[3][i] = xv[mapping[0]][i] - xv[mapping[1]][i] + xv[mapping[2]][i] - xv[mapping[3]][i];
		}
		newtonSolve2D(f, &u, &v);
		frac[mapping[0]] = (1 - u) * (1 - v);
		frac[mapping[1]] = u * (1 - v);
		frac[mapping[2]] = u * v;
		frac[mapping[3]] = (1 - u) * v;


	}
	break;
	default:
		printf("Interpolation not implemented for other types in 2D\n");
		break;
	}

}

void computeNodalWeight3D(double xv[8][3], double* xp, double frac[8], int nvert)
{
	int i, j, k, isolflag;
	double** lhs;
	double* rhs;
	double f[8][3];
	double u, v, w;
	double oneminusU, oneminusV, oneminusW, oneminusUV;

	switch (nvert)
	{
	case 4:
		//
		// tetrahedron
		//
		lhs = new double* [3];
		for (i = 0; i < 3; i++)
			lhs[i] = new double[3];
		rhs = new double[3];
		for (k = 0; k < 3; k++)
		{
			for (j = 0; j < 3; j++)
				lhs[j][k] = xv[k][j] - xv[3][j];
			rhs[k] = xp[k] - xv[3][k];
		}
		//
		// invert the 3x3 matrix
		//
		solvec(lhs, rhs, &isolflag, 3);
		//
		// check if the solution is not degenerate
		//
		if (isolflag)
		{
			for (k = 0; k < 3; k++) frac[k] = rhs[k];
			frac[3] = 1. - frac[0] - frac[1] - frac[2];
		}
		else
		{
			frac[0] = 1.0;
			frac[1] = frac[2] = frac[3] = 0;
		}
		for (i = 0; i < 3; i++) delete[]lhs[i];
		delete[] lhs;
		delete[] rhs;
		break;
	case 5:
		//
		// pyramid
		//
		for (j = 0; j < 3; j++)
		{
			f[0][j] = xv[0][j] - xp[j];
			f[1][j] = xv[1][j] - xv[0][j];
			f[2][j] = xv[3][j] - xv[0][j];
			f[3][j] = xv[4][j] - xv[0][j];
			//
			f[4][j] = xv[0][j] - xv[1][j] + xv[2][j] - xv[3][j];
			f[5][j] = xv[0][j] - xv[3][j];
			f[6][j] = xv[0][j] - xv[1][j];
			f[7][j] = -xv[0][j] + xv[1][j] - xv[2][j] + xv[3][j];
		}
		//
		newtonSolve3D(f, &u, &v, &w);
		oneminusU = 1.0 - u;
		oneminusV = 1.0 - v;
		oneminusW = 1.0 - w;
		//
		frac[0] = oneminusU * oneminusV * oneminusW;
		frac[1] = u * oneminusV * oneminusW;
		frac[2] = u * v * oneminusW;
		frac[3] = oneminusU * v * oneminusW;
		frac[4] = w;
		//
		break;
	case 6:
		//
		// prizm
		//
		for (j = 0; j < 3; j++)
		{
			f[0][j] = xv[0][j] - xp[j];
			f[1][j] = xv[1][j] - xv[0][j];
			f[2][j] = xv[2][j] - xv[0][j];
			f[3][j] = xv[3][j] - xv[0][j];
			//
			f[4][j] = 0;
			f[5][j] = xv[0][j] - xv[2][j] - xv[3][j] + xv[5][j];
			f[6][j] = xv[0][j] - xv[1][j] - xv[3][j] + xv[4][j];
			f[7][j] = 0.;
		}
		//
		newtonSolve3D(f, &u, &v, &w);
		//
		oneminusUV = 1.0 - u - v;
		oneminusU = 1.0 - u;
		oneminusV = 1.0 - v;
		oneminusW = 1.0 - w;
		//
		frac[0] = oneminusUV * oneminusW;
		frac[1] = u * oneminusW;
		frac[2] = v * oneminusW;
		frac[3] = oneminusUV * w;
		frac[4] = u * w;
		frac[5] = v * w;
		//
		break;
	case 8:
		//
		// hexahedra
		//
		for (j = 0; j < 3; j++)
		{
			f[0][j] = xv[0][j] - xp[j];
			f[1][j] = xv[1][j] - xv[0][j];
			f[2][j] = xv[3][j] - xv[0][j];
			f[3][j] = xv[4][j] - xv[0][j];
			//
			f[4][j] = xv[0][j] - xv[1][j] + xv[2][j] - xv[3][j];
			f[5][j] = xv[0][j] - xv[3][j] + xv[7][j] - xv[4][j];
			f[6][j] = xv[0][j] - xv[1][j] + xv[5][j] - xv[4][j];
			f[7][j] = -xv[0][j] + xv[1][j] - xv[2][j] + xv[3][j] +
				xv[4][j] - xv[5][j] + xv[6][j] - xv[7][j];
		}
		//
		newtonSolve3D(f, &u, &v, &w);
		//
		oneminusU = 1.0 - u;
		oneminusV = 1.0 - v;
		oneminusW = 1.0 - w;
		//
		frac[0] = oneminusU * oneminusV * oneminusW;
		frac[1] = u * oneminusV * oneminusW;
		frac[2] = u * v * oneminusW;
		frac[3] = oneminusU * v * oneminusW;
		frac[4] = oneminusU * oneminusV * w;
		frac[5] = u * oneminusV * w;
		frac[6] = u * v * w;
		frac[7] = oneminusU * v * w;
		//
		break;
	default:
		printf("Interpolation not implemented for polyhedra with %d vertices\n", nvert);
		break;
	}
}

void computeNodalWeights(int dim, double xv[8][3], double* xp, double* frac, int nvert)
{
	if (dim == 2)
	{
		computeNodalWeight2D(xv, xp, frac, nvert);
	}
	else if (dim == 3)
	{
		computeNodalWeight3D(xv, xp, frac, nvert);
	}
}


//Basic quickSort Algorithm, used to compare with the find_median method
int partition(double* A, int* iA, int p, int r)
{
	double x = A[r];
	int i = p - 1;
	double xtemp;
	int ixtemp;
	for (int j = p; j < r - 1; j++)
	{
		if (A[j] <= x)
		{
			i = i + 1;
			//exchange A[i] with A[j];
			xtemp = A[i];
			A[i] = A[j];
			A[j] = xtemp;

			ixtemp = iA[i];
			iA[i] = iA[j];
			iA[j] = ixtemp;
		}
	}
	//exchange A[i+1] with A[r];
	xtemp = A[i + 1];
	A[i + 1] = A[r];
	A[r] = xtemp;

	ixtemp = iA[i + 1];
	iA[i + 1] = iA[r];
	iA[r] = ixtemp;

	return i + 1;

}

void quick_sort(double* A, int* iA, int p, int r)
{
	if (p < r)
	{
		int q = partition(A, iA, p, r);
		quick_sort(A, iA, p, q - 1);
		quick_sort(A, iA, q + 1, r);
	}

}

void quick_sort_wrap(int* ix, double* x, int* n, double* xmed)
{
	quick_sort(x, ix, 0, *n - 1);
	if (*n % 2 == 1)//Odd
	{
		*xmed = x[(*n / 2)];
	}
	else
	{
		*xmed = 0.5 * (x[(*n / 2) - 1] + x[(*n / 2)]);
	}

}
