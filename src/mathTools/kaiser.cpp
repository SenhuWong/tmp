#include<math.h>
#include<iostream>
const double small_value = 1.0e-12;
//const double zero = 0.0;
const double half = 0.5;
const double one = 1.0;
#define TESTING true

void kaiser(double** a, int nrows, int n, double* eigenv, double* trace, double* sume, int& ier)
//a is the (nrows,n) matrix
//ier is:
// 0 no error;
// 1 n>nrows or n<1
// 2 Failed to converge within 10 iteration(for such a 2*2,a simple rotation should be fine)
{
	if (n<1 or n>nrows)
	{
		ier = 1;
		return;
	}
	ier = 0;
	//Get tolerance and trace
	*trace = 0;
	double eps;
	double ss = 0;
	for (int i = 0; i < n; i++)
	{
		*trace += a[i][i];
		for (int j = 0; j < n; j++)
		{
			ss += a[i][j] * a[i][j];
		}
	}
	eps = small_value * ss / n;
	//std::cout << "EPS is " << eps << "ss is " << ss << '\n';
	double halfP;// Xi*Yi
	double p;
	double absP;
	double q;//Xi*Xi - Yi*Yi
	double absQ;
	//This is a single rotation
	int iter_count;
	int max_iter_count = 100;
	int skip_count; 
	double TAN;
	double COS;
	double SIN;
	double CTN;
	double temp;
	//std::cout << "Matrix is:\n";
	//std::cout << a[0][0] << "\t" << a[0][1] << '\n';
	//std::cout << a[1][0] << "\t" << a[1][1] << '\n';
	for (iter_count = 0; iter_count < max_iter_count; iter_count++)
	{
		skip_count = 0;
		for (int i = 0; i < n - 1; i++)
		{
			for (int j = i + 1; j < n; j++)
			{

				//Calculate the rotation
				halfP = 0;
				q = 0;
				for (int k = 0; k < n; k++)
				{
					halfP += a[k][i] * a[k][j];
					q += (a[k][i] + a[k][j]) * (a[k][i] - a[k][j]);
				}


				//
				p = 2 * halfP;
				absP = abs(p);
				//If p is very small, then vectors are almost orthogonal
				//If q > 0, correct ordering then skip the rotation
				if (absP < eps and q >= 0)
				{
					//std::cout << "absP is " << absP << " and q is " << q << '\n';
					skip_count += 1;
					continue;
				}
				//Rotation needed ,now use the formula
				absQ = abs(q);
				if (absP <= absQ) //tan(2a) <= 1;
				{
					TAN = absP / absQ;//get cos(2a) and sin(2a)
					COS = 1 / sqrt(1 + TAN * TAN);
					SIN = TAN * COS;
				}
				else
				{
					CTN = absQ / absP;
					SIN = 1 / sqrt(1 + CTN * CTN);
					COS = SIN * CTN;
				}
				COS = sqrt((1 + COS) / 2);//get cos(a) and sin(a)
				SIN = SIN / (2 * COS);
				//Value of p and q decides the final rotation
				if (q < 0)
				{
					temp = COS;
					COS = SIN;
					SIN = temp;
				}
				if (p < 0)
				{
					SIN = -SIN;
				}
				//Do the rotation
				for (int k = 0; k < n; k++)
				{
					temp = a[k][i];
					a[k][i] = temp * COS + a[k][j] * SIN;
					a[k][j] = -temp * SIN + a[k][j] * COS;
				}
			}
		}
		//Check if we stop the rotation 
		if (skip_count == 0.5* n * (n - 1))
		{
			//std::cout << "Breaking out of this at iter " << iter_count << '\n';
			break;//No more need for any left iteration
		}
		else
		{
			//std::cout << "Skip count is:" << skip_count<<" at iter "<< iter_count << '\n';
		}
	}
	if (iter_count == max_iter_count)
	{
		ier = 2;
	}
	*sume = 0;
	//The ith column's eigenvalue
	for (int i = 0; i < n; i++)
	{
		temp = 0;
		for (int j = 0; j < n; j++)
		{
			temp += a[j][i] * a[j][i];
		}
		eigenv[i] = sqrt(temp);
		*sume += eigenv[i];
	}
	for (int i = 0; i < n; i++)
	{
		if (eigenv[i] > 0)
		{
			for (int j = 0; j < n; j++)
			{
				a[j][i] = a[j][i] / eigenv[i];
			}
		}
		else
		{
			for (int j = 0; j < n; j++)
			{
				a[j][i] = 0;
			}
		}
	}

}

void kaiser_wrap(double** matrixA, const int& nrows,const int& n , double* eigenv, double* trace, double* sume, int* ier)
{
	kaiser(matrixA, nrows, n, eigenv, trace, sume, *ier);
}