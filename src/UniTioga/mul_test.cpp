#include<random>
#include<iomanip>
#include<string>
#include<fstream>
#include<iostream>
#include<set>
const int ODDVOL = 5;
const int EVENVOL = 8;
const double TOL = 1.0e-10;
#include"tioga.h"
//#include"MeshBlock.h"
void find_median_wrap(int* ix, double* x, int n);
//%TODO::The trouble of using [8][3] here only relies on the testing end, so if in later code this could be modified to be dimension and nvert specific pointer
double computeVolume(int dim, double xv[8][3], int nvert, bool counter_clock);
void computeNodalWeight2D(double xv[8][3], double* xp, double* frac, int nvert);
static int fail_count = 0;
void findOBB(int dim, double* x, double xc[3], double dxc[3], double vec[3][3], int nnodes);
void output(int* ix, double* x, int n)
{
	for (int i = 0; i < n; i++)
	{
		std::cout << std::setw(10) << ix[i];
	}
	std::cout << '\n';
	for (int i = 0; i < n; i++)
	{
		std::cout << std::setw(10) << x[i];
	}
	std::cout << '\n';
}

void assertion(int* ix, double* x, int n)
{
	fail_count = 0;
	int imid;
	if (n % 2 == 0)
	{
		imid = n / 2 - 1;
	}
	else
	{
		imid = n / 2;
	}
	bool valid = true;
	for (int i = 0;i < imid; i++)
	{
		if (x[i] > x[imid])
		{
			fail_count++;
			valid = false;
			//output(ix, x, n);
		}
	}
	for (int i = imid + 1; valid and i < n; i++)
	{
		if (x[i] < x[imid])
		{
			fail_count++;
			valid = false;
			//output(ix, x, n);
		}
	}
	if (valid == false)
	{
		fail_count += 1;
		std::cout << "Fail count is: " << fail_count << '\n';
		output(ix, x, n);
	}

}

void checkEven(int times,double* xEven, int *indexEven,std::default_random_engine& engine,std::uniform_real_distribution<double>& distr)
{
	for (int k = 0; k < times; k++)
	{
		for (int i = 0; i < EVENVOL; i++)
		{
			xEven[i] = distr(engine);
		}
		find_median_wrap(indexEven, xEven, EVENVOL);
		assertion(indexEven, xEven, EVENVOL);
	}
}
void checkOdd(int times, double* xOdd, int* indexOdd, std::default_random_engine& engine, std::uniform_real_distribution<double>& distr)
{
	for (int k = 0; k < times; k++)
	{
		for (int i = 0; i < ODDVOL; i++)
		{
			xOdd[i] = distr(engine);
		}
		find_median_wrap(indexOdd, xOdd, ODDVOL);
		assertion(indexOdd, xOdd, ODDVOL);
	}

}
void median_test()
{
	std::random_device rd;
	std::default_random_engine engine(rd());
	std::uniform_real_distribution<double> distr(100, 200);
	int indexEven[EVENVOL];
	int indexOdd[ODDVOL];
	double xEven[EVENVOL];
	double xOdd[ODDVOL];
	for (int i = 0; i < EVENVOL; i++)
	{
		indexEven[i] = i;
	}
	for (int i = 0; i < ODDVOL; i++)
	{
		indexOdd[i] = i;
	}
	checkEven(1000, xEven, indexEven, engine, distr);
	std::cout << "End of checking Even(Enter)\n";
	std::cin.get();
	checkOdd(1000, xOdd, indexOdd, engine, distr);
	std::cout << "End of checking Odd(Enter)\n";
	std::cin.get();

}

void volume_test()
{
	double triangle_coordinates[3][3] = 
	{ 
		0,1,0,
		1,1,0,
		1,0,0 
	};
	double result = computeVolume(2, triangle_coordinates, 3, false);
	std::cout << "Volume is of triangle is:\n" << result << "\n";
	double ABCD[4][3] =
	{
		0,0,0,
		1,0,0,
		2,3,0,
		0,3,0
	};
	double resultABCD = computeVolume(2, ABCD, 4, false);
	double ACBD[4][3] =
	{
		0,0,0,
		2,3,0,
		1,0,0,
		0,3,0
	};
	double resultACBD = computeVolume(2, ACBD, 4, false);
	double ADCB[4][3] =
	{
		0,0,0,
		0,3,0,
		2,3,0,
		1,0,0,
	};
	double resultADCB = computeVolume(2, ADCB, 4, false);
	std::cout << "Volume of quadraliteral is :\n";
	std::cout << "resultABCD\t" << resultABCD << '\n';
	std::cout << "resultACBD\t" << resultACBD << '\n';
	std::cout << "resultADCB\t" << resultADCB << '\n';
}

#include"utils/vector3d.h"
bool real_TriangleInclusionTest(double *vp, double xv[3][3])
{
	//std::cout << "Real checking\n";
	Vector3D v1 = Vector3D(xv[0][0] - xv[1][0], xv[0][1] - xv[1][1], 0);
	Vector3D vp1 = Vector3D(vp[0] - xv[1][0], vp[1] - xv[1][1], 0);
	Vector3D v2 = Vector3D(xv[1][0] - xv[2][0], xv[1][1] - xv[2][1], 0);
	Vector3D vp2 = Vector3D(vp[0] - xv[2][0], vp[1] - xv[2][1], 0);
	Vector3D v3 = Vector3D(xv[2][0] - xv[0][0], xv[2][1] - xv[0][1], 0);
	Vector3D vp3 = Vector3D(vp[0] - xv[0][0], vp[1] - xv[0][1], 0);
	
	Vector3D v1xvp = v1.cross_product(vp1);
	Vector3D v2xvp = v2.cross_product(vp2);
	Vector3D v3xvp = v3.cross_product(vp3);
	if (v1xvp[2] <= 0 and v2xvp[2] <= 0 and v3xvp[2] <= 0)
	{
		return true;
	}
	if (v1xvp[2] >= 0 and v2xvp[2] >= 0 and v3xvp[2] >= 0)
	{
		return true;
	}
	return false;
}

//
void triangleTest(const std::string& filename,int nvert, int times, double coord[4][3], std::default_random_engine& engine, std::uniform_real_distribution<double>& distr)
{
	std::cout << "Testing Barycentric Coordinates for " << filename << ", vert count is " << nvert << "\n";
	std::ofstream fout;
	fout.open(filename);
	if (!fout.is_open())
	{
		std::cout << "Error at opening " << filename << '\n';
		return;
	}
	double* vp = new double[2];
	double xv[8][3];
	
	double* frac = new double[nvert];

	for (int i = 0; i < nvert; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			xv[i][j] = coord[i][j];
			fout << xv[i][j] << '\t';
		}
		fout << '\n';
	}
	//examine the frac of nodes
	for (int i = 0; i < nvert; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			vp[j] = xv[i][j];
		}
		computeNodalWeight2D(xv, vp, frac, nvert);
		std::cout << i << ": ";
		for (int j = 0; j < nvert; j++)
		{
			std::cout << frac[j] << '\t';
		}
		std::cout << '\n';
	}

	bool inner = true;
	double recomputed_vp[2];
	for (int i = 0; i < times; i++)
	{
		inner = true;
		recomputed_vp[0] = recomputed_vp[1] = 0;
		vp[0] = distr(engine);
		vp[1] = distr(engine);
		computeNodalWeight2D(xv, vp, frac, nvert);
		fout << vp[0] << '\t' << vp[1] << '\t';
		for (int j = 0; j < nvert; j++)
		{
			inner = inner and frac[j] >= 0 and frac[j] <= 1;
			recomputed_vp[0] += frac[j] * xv[j][0];
			recomputed_vp[1] += frac[j] * xv[j][1];
			fout << frac[j] << '\t';
		}
		bool real_inner = real_TriangleInclusionTest(vp, xv);
		if (inner == real_inner)
		{

		}
		else
		{
			throw std::runtime_error("Not right inclusion test\n");
		}
		if (inner)
		{
			if ((abs(recomputed_vp[0] - vp[0]) > TOL) or (abs(recomputed_vp[1] - vp[1]) > TOL))
			{
				std::cout << filename << "Not correspondant ,delta is " << recomputed_vp[0] - vp[0] << "," << recomputed_vp[1] - vp[1] << " \n";
			}

			fout << '\n';
		}
	}
	fout.close();
	return;
	
}
void barycentric_test()
{
	std::random_device rd;
	std::default_random_engine engine(rd());
	std::uniform_real_distribution<double> distr(-4, 4);
	double triangle_coordinates[3][3] =
	{
		0,1,0,
		1,1,0,
		1,0,0
	};
	triangleTest("triangle", 3, 1000, triangle_coordinates, engine, distr);
	double triangle_coordinates2[3][3] =
	{
		0,1,0,
		1,0,0,
		1,1,0
	};
	triangleTest("Inverse", 3, 1000, triangle_coordinates2, engine, distr);
	return;
	double ABCD[4][3] =
	{
		0,0,0,
		1,0,0,
		2,3,0,
		0,3,0
	};
	triangleTest("ABCD", 4, 1000, ABCD, engine, distr);


	double ACBD[4][3] =
	{
		0,0,0,
		2,3,0,
		1,0,0,
		0,3,0
	};
	triangleTest("ACBD", 4, 1000, ACBD, engine, distr);
	double ADCB[4][3] =
	{
		0,0,0,
		0,3,0,
		2,3,0,
		1,0,0,
	};
	triangleTest("ADCB", 4, 1000, ADCB, engine, distr);
	double ABDC[4][3] =
	{
		0,0,0,
		1,0,0,
		0,3,0,
		2,3,0
	};
	triangleTest("ABDC", 4, 1000, ABDC, engine, distr);
	double ACDB[4][3] =
	{
		0,0,0,
		2,3,0,
		0,3,0,
		1,0,0
	};
	triangleTest("ACDB", 4, 1000, ACDB, engine, distr);
	double ADBC[4][3] =
	{
		0,0,0,
		0,3,0,
		1,0,0,
		2,3,0
	};
	triangleTest("ADBC", 4, 1000, ADBC, engine, distr);
	

}

void readGrdTest(const std::string& filename)
{
	int dimension = 2;
	int nnodes;
	int nedges;
	int ncells;
	std::ifstream fin;
	fin.open(filename);
	if(!fin.is_open())
	{
		throw std::runtime_error("Failure at readGrdTest for opening file.");
		return;
	}
	fin >> nnodes >> nedges >> ncells;
	double* x = new double[dimension * nnodes];
	int* ibl = new int[nnodes];
	int temp[4];
	for (int i = 0; i < nnodes; i++)
	{
		ibl[i] = 0;
		for (int j = 0; j < dimension; j++)
		{
			fin >> x[dimension * i + j];
		}
	}
	std::set<int, std::less<int>> wbc;
	std::set<int, std::less<int>> obc;
	
	//I only need to know which are wbc and which be obc, so the connection is not what i  would care
	for (int i = 0; i < nedges; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			fin >> temp[j];
		}
		if (temp[3] == -1)//wallboundary
		{
			wbc.insert(temp[0]);
			wbc.insert(temp[1]);

		}
		else if (temp[3] == -2)//overboundary
		{
			obc.insert(temp[0]);
			obc.insert(temp[1]);
		}
	}
	//%TODO::Here is a special case where only triangle exists
	int** vconn = new int* [1];
	vconn[0] = new int[3*ncells];
	int nv[1] = { 3 };
	int nc[1] = { ncells };
	for (int i = 0; i < ncells; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			fin >> vconn[0][3 * i + j];
		}
	}
	int nwbc = wbc.size();
	int nobc = obc.size();
	int* wbc_array = new int[nwbc];
	int* obc_array = new int[nobc];
	auto iter = wbc.begin();
	for (int k =0;k<nwbc;k++)
	{
		wbc_array[k] = *iter;
		iter++;
	}
	iter = obc.begin();
	for (int k = 0; k < nobc; k++)
	{
		obc_array[k] = *iter;
		iter++;
	}
	MeshBlock* mb = new MeshBlock[1];
	mb->setData(2,1,nnodes,x,ibl,nwbc,nobc,wbc_array,obc_array,1,nv,nc,vconn);
	mb->preprocess();
	delete[] mb;
	return;
	mb->writeOBB("OBB");
	mb->writeGridFile("whatever");
}

void kaiser_test(int dim,int times)
{
	//2d only
	std::default_random_engine engine;
	std::uniform_real_distribution<double> ten2twen(10, 20);
	std::uniform_real_distribution<double> one2two(1, 2);
	double* xs = new double[2 * times];
	double xmax[2] = {-1000000,-10000000};
	double xmin[2] = {1000000,1000000};
	for (int i = 0; i < times; i++)
	{
		xs[2 * i] = ten2twen(engine);
		xs[2 * i + 1] = 2 * xs[2 * i] + one2two(engine);
		xmax[0] = TIOGA_MAX(xmax[0], xs[2 * i]);
		xmax[1] = TIOGA_MAX(xmax[1], xs[2 * i + 1]);
		xmin[0] = TIOGA_MIN(xmin[0], xs[2 * i]);
		xmin[1] = TIOGA_MIN(xmin[1], xs[2 * i + 1]);
	}
	std::cout << "xmax " << xmax[0] << '\t' << xmax[1] << '\n';
	std::cout << "xmin " << xmin[0] << '\t' << xmin[1] << '\n';
	OBB* obb = new OBB[1];
	findOBB(2, xs, obb->xc, obb->dxc, obb->vec, times);
	MeshBlock* mb = new MeshBlock[1];
	mb->writeOBB2("random", obb);
}
using namespace TIOGA;

void tg_readGrdTest(const std::string& filename,int cur_proc,tioga*tg)
{

	int dimension = 2;
	int nnodes;
	int nedges;
	int ncells;
	std::ifstream fin;

	fin.open(filename);
	if (!fin.is_open())
	{
		throw std::runtime_error("Failure at readGrdTest for opening file.");
		return;
	}
	fin >> nnodes >> nedges >> ncells;
	double* x = new double[dimension * nnodes];
	int* ibl = new int[nnodes];
	int temp[4];
	for (int i = 0; i < nnodes; i++)
	{
		ibl[i] = 0;
		for (int j = 0; j < dimension; j++)
		{
			fin >> x[dimension * i + j];
		}
	}
	if (cur_proc == 0)
	{
		double temp;
		for (int i = 0; i < nnodes; i++)
		{
			temp = x[2 * i];
			x[2 * i] = temp * sqrt(3) / 2 + x[2 * i + 1] / 2;
			x[2 * i + 1] = temp * (-0.5) + x[2 * i + 1] * sqrt(3) / 2;

		}
		for (int i = 0; i < dimension * nnodes; i++)
		{
			x[i] += 6;
		}
	}
	std::set<int, std::less<int>> wbc;
	std::set<int, std::less<int>> obc;

	//I only need to know which are wbc and which be obc, so the connection is not what i  would care
	for (int i = 0; i < nedges; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			fin >> temp[j];
		}
		if (temp[3] == -1)//wallboundary
		{
			wbc.insert(temp[0]);
			wbc.insert(temp[1]);

		}
		else if (temp[3] == -2)//overboundary
		{
			obc.insert(temp[0]);
			obc.insert(temp[1]);
		}
	}
	//%TODO::Here is a special case where only triangle exists
	int** vconn = new int* [1];
	vconn[0] = new int[3 * ncells];
	int* nv = new int[1];
	nv[0] = 3;
	int* nc = new int[1];
	nc[0] = ncells;

	for (int i = 0; i < ncells; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			fin >> vconn[0][3 * i + j];
		}
	}
	int nwbc = wbc.size();
	int nobc = obc.size();
	int* wbc_array = new int[nwbc];
	int* obc_array = new int[nobc];
	auto iter = wbc.begin();
	for (int k = 0; k < nwbc; k++)
	{
		wbc_array[k] = *iter;
		iter++;
	}
	iter = obc.begin();
	for (int k = 0; k < nobc; k++)
	{
		obc_array[k] = *iter;
		iter++;
	}
	//MeshBlock* mb = new MeshBlock[1];
	tg->registerGridData(cur_proc+1, nnodes, x, ibl, nwbc, nobc, wbc_array, obc_array, 1, nv, nc, vconn);
	std::cout << "GridData registered\n";
	//delete[] mb;
	return;
	//mb->writeOBB("OBB");
	//mb->writeGridFile("whatever");
}


void random_baryCentric_test(const std::string& filename,bool random = true)
{
	int times = 1000000;
	int dim = 2;
	int nfrac = 3;

	std::ofstream fout;
	fout.open(filename);
	if (!fout.is_open())
	{
		throw std::runtime_error("Can't open file\n");
		return;
	}
	std::default_random_engine engine;
	std::uniform_real_distribution<double> center(-2,2);
	std::uniform_real_distribution<double> range(-8, 8);
	std::uniform_real_distribution<double> scatter(-10, 10);
	double xv[3][3];
	xv[0][0] = -7.5;
	xv[0][1] = -5;
	xv[1][0] = 2.5;
	xv[1][1] = -2.5;
	xv[2][0] = -3;
	xv[2][1] = 6;
	double xp[2];
	double frac[3];
	if (random)
	{
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				xv[i][j] = center(engine) + range(engine);
			}
		}
	}
	fout << times << " " << dim << " " << nfrac << '\n';
	for (int i = 0; i < 3; i++)
	{
		fout << xv[i][0] << " " << xv[i][1] << '\n';
	}
	for (int i = 0; i < times; i++)
	{
		xp[1] = scatter(engine);
		xp[0] = scatter(engine);
		computeNodalWeight2D(xv, xp, frac, 3);
		bool valid = true;
		for (int j = 0; j < 3; j++)
		{
			if (frac[j] < 0 or frac[j]>1)
			{
				valid = false;
			}
		}
		if (!valid) continue;
		fout << xp[0] << " " << xp[1] << " ";
		for (int j = 0; j < 3; j++)
		{
			fout << frac[j] << " ";
		}
		fout << '\n';
	}
	fout.close();
	



}


void volume_test3D()
{
	//Tetrahedron
	double tetra_coord[4][3] =
	{
		0,0,0,
		1,0,0,
		1,2,0,
		0,1,1
	};
	double tetra_area =computeVolume(3, tetra_coord, 4, false);
	std::cout << "tetra vol " << tetra_area << '\n';
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			tetra_coord[i][j] += 10;
		}
	}
	tetra_area = computeVolume(3, tetra_coord, 4, false);
	std::cout << "tetra vol " << tetra_area << '\n';
	//pyramid
	double pyramid_coord[5][3] =
	{
		0,0,0,
		1,0,0,
		1,2,0,
		0,1,0,
		0.5,0.5,1
	};
	double pyramid_area = computeVolume(3, pyramid_coord, 5, false);
	std::cout << "pyramid_area " << pyramid_area << '\n';
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			pyramid_coord[i][j] += 5;
		}
	}
	pyramid_area = computeVolume(3, pyramid_coord, 5, false);
	std::cout << "pyramid_area " << pyramid_area << '\n';
	//prism
	double prism_coord[6][3] =
	{
		0,0,0,
		1,0,0,
		1,2,0,
		0.5,0,1,
		1,0,1,
		1,1,1
	};
	double prism_area = computeVolume(3, prism_coord, 6, false);
	std::cout << "prism area " << prism_area << '\n';
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			prism_coord[i][j] += 5;
		}
	}
	pyramid_area = computeVolume(3, prism_coord, 6, false);
	std::cout << "prism area " << pyramid_area << '\n';
	//hexahedra
	double hexahedra_coord[8][3] =
	{
		0,0,0,
		1,0,0,
		1,1,0,
		0,1,0,
		0,0,1,
		1,0,1,
		1,2,1,
		0,1,1
	};
	double hexahedra_coord0[8][3] =
	{
		0,0,0,
		1,0,0,
		1,1,0,
		0,1,0,
		0,0,1,
		1,0,1,
		1,1,1,
		0,1,1
	};
	double hexahedra_coord1[5][3] =
	{
		1,1,0,
		0,1,0,
		0,1,1,
		1,1,1,
		1,2,1
	};
	//double decompose = computeVolume(3, hexahedra_coord0, 8, false)+computeVolume(3, hexahedra_coord1, 5, false);
	double hexah_area = computeVolume(3, hexahedra_coord, 8, false);
	std::cout << "hexah_area " << hexah_area <<"while decompose is " << '\n';
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			hexahedra_coord[i][j] += 5;
		}
	}
	hexah_area = computeVolume(3, hexahedra_coord, 8, false);
	std::cout << "hexah_area " << hexah_area << '\n';
}

//int main()
//{
//	volume_test3D();
//}
//void test_iso_search()
//{
//	tioga* tg = new tioga[1];
//	tg_readGrdTest("naca0012.grd", 0, tg);
//	std::ifstream fin;
//	fin.open("Queries.dat");
//	int nsearches;
//	fin >> nsearches;
//	double* xsearches = new double[nsearches * 2];
//	for (int i = 0; i < nsearches; i++)
//	{
//		fin >> xsearches[2 * i] >> xsearches[2 * i + 1];	
//	}
//	tg->test_from_search(xsearches, nsearches);
//
//}
//
//int main()
//{
//	int mexclude = 1;
//	//test_iso_search();
//	MPI_Init(NULL, NULL);
//	int cur_proc, num_proc;
//	MPI_Comm_rank(MPI_COMM_WORLD, &cur_proc);
//	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
//	////
//	tioga* tg = new tioga[1];
//	tg->setCommunicator(MPI_COMM_WORLD, cur_proc, num_proc);
//	////
//	tg_readGrdTest("naca0012.grd", cur_proc, tg);
//	//
//	////
//	tg->setMexclude(&mexclude);
//	tg->profile();
//	tg->performConnectivity();
//	
//	delete[] tg;
//	MPI_Finalize();
//	//median_test();
//	//volume_test();
//	//random_baryCentric_test("random.dat",false);
//	//barycentric_test();
//	//readGrdTest("naca0012.grd");
//	//kaiser_test(2,1000);
//	//testTioga_exchangeBoxes();
//}