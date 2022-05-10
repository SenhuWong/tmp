#include <stdio.h>
#include <stdlib.h>
#include <time.h>


struct cell
{
	double x, y, z;
};

struct tri
{
	float normal[3];
	float p0[3];
	float p1[3];
	float p2[3];
	unsigned short source;
};

void get3minmax(double v[][3],double max[],double min[])
{
	int i;
	for(i=0;i<3;i++)
	{
		max[i]=min[i]=v[0][i];
		if(v[1][i]<min[i])   min[i]=v[1][i];
		if(v[1][i]>max[i])   max[i]=v[1][i];
		if(v[2][i]<min[i])   min[i]=v[2][i];
		if(v[2][i]>max[i])   max[i]=v[2][i];
	}
}

void getminmax(double p0, double p1,double *max,double *min)
{
	if(p0>p1)
	{
		*max=p0;
		*min=p1;
	}
	else
	{
		*max=p1;
		*min=p0;
	}
}

int IsSeparatingAxis(double p1,double p2,double r)
{
	double min,max;
	getminmax(p1,p2,&max,&min);
	//printf("max=%f,min=%f   r=%f\n",max,min,r);
	//getchar();
	if((max<-r)||(min>r))
		return 1;
	else 
		return 0;
}

int TestAxisCrossEdge(double v[][3],double e[],double h)
{
	double fex,fey,fez;
	double p0,p2;
	double r;

	
	fex=fabs(e[0]);
	fey=fabs(e[1]);
	fez=fabs(e[2]);
	p0=-e[2]*v[0][1]+e[1]*v[0][2];
	p2=-e[2]*v[2][1]+e[1]*v[2][2];
	r=fez*h+fey*h;
	//printf("r1=%f\n", r);
	if (IsSeparatingAxis(p0, p2, r)==1)
	{
        return 0;
	}

	p0=e[2]*v[0][0]-e[0]*v[0][2];
	p2=e[2]*v[2][0]-e[0]*v[2][2];
	r=fez*h+fex*h;

	//getchar();
	if (IsSeparatingAxis(p0, p2, r)==1)
	{
		return 0;
	}
	//printf("1\n");
	//getchar();
    p0=-e[1]*v[0][0]+e[0]*v[0][1];
	p2=-e[1]*v[2][0]+e[0]*v[2][1];
	r=fey*h+fex*h;
	//printf("r3=%f\n", r);
	if (IsSeparatingAxis(p0, p2, r)==1)
	{
		return 0;
	}

	return 1;
}

//0 separate 1 overlap
int intersect_judge(struct cell CartC,struct tri triangle)  //true false
{
	//CABA_point_t CartC；   //笛卡尔网格中点
	//CABA_point_t triP[3];  //待判定三角形的三个顶点
	int i,j;
	double h=0.5;   // 半边长
	double d;
	double t;
	double v[3][3];  //Vi-C
	double e0[3],e1[3],e2[3];
	double u[3][3];
	double max[3],min[3];
	double normal[3];
	double maxpoint[3],minpoint[3];
	
	v[0][0] = (double)triangle.p0[0] - CartC.x;
	v[1][0] = (double)triangle.p1[0] - CartC.x;
	v[2][0] = (double)triangle.p2[0] - CartC.x;
	v[0][1] = (double)triangle.p0[1] - CartC.y;
	v[1][1] = (double)triangle.p1[1] - CartC.y;
	v[2][1] = (double)triangle.p2[1] - CartC.y;
	v[0][2] = (double)triangle.p0[2] - CartC.z;
	v[1][2] = (double)triangle.p1[2] - CartC.z;
	v[2][2] = (double)triangle.p2[2] - CartC.z;

	/* 9 axis separation edge check*/
	for(i=0;i<3;i++)
	{
		e0[i]=v[1][i]-v[0][i]; //p1-p0
		e1[i]=v[2][i]-v[1][i];
		e2[i]=v[0][i]-v[2][i]; //  f0 f1 f2
	}
	if(TestAxisCrossEdge(v,e0,h)==0)   //section iii
	{
		return 0;
	}
	for(i=0;i<3;i++)
	{
		u[0][i]=v[1][i];
		u[1][i]=v[2][i];
		u[2][i]=v[0][i];
	}
	if (TestAxisCrossEdge(u, e1, h)==0)
	{
		return 0;
	}
	for(i=0;i<3;i++)
	{
		u[0][i]=v[2][i];
		u[1][i]=v[0][i];
		u[2][i]=v[1][i];
	}
	if (TestAxisCrossEdge(u, e2, h)==0)
	{
		return 0;
	}
	/* 3 axis range edge check*/
	get3minmax(v,max,min);
	if((min[0]>h)||(max[0]<-h))
	{
		return 0;
	}
	if((min[1]>h)||(max[1]<-h))
	{
		return 0;
	}
	if((min[2]>h)||(max[2]<-h))
	{
		return 0;
	}
	/* normal vector axis check*/
	normal[0]=(v[1][1]-v[0][1])*(v[2][2]-v[0][2])-(v[1][2]-v[0][2])*(v[2][1]-v[0][1]);
	normal[1]=(v[1][2]-v[0][2])*(v[2][0]-v[0][0])-(v[1][0]-v[0][0])*(v[2][2]-v[0][2]);
	normal[2]=(v[1][0]-v[0][0])*(v[2][1]-v[0][1])-(v[1][1]-v[0][1])*(v[2][0]-v[0][0]);
	d = 0.0;
	d=-1.0*(normal[0]*v[0][0]+normal[1]*v[0][1]+normal[2]*v[0][2]);
	for(i=0;i<3;i++)
	{
		if(normal[i]>1.0e-16)
		{
			minpoint[i]=-h;
			maxpoint[i]=h;
		}
		else
		{
			minpoint[i]=h;
			maxpoint[i]=-h;
		}	
	}
	t = 0.0;
	for(i=0;i<3;i++)
	{
		t=t+normal[i]*minpoint[i];
	}
	if(t+d>1.0e-16)
	{
		return 0;
	}
	t=0.0;
	for(i=0;i<3;i++)
	{
		t=t+normal[i]*maxpoint[i];
	}
	//printf("maxpoint[0]=%f maxpoint[1]=%f  maxpoint[2]=%f  normal[0]=%f normal[1]=%f norma[2]=%f t=%f d=%f\n", maxpoint[0],maxpoint[1],maxpoint[2],normal[0], normal[1], normal[2], t, d);
	if(t+d<=-(1.0e-16))
	{
		return 0;
	}
	//printf("6\n");
	//getchar();
	return 1;
	
}

 int main(void)
{
	 FILE *fout,*fin;
	int i ,j;
	int result;
	int num = 10000;
	double a = 2.0;
	struct cell cart;
	struct tri triangle;
	char out[20];
	char in[20];
	int size;
	cart.x = 0.5;
	cart.y = 0.5; 
	cart.z = 0.5;
	sprintf_s(out, sizeof(out), "data.txt");
	sprintf_s(in, sizeof(in), "selta.stl");
	fopen_s(&fout, out, "w");
	fopen_s(&fin, in, "rb");

	fseek(fin, 0, SEEK_END);
	size = ftell(fin);
	//printf("size=%d\n",size);
	//getchar();
	//rewind(fin);
	fseek(fin,80L,0);
	fread(&num, sizeof(int), 1, fin);
	printf("num=%d\n", num);
	int standardBinaryFileSize = 80 + sizeof(int)+num * 50;
	if (size == standardBinaryFileSize)
	{
		printf("standard stl file\n");
	}
	getchar();
	for (i = 0; i < num; i++)
	{
		fread(triangle.normal, sizeof(float), 3, fin);
		//fseek(fin, 12L, 1);
		//for (j = 0; j < 3; j++)
		{
			//fseek(fin, 12L, 1);
			
			fread(triangle.p0, sizeof(float), 3, fin);
			fread(triangle.p1, sizeof(float), 3, fin);
			fread(triangle.p2, sizeof(float), 3, fin);
			//printf("triangle.x=%f y=%f z=%f\n", triangle.p0[0],triangle.p0[1],triangle.p0[2]);
	
		}
		fread(&(triangle.source), sizeof(unsigned short), 1, fin);
		//getchar();
		
		result = intersect_judge(cart, triangle);
		fprintf(fout, "%lf   %lf   %lf   %lf   %lf   %lf   %lf   %lf   %lf   %d\n", triangle.p0[0], triangle.p0[1], triangle.p0[2], triangle.p1[0], triangle.p1[1], triangle.p1[2], triangle.p2[0], triangle.p2[1], triangle.p2[2], result);
	}
	fclose(fout);
	fclose(fin);
	getchar();
	return 0;
}