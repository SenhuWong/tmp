#pragma once
#include<iostream>

class Brick
{
protected:
	int nodes[8] = {-1,-1,-1,-1,-1,-1,-1,-1};
public:
	int size = 0;
	void setSize(int isize)
	{
		if(isize > 8 or isize <=0)
		{
			std::cout<<"Wrong size "<<isize<<"\n";
		}
		size = isize;
	}
	int sizeIs()
	{
		return size;
	}
	void setNode(int inode, int index)
	{
		if(index<0 or index >=8)
		{
			std::cout<<"INdex error: ind is "<<index<<" inode is"<<inode <<'\n';
		}
		else
		{

		nodes[index] = inode;
		}
	}
	int nodeIs(int index)
	{
		return nodes[index];
	}
	int* begin()
	{
		return nodes;
	}

	virtual void insert(int* start, double* coordinates, int nvert)
	{
		
	}


};

class Brick4Cobalt2D:public Brick
{
	bool closed = false;
public:
	Brick4Cobalt2D(){}
	Brick4Cobalt2D(const Brick4Cobalt2D& another)
	{
		size = another.size;
		closed = true;
		for (int i = 0; i < 4; i++)
		{
			nodes[i] = another.nodes[i];
		}
	}
	Brick4Cobalt2D& operator=(const Brick4Cobalt2D& another);
	void insert(int* start, double* coordinates, int nvert);
	bool checkRepeat();
};

class Brick4Cobalt3D:public Brick
{
	int nface = 0;
	int nnodes[6] = { -1,-1,-1,-1,-1,-1 };
	int eachNodesDir[6] = {0};
	int* eachNodes[6] = {nullptr};
public:
	Brick4Cobalt3D() {}
	Brick4Cobalt3D(const Brick4Cobalt3D& another)
	{
		for (int i = 0; i < 8; i++)
		{
			nodes[i] = another.nodes[i];
		}
		size = another.size;
	}
	Brick4Cobalt3D& operator=(const Brick4Cobalt3D& another)
	{
		for (int i = 0; i < 8; i++)
		{
			nodes[i] = another.nodes[i];
		}
		size = another.size;
		return *this;
	}

	void insert(int* start,int dir, double* coordinates,int nvert)//Totally could be 3 or 4 vert 
	{
		eachNodes[nface] = start;
		eachNodesDir[nface] = dir;//Dir is used to indicate if this face is pointing out of the cell (1/-1);
		nnodes[nface] = nvert;
		nface++;
	}
	bool checkRepeat();
	void reOrdering();
};