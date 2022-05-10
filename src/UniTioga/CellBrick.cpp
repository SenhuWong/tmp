#include"CellBrick.h"

Brick4Cobalt2D& Brick4Cobalt2D::operator=(const Brick4Cobalt2D& another)
{
	for (int i = 0; i < 4; i++)
	{
		nodes[i] = another.nodes[i];
	}
	size = another.size;
	return *this;
}

void Brick4Cobalt2D::insert(int* start, double* cellId, int nvert = 2)//In 2D edge is a line, nvert is set to 2 by default
{
	int cellId_int = *(int*)(cellId);
	//std::cout<<'\t'<<(cellId_int)<<'\n';
	bool recording = false;
	// if(cellId_int==16)
	// {
	// 	recording = true;
	// }
	//std::cout << "size is " << this->sizeIs() << '\t' << "closed ? " << closed << '\n';
	if (closed)//If already closed, do nothing
	{
		if(recording)
		{
			std::cout<<"Already filled\n";
		}
		return;
	}
	else if (size == 0)//If it is initialized, copy into first 2.
	{
		for (int i = 0; i < 2; i++)
		{
			nodes[i] = start[i];
		}
		if(recording)
		{
			std::cout<<"size==0, now become:\n";
			for(int i = 0;i<4;i++)
			{
				std::cout<<nodes[i]<<'\t';
			}
			std::cout<<'\n';
		}
		size = 2;
		return;
	}
	//There are already 2 or 3 vertices inside
	else
	{
		int indexes[4];//Common nodes(if any)
		int k = 0;
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < nvert; j++)
			{
				if (nodes[i] == start[j])
				{
					indexes[k++] = i;//i is the common node index in size
					indexes[k++] = j;//j is the common node index in start
				}
			}
		}
		if (k == 0)//No common vertice found, they are two 'parallel' line of a quadrilateral
		{
			if(recording)
			{
				std::cout<<"A parallel edge is added ,do nothing\n";
				for(int i = 0;i<4;i++)
				{
					std::cout<<nodes[i]<<'\t';
				}
				std::cout<<'\n';
			}
			//Do nothing
			/*for (int i = 0; i < 2; i++)
			{
				nodes[i + size] = start[1 - i];
			}
			size += 2;
			closed = true;*/
		}
		else if (k == 2)//Only one vertice in common
		{
			if(indexes[0]<2)//At the first added side.we could directly use it
			{
				nodes[3-indexes[0]] = start[1-indexes[1]];
			}
			else//At the later added side, we could directly
			{
				nodes[5-indexes[0]] = start[1-indexes[1]];

			}
			size += 1;
			if(recording)
			{
				std::cout<<"One vertice in common\n";
				std::cout<<indexes[0]<<" "<<indexes[1]<<'\n';
				for(int i = 0;i<4;i++)
				{
					std::cout<<nodes[i]<<'\t';
				}
				std::cout<<'\n';
			}
			if (size == 4)
			{
				closed = true;
				return;
			}

		}
		else
		{
			if (size == 3)
			{
				if (nodes[2] == -1 and nodes[3] != -1)
				{
					nodes[2] = nodes[3];
					nodes[3] = -1;
				}
			}
			if(recording)
			{
				std::cout<<"A triangle is found ,make it tight \n";
				for(int i = 0;i<4;i++)
				{
					std::cout<<nodes[i]<<'\t';
				}
				std::cout<<'\n';
			}
			closed = true;
		}

	}
}


bool Brick4Cobalt2D::checkRepeat()
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if (nodes[i] == nodes[j] and i != j)
			{
				return true;
			}
		}
	}
	return false;
}


bool Brick4Cobalt3D::checkRepeat()
{
	for (int i = 0; i < 8; i++)
	{
		for (int j = i + 1; j < 8; j++)
		{
			if (nodes[i] == nodes[j] and nodes[i]!=-1)
			{
				std::cout <<"size is "<<size<<" "<< nodes[i] << '\t' << nodes[j] << "are the same\n";
				return true;
			}
		}
	}
	return false;
}

void Brick4Cobalt3D::reOrdering()
{
	int nTri = 0;
	int nQuadri = 0;
	for (int i = 0; i < nface; i++)
	{
		if (nnodes[i] == 3)
		{
			nTri++;
		}
		else if (nnodes[i] == 4)
		{
			nQuadri++;
		}
		else
		{
			std::cout << "Error at reordering: ";
			throw std::runtime_error("face is not a triangle or quadrilateral\n");
		}
	}
	//(4,0) tetra
//(4,1) pyramid
//(2,3) prism
//(0,6) brick
	if (nTri == 4 and nQuadri == 0)
	{
		// std::cout << "TetraHedron\n";
		size = 4;
		for (int i = 0; i < 3; i++)
		{
			nodes[i] = eachNodes[0][i];
		}
		for (int i = 0; i < 3; i++)
		{
			nodes[size - 1] = eachNodes[1][i];
			int sameCount = 0;
			for (int j = 0; j < 3; j++)
			{
				if (nodes[size - 1] == nodes[j])
				{
					sameCount++;
					break;
				}
			}
			if (sameCount == 0)//There is no same vert
			{
				break;
			}
		}
	}
	else if (nTri == 4 and nQuadri == 1)//Pyramid
	{
		// std::cout << "Pyramid\n";
		size = 5;
		int bottom = 0;
		for (int i = 0; i < 5; i++)
		{
			if (nnodes[i] == 4)
			{
				bottom = i;
				break;
			}
		}
		for (int i = 0; i < 4; i++)
		{
			nodes[i] = eachNodes[bottom][i];
		}
		for (int i = 0; i < 3; i++)
		{
			nodes[size - 1] = eachNodes[(bottom + 1) % nface][i];
			int sameCount = 0;
			for (int j = 0; j < 4; j++)
			{
				if (nodes[size - 1] == nodes[j])
				{
					sameCount++;
					break;
				}
			}
			if (sameCount == 0)
			{
				break;
			}
		}
	}
	else if (nTri == 2 and nQuadri == 3)
	{
		// std::cout<<"Column\n";
		size = 6;
		int temp_count_tri = 0;
		int temp_count_quad = 0;
		int QuadIndex[3];
		int TriIndex[2];
		for (int i = 0; i < nface; i++)
		{
			if (nnodes[i] == 3)
			{
				TriIndex[temp_count_tri++] = i;
			}
			else if (nnodes[i] == 4)
			{
				QuadIndex[temp_count_quad++] = i;
			}
			else
			{
				std::cout << "Something wrong\n";
				throw std::runtime_error("Something wrong\n");
			}
		}
		//Fill in the first 3
		for (int i = 0; i < 3; i++)
		{
			nodes[i] = eachNodes[TriIndex[0]][i];
		}
		//Fill in the next 3 with Quadrilateral
		int indexes[4];
		int mappedindexes[2];
		int count = 0;
		for (int i = 0; i < 2; i++)//Each face
		{
			count = 0;
			for (int j = 0; j < 4; j++)//Each Vert in the face
			{
				int cur_index = eachNodes[QuadIndex[i]][j];
				for (int k = 0; k < 3; k++)//Each Vert in the first 3
				{
					if (nodes[k] == cur_index)
					{
						indexes[count++] = k;
						indexes[count++] = j;
					}
				}
			}
			if (count != 4)
			{
				std::cout << "Something wrong count is not 4\n";
			}
			int mini, maxi;
			mini = std::min(indexes[1], indexes[3]);
			maxi = std::max(indexes[1], indexes[3]);
			if (mini == 0 and maxi == 1)
			{
				mappedindexes[0] = 3 - indexes[1];
				mappedindexes[1] = 3 - indexes[3];
			}
			else if (mini == 2 and maxi == 3)
			{
				mappedindexes[0] = 3 - indexes[1];
				mappedindexes[1] = 3 - indexes[3];
			}
			else if (mini == 0 and maxi == 3)
			{
				mappedindexes[0] = indexes[1] == 0 ? 1 : 2;
				mappedindexes[1] = indexes[3] == 0 ? 1 : 2;
			}
			else if (mini == 1 and maxi == 2)
			{
				mappedindexes[0] = indexes[1] == 1 ? 0 : 3;
				mappedindexes[1] = indexes[3] == 1 ? 0 : 3;
			}
			else
			{
				std::cout << "Something I don't know happens\n";
			}
			for (int l = 0; l < 2; l++)
			{
				if (nodes[indexes[2 * l] + 3] != -1)
				{
					//Do nothing
				}
				else
				{
					nodes[indexes[2 * l] + 3] = eachNodes[QuadIndex[i]][mappedindexes[l]];
				}
			}
		}


	}
	else if (nTri == 0 and nQuadri == 6)
	{
		// std::cout<<"Brick\n";
		size = 8;
		for (int i = 0; i < 4; i++)
		{
			nodes[i] = eachNodes[0][i];
		}
		//We don't know if it is 0-3 1-2 or 0-1 2-3
		int indexes[4];
		int mappedindexes[2];
		int count = 0;
		for (int i = 1; i < 6; i++)//Each face
		{
			count = 0;
			for (int j = 0; j < 4; j++)//Each vert in this face
			{
				for (int k = 0; k < 4; k++)//Each first 4 
				{
					if (eachNodes[i][j] == nodes[k])
					{
						indexes[count++] = k;//Index inside first 4 nodes
						indexes[count++] = j;//Index in this face
					}
				}
			}
			if (count == 0)
			{
				continue;
			}
			int mini, maxi;
			mini = std::min(indexes[1], indexes[3]);
			maxi = std::max(indexes[1], indexes[3]);
			if (mini == 0 and maxi == 1)
			{
				mappedindexes[0] = 3 - indexes[1];
				mappedindexes[1] = 3 - indexes[3];
			}
			else if (mini == 2 and maxi == 3)
			{
				mappedindexes[0] = 3 - indexes[1];
				mappedindexes[1] = 3 - indexes[3];
			}
			else if (mini == 0 and maxi == 3)
			{
				mappedindexes[0] = indexes[1] == 0 ? 1 : 2;
				mappedindexes[1] = indexes[3] == 0 ? 1 : 2;
			}
			else if (mini == 1 and maxi == 2)
			{
				mappedindexes[0] = indexes[1] == 1 ? 0 : 3;
				mappedindexes[1] = indexes[3] == 1 ? 0 : 3;
			}
			else
			{
				std::cout << "Something I don't know happens\n";
			}
			for (int l = 0; l < 2; l++)
			{
				if (nodes[indexes[2 * l] + 4] != -1)
				{
					if (nodes[indexes[2 * l] + 4] != eachNodes[i][mappedindexes[l]])
					{

						std::cout << "Not cohesient to the existing one\n";
						std::cout << nodes[indexes[2 * l] + 4] << '\t' << eachNodes[i][mappedindexes[l]];
						std::cout << "nodes are :";
						for (int k = 0; k < 8; k++)
						{
							std::cout << nodes[k] << ',';
						}
						std::cout << '\n';
						std::cout << "starts are :";
						for (int k = 0; k < 4; k++)
						{
							std::cout << eachNodes[i][k] << ',';
						}
						std::cout << '\n';
						std::cin.get();
					}
				}
				else
				{
					nodes[indexes[2 * l] + 4] = eachNodes[i][mappedindexes[l]];
					if (checkRepeat())
					{
						std::cout << nodes[indexes[2 * l] + 4] << '\t' << eachNodes[i][mappedindexes[l]];
						std::cout << "nodes are :";
						for (int k = 0; k < 8; k++)
						{
							std::cout << nodes[k] << ',';
						}
						std::cout << '\n';
						for (int k = 0; k < 2; k++)
						{
							std::cout <<"mappedindexes" << mappedindexes[k] << '\n';
						}
						for (int k = 0; k < 4; k++)
						{
							std::cout << indexes[k] << '\t';
						}
						std::cout << '\n';

						std::cout << '\n';
						std::cout << "starts are :";
						for (int k = 0; k < 4; k++)
						{
							std::cout << eachNodes[i][k] << ',';
						}
						std::cout << '\n';
						std::cin.get();
					}
				}
			}
		}
	}
	else
	{
		std::cout << "Something wrong with this brick with nTri and nQuadri being " << nTri << " " << nQuadri << '\n';
	}
}