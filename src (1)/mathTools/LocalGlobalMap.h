#pragma once

#include<set>
class LocalGlobalMap
{
public:
	std::pair<int, int>* maps;
	int size;
	LocalGlobalMap();
	LocalGlobalMap(LocalGlobalMap& another);
	//iter must be starting from 1.
	LocalGlobalMap(int Ssize, std::set<int, std::less<int>>::iterator begin, std::set<int, std::less<int>>::iterator end);
	~LocalGlobalMap()
	{
		delete[] maps;
	}
	//mapping at both directions should be starting from 1,so does the construction.
	int local2global(int to_find_local)//starting from 1
	{
		int result = maps[to_find_local - 1].first;
		return result;
	}
	int global2local(int to_find);
	std::pair<int, int>* begin()
	{
		return maps;
	}
	int sizeIs()
	{
		return size;
	}
};