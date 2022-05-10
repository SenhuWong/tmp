#include "LocalGlobalMap.h"

LocalGlobalMap::LocalGlobalMap() {}
LocalGlobalMap::LocalGlobalMap(LocalGlobalMap& another)
{
	size = another.size;
	maps = new std::pair<int, int>[size];
	for (int i = 0; i < size; i++)
	{
		maps[i] = another.maps[i];
	}
}
//iter must be starting from 1.
LocalGlobalMap::LocalGlobalMap(int Ssize, std::set<int, std::less<int>>::iterator begin, std::set<int, std::less<int>>::iterator end)
{
	size = Ssize;
	int count = 0;
	maps = new std::pair<int, int>[size];
	for (auto iter = begin; iter != end; iter++)
	{
		//std::cout<<*iter<<'\t'<<count+1<<'\n';
		maps[count] = std::pair<int, int>(*iter, count + 1);
		++count;
	}
}

int LocalGlobalMap::global2local(int to_find)
{
	int mid;
	int min = 0;
	int max = size - 1;

	bool signal = true;
	if (maps[min].first > to_find or maps[max].first < to_find)
	{
		return -1;
	}
	else if (maps[min].first == to_find)
	{
		return maps[min].second;
	}
	else if (maps[max].first == to_find)
	{
		return maps[max].second;
	}
	while (signal)
	{
		//std::cout<<"max"<<max<<"min"<<min<<'\n';
		mid = static_cast<int>((min + max) / 2);
		if (maps[mid].first > to_find)
		{
			max = mid;
		}
		else if (maps[mid].first < to_find)
		{
			min = mid;
		}
		else
		{
			return maps[mid].second;
		}
		if (max - min <= 1)
		{
			//std::cout<<max<<'\t'<<min<<'\n';
			return -1;
		}
	}
	return -1;
}