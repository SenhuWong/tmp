#include<list>
#include<iostream>
#include"codetypes.h"
void test_list()
{
    std::list<INTERPLIST2*> intptr_list;
    for(int i=0;i<100;i++)
    {
        std::cout<<"i is "<<i<<" out of "<<100<<'\n';
        intptr_list.push_back(new INTERPLIST2[1]);
    }
}