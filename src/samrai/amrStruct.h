#include <cmath>
#include <iostream>
template <typename T, int porder, int dim>
class amrStruct
{
    T components[dim == 2 ? (porder + 1) * (porder + 1) : (porder + 1) * (porder + 1) * (porder + 1)];
    amrStruct()
    {
        
    }
    
    amrStruct(const amrStruct &another)
    {
        for (int i = 0; i < dim == 2 ? (porder + 1) * (porder + 1) : (porder + 1) * (porder + 1) * (porder + 1); i++)
        {
            components[i] = another.components[i];
        }
    }
    amrStruct &operator=(const amrStruct &another)
    {
        for (int i = 0; i < dim == 2 ? (porder + 1) * (porder + 1) : (porder + 1) * (porder + 1) * (porder + 1); i++)
        {
            components[i] = another.components[i];
        }
        return (*this);
    }
};

