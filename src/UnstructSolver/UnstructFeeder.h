#pragma once
#include "UnstructBlock.h"
#include "../reader/Reader.h"
class UnstructFeeder: public virtual Reader
{
public:
    UnstructBlock2D* UBs2 = NULL;
    UnstructBlock3D* UBs3 = NULL;
    UnstructFeeder()
    {

    }
    ~UnstructFeeder()
    {

    }

    void initialize()
    {
        if(d_dim==2)
        {
            if(UBs2) 
                delete[] UBs2;
            
            UBs2 = new UnstructBlock2D[d_nmesh];
        }
        else if(d_dim==3)
        {
            if(UBs3)
                delete[] UBs3;
            UBs3 = new UnstructBlock3D[d_nmesh];
        }
    }
    
    virtual void readForUnstruct() = 0;
};