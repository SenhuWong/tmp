//2 Basic interface to provide tioga with grid and boundary information
//#include "tioga.h"
#pragma once
#include <vector>
#include "../reader/Reader.h"
#include "BlockInfo.h"
class TiogaFeeder: public virtual Reader
{
public:
    std::vector<BlockInfo> bis; //in SAMRAI FEEDER
public:
    TiogaFeeder();
    ~TiogaFeeder();
    void initialize()
    {
        bis.clear();
        bis.resize(d_nmesh);
    }
    virtual void readForTioga() = 0;
    //resize and filling in for bis
    //Tioga should take this in.
    //virtual void registerTioga(TIOGA::tioga* tg);
    //virtual void registerTiogaBoundary(TIOGA::tioga* tg);
};