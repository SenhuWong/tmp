#pragma once
#include<time.h>
#include<chrono>
#include<iostream>
class simpleTimer
{
    std::chrono::time_point<std::chrono::steady_clock> starttime,endtime;
public:
    simpleTimer()
    {
        starttime = std::chrono::steady_clock::now();
    }

    ~simpleTimer()
    {
        endtime = std::chrono::steady_clock::now();

        auto duration = endtime - starttime;
        std::cout<<"Time cost:"<< duration.count()/1000000000<<" s"<<'\n';
    }
};