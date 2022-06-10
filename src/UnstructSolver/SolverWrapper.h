#pragma once
#include "Euler2D.h"
#include "RungeKuttaStrategy.h"
#include "LUSGS_Integrator.h"
#include "ViscousFlow2D.h"
#include <mpi/mpi.h>
class WrapSolver
{
public:
    
    UnstructTopologyHolder* d_topology_holder = NULL;
    TopologyHolderStrategy* d_flow_strategy = NULL;
    TimeIntegrator* d_advance_strategy = NULL;

    //These 2 should be eliminated and control is given to Orchestra.
    int d_maxStep;
    double d_minRes;
    

    
public:
    void setTopology(UnstructFeeder* fder)
    {
        d_topology_holder = new UnstructTopologyHolder(fder);
        d_topology_holder->arrange_communication();
    }
    void setParams(const std::string&flowT,const std::string&advanceT,const std::string& solveT, 
                    int MaxStep,double MinRes,double ma, double re,double* aoa)
    {
        if(flowT=="Euler2D")
        {
            d_flow_strategy = new Euler2D();
        }
        else if(flowT =="ViscousFlow2D")
        {
            d_flow_strategy = new ViscousFlow2D();
        }
        else
        {
            throw std::runtime_error("Undefined Flow Type\n");
        }

        if(advanceT=="RungeKutta5")
        {
            d_advance_strategy = new RungeKuttaStrategy(5);
        }
        else if(advanceT=="LUSGS")
        {
            d_advance_strategy = new LUSGSIntegrator();
        }
        else
        {
            throw std::runtime_error("Undefined Advance Type\n");
        }        
        d_flow_strategy->set_fs_variable2(ma,*aoa,1.225,101325,re);
    }
    
    void init(UnstructFeeder* fder)
    {
        d_topology_holder->initialize_integrator(fder);
        d_topology_holder->arrange_communication();

        d_flow_strategy->initializeStrategy(d_topology_holder);

        d_advance_strategy->initializeStrategy(d_topology_holder,d_flow_strategy);

        d_advance_strategy->initializeData();
    }


    void Loop()
    {
        for(int i = 0;i<d_maxStep;i++)
        {
            d_advance_strategy->singleStep(i);
            if(i%5000==0)
            {
                std::cout<<i<<'\n';
                int cur_proc;
                d_flow_strategy->AllCellCommunication(d_flow_strategy->getU());
                MPI_Comm_rank(MPI_COMM_WORLD,&cur_proc);
                for(int j = 0;j<d_topology_holder->d_nmesh;j++)
                {
                    d_flow_strategy->writeCellData("CellDataParallelLUSGS"+std::to_string(i+1),cur_proc,j,d_flow_strategy->getU());
                    std::string Cp_name = "OneAndOnlyLegendaryCp_"+std::to_string(j+1);
                    d_flow_strategy->outPutCp(Cp_name,j,i);
                }
            }
        }
    }




};
