#pragma once
#include "TopologyHolderStrategy.h"
#include "UnstructIntegrator.h"
#include "toolBox/edge3d_int.h"
#include "math.h"
class BoundaryStrategy
{
private:
    int d_nmesh = -1;
    int d_NEQU = -1;
    int d_dim = -1;
    // Flow params of freestream 
    double fsnd_density;
    double fsnd_pressure;   
    double fsnd_velocity_components[3];
    double fsnd_soundSpeed;
    double fs_Temperature;
    double fs_Ma;
    double fs_Re;
    double fs_Pr;
    double fs_r;
    double fs_mu;
    const double Cp = 1005.0;
    //
    UnstructTopologyHolder *d_hder = nullptr;
    TopologyHolderStrategy *d_hder_strategy = nullptr;
public:
    BoundaryStrategy(UnstructTopologyHolder *hder, TopologyHolderStrategy *hder_strategy)
        :d_nmesh(hder_strategy->getNMesh()),
        d_NEQU(hder_strategy->getNEquation()),
        d_dim(hder_strategy->getNDim()),
        d_hder(hder), d_hder_strategy(hder_strategy)
    {
        fsnd_density = d_hder_strategy->fsnd_density;
        std::cout<<"BoundaryStrategy::fsnd_density "<<fsnd_density<<'\n';
        fsnd_pressure = d_hder_strategy->fsnd_pressure;
        std::cout<<"BoundaryStrategy::fsnd_pressure "<<fsnd_pressure<<'\n';
        
        for(int i = 0;i<d_dim;i++)
        {
            fsnd_velocity_components[i] = d_hder_strategy->fsnd_velocity_components[i];
        }
        fsnd_soundSpeed = d_hder_strategy->fsnd_soundSpeed;
        std::cout<<"BoundaryStrategy::fsnd_soundSpeed "<<fsnd_soundSpeed<<'\n';
        fs_Temperature = d_hder_strategy->fs_Temperature;
        std::cout<<"BoundaryStrategy::fs_Temperature "<<fs_Temperature<<'\n';
        fs_Ma = d_hder_strategy->fs_mach;
        std::cout<<"BoundaryStrategy::fs_Ma "<<fs_Ma<<'\n';
        fs_Re = d_hder_strategy->fs_Re;
        std::cout<<"BoundaryStrategy::fs_Re "<<fs_Re<<'\n';
        fs_Pr = d_hder_strategy->fs_Pr;
        std::cout<<"BoundaryStrategy::fs_Pr "<<fs_Pr<<'\n';
        fs_r = powf64(fs_Pr,1.0/3.0);
        std::cout<<"BoundaryStrategy::fs_r "<<fs_r<<'\n';
        fs_mu = d_hder_strategy->fs_mu;
        std::cout<<"BoundaryStrategy::fs_mu "<<fs_mu<<'\n';
    }
    //Boundary at wall should be speciall treated when viscosity yields its effect.
    //And also mu should be calculated.
    void HandleBoundary()
    {
        if(d_hder_strategy->isInvicid())
        {
            HandleInvicidBoundary();
        }
        else
        {
            HandleViscousBoundary();
        }

    }

    void HandleViscousBoundary()
    {
        const double Gamma =1.4;
        double** U_edge = d_hder_strategy->getUEdge();
        double** U = d_hder_strategy->getU();

        // For far field boundary
        double cc = 110.4/fs_Temperature;
        //

        GeomElements::vector3d<2, double> Leftvelocity;
        GeomElements::vector3d<2, double> Rightvelocity;

        GeomElements::vector3d<2, double> normalVec;

        double normalVelocityComponentLeft;
        GeomElements::vector3d<2, double> tangenVelocityLeft;
        double tangenVelocityComponentLeftSquare;

        double normalVelocityComponentRight;
        GeomElements::vector3d<2, double> tangenVelocityRight;

        for (int i = 0; i < d_nmesh; i++)
        {
            auto &curBlk = d_hder->blk2D[i];
            for (int k = 0; k < d_hder->nEdges(i); k++)
            {
                auto &curEdge = curBlk.d_localEdges[k];
                int lC = curEdge.lCInd();
                int rC = curEdge.rCInd();

                if (rC >= 0) // Normal
                {
                    for (int j = 0; j < d_NEQU; j++)
                    {
                        U_edge[i][j+d_NEQU*k] = (U[i][j+d_NEQU*lC] + U[i][j+d_NEQU*rC]) / 2;
                    }
                }
                else if (rC == GeomElements::edge3d<2>::BoundaryType::WALL) // Wall
                {
                    U_edge[i][0+d_NEQU*k] = U[i][0+d_NEQU*lC];
                    U_edge[i][1+d_NEQU*k] = U[i][1+d_NEQU*lC];
                    for(int j = 0;j<d_dim;j++)
                    {
                        U_edge[i][2+j+d_NEQU*k] = 0;
                    }
                }
                else if (rC == GeomElements::edge3d<2>::BoundaryType::FARFIELD) // Far field
                {
                    // This is the derivation
                    for(int j = 0;j<d_dim;j++)
                    {
                        Leftvelocity[j] = U[i][2+j+d_NEQU*lC];
                    }
                    normalVec = curEdge.normal_vector();
                    normalVelocityComponentLeft = Leftvelocity.dot_product(normalVec);
                    tangenVelocityLeft = Leftvelocity - normalVec * normalVelocityComponentLeft;

                    Rightvelocity[0] = fsnd_velocity_components[0];
                    Rightvelocity[1] = fsnd_velocity_components[1];

                    normalVec = curEdge.normal_vector();
                    normalVelocityComponentRight = Rightvelocity.dot_product(normalVec);
                    tangenVelocityRight = (Rightvelocity - normalVec * normalVelocityComponentRight);

                    GeomElements::vector3d<2, double> velo = Leftvelocity;
                    double rho_left = U[i][0+d_NEQU*lC];
                    double p_left = U[i][1+d_NEQU*lC];

                    if (fabsf64(normalVelocityComponentRight) < fsnd_soundSpeed) // Sub sonic
                    {
                        double Vnb = (normalVelocityComponentLeft + normalVelocityComponentRight + 2 * (sqrt(Gamma * p_left / rho_left) - fsnd_soundSpeed) / (Gamma - 1)) / 2;
                        double Cb = (normalVelocityComponentLeft - normalVelocityComponentRight + 2 * (fsnd_soundSpeed + sqrt(Gamma * p_left / rho_left)) / (Gamma - 1)) * (Gamma - 1) / 4;
                        if (normalVelocityComponentRight >= 0) // Out
                        {
                        // double Vnb = (normalVelocityComponentLeft + normalVelocityComponentRight + 2 * (sqrt(Gamma * p_left / rho_left) - fsnd_soundSpeed) / (Gamma - 1)) / 2;
                        // double Cb = (normalVelocityComponentLeft - normalVelocityComponentRight + 2 * (sqrt(Gamma * p_left / rho_left) + fsnd_soundSpeed) / (Gamma - 1)) * (Gamma - 1) / 4;
                        //   s should equal to Soo
                            double Sb = (p_left / pow(rho_left, Gamma));
                            GeomElements::vector3d<2, double> EdgeVelocity = tangenVelocityLeft + normalVec * Vnb;
                            U_edge[i][0+d_NEQU*k] = pow((Cb * Cb / (Gamma * Sb)), (1.0 / (Gamma - 1)));
                            U_edge[i][1+d_NEQU*k] = U_edge[i][0+d_NEQU*k] * Cb * Cb / Gamma;
                            U_edge[i][2+d_NEQU*k] = EdgeVelocity[0];
                            U_edge[i][3+d_NEQU*k] = EdgeVelocity[1];
                        }
                        else
                        {
                            // double Vnb = (normalVelocityComponentLeft + normalVelocityComponentRight + 2 * (fsnd_soundSpeed - sqrt(Gamma * p_left / rho_left)) / (Gamma - 1)) / 2;
                        // double Cb = (normalVelocityComponentRight - normalVelocityComponentLeft + 2 * (sqrt(Gamma * p_left / rho_left) + fsnd_soundSpeed) / (Gamma - 1)) * (Gamma - 1) / 4;
                            double Sb = (fsnd_pressure / pow(fsnd_density, Gamma));
                            GeomElements::vector3d<2, double> Edgevelocity = tangenVelocityRight + normalVec * Vnb;
                            U_edge[i][0+d_NEQU*k] = pow((Cb * Cb / (Gamma * Sb)), (1.0 / (Gamma - 1)));
                            U_edge[i][1+d_NEQU*k] = U_edge[i][0+d_NEQU*k] * Cb * Cb / Gamma;
                            U_edge[i][2+d_NEQU*k] = Edgevelocity[0];
                            U_edge[i][3+d_NEQU*k] = Edgevelocity[1];
                        }
                    }
                    else
                    {
                        if (normalVelocityComponentRight >= 0)
                        {
                            U_edge[i][0+d_NEQU*k] = U[i][0+d_NEQU*lC];
                            U_edge[i][1+d_NEQU*k] = U[i][1+d_NEQU*lC];
                            U_edge[i][2+d_NEQU*k] = U[i][2+d_NEQU*lC];
                            U_edge[i][3+d_NEQU*k] = U[i][3+d_NEQU*lC];
                        }
                        else
                        {
                            U_edge[i][0+d_NEQU*k] = fsnd_density;
                            U_edge[i][1+d_NEQU*k] = fsnd_pressure;
                            U_edge[i][2+d_NEQU*k] = fsnd_velocity_components[0];
                            U_edge[i][3+d_NEQU*k] = fsnd_velocity_components[1];
                        }
                    }
                }
                else
                {
                    throw std::runtime_error("Undefined Boundary Type\n");
                }
            }
        }

    }
    
    void HandleInvicidBoundary()
    {
        const double Gamma =1.4;
        double** U_edge = d_hder_strategy->getUEdge();
        double** U = d_hder_strategy->getU();
        //Establish a method to pass in freestream variables.
        double fs_primVars[6];
        d_hder_strategy->getFreeStreamVars(fs_primVars);
        double fsnd_density = fs_primVars[0];
        double fsnd_pressure = fs_primVars[1];
        
        double fsnd_velocity_components[3];
        for(int i = 2;i<d_NEQU;i++)
        {
            fsnd_velocity_components[i-2] = fs_primVars[i];
        }
        double fsnd_soundSpeed = fs_primVars[d_NEQU];
        
        
        GeomElements::vector3d<2, double> Leftvelocity;
        GeomElements::vector3d<2, double> Rightvelocity;

        GeomElements::vector3d<2, double> normalVec;

        double normalVelocityComponentLeft;
        GeomElements::vector3d<2, double> tangenVelocityLeft;

        double normalVelocityComponentRight;
        GeomElements::vector3d<2, double> tangenVelocityRight;

        for (int i = 0; i < d_nmesh; i++)
        {
            auto &curBlk = d_hder->blk2D[i];
            for (int k = 0; k < curBlk.nEdges(); k++)
            {
                auto &curEdge = curBlk.d_localEdges[k];
                int lC = curEdge.lCInd();
                int rC = curEdge.rCInd();

                if (rC >= 0) // Normal
                {
                    for (int j = 0; j < d_NEQU; j++)
                    {
                        U_edge[i][j+d_NEQU*k] = (U[i][j+d_NEQU*lC] + U[i][j+d_NEQU*rC]) / 2;
                    }
                }
                else if (rC == GeomElements::edge3d<2>::BoundaryType::WALL) // Wall
                {
                    Leftvelocity[0] = U[i][2+d_NEQU*lC];
                    Leftvelocity[1] = U[i][3+d_NEQU*lC];
                    normalVec = curEdge.normal_vector();
                    normalVelocityComponentLeft = Leftvelocity.dot_product(normalVec);
                    tangenVelocityLeft = Leftvelocity - normalVec * normalVelocityComponentLeft;
                    if (tangenVelocityLeft.dot_product(normalVec) > 1.0e-7)
                    {
                        std::cout<<tangenVelocityLeft[0]<<","<<tangenVelocityLeft[1]<<'\n';
                        std::cout<<normalVec[0] <<","<<normalVec[1]<<'\n';
                        std::cout<<tangenVelocityLeft.dot_product(normalVec)<<'\n';
                        throw std::runtime_error("Something wrong with Vector computation\n");
                    }
                    U_edge[i][0+d_NEQU*k] = U[i][0+d_NEQU*lC];
                    U_edge[i][1+d_NEQU*k] = U[i][1+d_NEQU*lC];
                    U_edge[i][2+d_NEQU*k] = tangenVelocityLeft[0];
                    U_edge[i][3+d_NEQU*k] = tangenVelocityLeft[1];
                }
                else if (rC == GeomElements::edge3d<2>::BoundaryType::FARFIELD) // Far field
                {
                    // This is the derivation
                    Leftvelocity[0] = U[i][2+d_NEQU*lC];
                    Leftvelocity[1] = U[i][3+d_NEQU*lC];

                    normalVec = curEdge.normal_vector();
                    normalVelocityComponentLeft = Leftvelocity.dot_product(normalVec);
                    tangenVelocityLeft = Leftvelocity - normalVec * normalVelocityComponentLeft;

                    Rightvelocity[0] = fsnd_velocity_components[0];
                    Rightvelocity[1] = fsnd_velocity_components[1];

                    normalVec = curEdge.normal_vector();
                    normalVelocityComponentRight = Rightvelocity.dot_product(normalVec);
                    tangenVelocityRight = (Rightvelocity - normalVec * normalVelocityComponentRight);

                    GeomElements::vector3d<2, double> velo = Leftvelocity;
                    double rho_left = U[i][0+d_NEQU*lC];
                    double p_left = U[i][1+d_NEQU*lC];

                    if (fabsf64(normalVelocityComponentRight) < fsnd_soundSpeed) // Sub sonic
                    {
                        double Vnb = (normalVelocityComponentLeft + normalVelocityComponentRight + 2 * (sqrt(Gamma * p_left / rho_left) - fsnd_soundSpeed) / (Gamma - 1)) / 2;
                        double Cb = (normalVelocityComponentLeft - normalVelocityComponentRight + 2 * (fsnd_soundSpeed + sqrt(Gamma * p_left / rho_left)) / (Gamma - 1)) * (Gamma - 1) / 4;
                        if (normalVelocityComponentRight >= 0) // Out
                        {
                        // double Vnb = (normalVelocityComponentLeft + normalVelocityComponentRight + 2 * (sqrt(Gamma * p_left / rho_left) - fsnd_soundSpeed) / (Gamma - 1)) / 2;
                        // double Cb = (normalVelocityComponentLeft - normalVelocityComponentRight + 2 * (sqrt(Gamma * p_left / rho_left) + fsnd_soundSpeed) / (Gamma - 1)) * (Gamma - 1) / 4;
                        //   s should equal to Soo
                            double Sb = (p_left / pow(rho_left, Gamma));
                            GeomElements::vector3d<2, double> EdgeVelocity = tangenVelocityLeft + normalVec * Vnb;
                            U_edge[i][0+d_NEQU*k] = pow((Cb * Cb / (Gamma * Sb)), (1.0 / (Gamma - 1)));
                            U_edge[i][1+d_NEQU*k] = U_edge[i][0+d_NEQU*k] * Cb * Cb / Gamma;
                            U_edge[i][2+d_NEQU*k] = EdgeVelocity[0];
                            U_edge[i][3+d_NEQU*k] = EdgeVelocity[1];
                        }
                        else
                        {
                            // double Vnb = (normalVelocityComponentLeft + normalVelocityComponentRight + 2 * (fsnd_soundSpeed - sqrt(Gamma * p_left / rho_left)) / (Gamma - 1)) / 2;
                        // double Cb = (normalVelocityComponentRight - normalVelocityComponentLeft + 2 * (sqrt(Gamma * p_left / rho_left) + fsnd_soundSpeed) / (Gamma - 1)) * (Gamma - 1) / 4;
                            double Sb = (fsnd_pressure / pow(fsnd_density, Gamma));
                            GeomElements::vector3d<2, double> Edgevelocity = tangenVelocityRight + normalVec * Vnb;
                            U_edge[i][0+d_NEQU*k] = pow((Cb * Cb / (Gamma * Sb)), (1.0 / (Gamma - 1)));
                            U_edge[i][1+d_NEQU*k] = U_edge[i][0+d_NEQU*k] * Cb * Cb / Gamma;
                            U_edge[i][2+d_NEQU*k] = Edgevelocity[0];
                            U_edge[i][3+d_NEQU*k] = Edgevelocity[1];
                        }
                    }
                    else
                    {
                        if (normalVelocityComponentRight >= 0)
                        {
                            U_edge[i][0+d_NEQU*k] = U[i][0+d_NEQU*lC];
                            U_edge[i][1+d_NEQU*k] = U[i][1+d_NEQU*lC];
                            U_edge[i][2+d_NEQU*k] = U[i][2+d_NEQU*lC];
                            U_edge[i][3+d_NEQU*k] = U[i][3+d_NEQU*lC];
                        }
                        else
                        {
                            U_edge[i][0+d_NEQU*k] = fsnd_density;
                            U_edge[i][1+d_NEQU*k] = fsnd_pressure;
                            U_edge[i][2+d_NEQU*k] = fsnd_velocity_components[0];
                            U_edge[i][3+d_NEQU*k] = fsnd_velocity_components[1];
                        }
                    }
                }
                else
                {
                    throw std::runtime_error("Undefined Boundary Type\n");
                }
            }
        }
    } 
};