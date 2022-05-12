#include "TopologyHolderStrategy.h"
#include "UnstructIntegrator.h"
#include "toolBox/edge3d_int.h"
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
        fsnd_pressure = d_hder_strategy->fsnd_pressure;
        for(int i = 0;i<d_dim;i++)
        {
            fsnd_velocity_components[i] = d_hder_strategy->fsnd_velocity_components[i];
        }
        fsnd_soundSpeed = d_hder_strategy->fsnd_soundSpeed;
        fs_Temperature = d_hder_strategy->fs_Temperature;
        fs_Ma = d_hder_strategy->fs_mach;
        fs_Re = d_hder_strategy->fs_Reynolds_number;
        fs_Pr = d_hder_strategy->fs_Prandtl;

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

        double** mu = d_hder_strategy->getMu();
        double** mu_edge = d_hder_strategy->getMuEdge();
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
                    mu_edge[i][k] = (mu[i][lC] + mu[i][rC]) / 2;
                }
                else if (rC == GeomElements::edge3d<2>::BoundaryType::WALL) // Wall
                {
                    for(int j = 0;j<d_dim;j++)
                    {
                        Leftvelocity[i] = U[i][2+j+d_NEQU*lC];
                    }
                    normalVec = curEdge.normal_vector();
                    normalVelocityComponentLeft = Leftvelocity.dot_product(normalVec);
                    tangenVelocityLeft = Leftvelocity - normalVec * normalVelocityComponentLeft;
                    tangenVelocityComponentLeftSquare = tangenVelocityLeft.L2Square();
                    if (tangenVelocityLeft.dot_product(normalVec) > 1.0e-11)
                    {
                        throw std::runtime_error("Something wrong with Vector computation\n");
                    }
                    double T = U_edge[i][1+d_NEQU*k]/U_edge[i][0+d_NEQU*k];
                    double Tw = T +0.5* fs_r *tangenVelocityComponentLeftSquare/Cp;
                    double Pw = U_edge[i][1+d_NEQU*k];
                    double Dw = Pw/Tw;
                    double muw = (1.0+cc)/(Tw+cc)*powf64(Tw,1.5)*sqrtf64(Gamma)*fs_Ma/fs_Re;

                    U_edge[i][0+d_NEQU*k] = Dw;
                    U_edge[i][1+d_NEQU*k] = Pw;
                    for(int j = 2;j<d_dim+2;j++)
                    {
                        U_edge[i][j+d_NEQU*k] = 0;
                    }
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

                    if (abs(normalVelocityComponentRight) < fsnd_soundSpeed) // Sub sonic
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
                            double T = U_edge[i][1+d_NEQU*k]/U_edge[i][0+d_NEQU*k];
                            mu_edge[i][k] = (1.0+cc)/(T+cc)*powf64(T,1.5)*sqrtf64(Gamma)*fs_Ma/fs_Re;
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
                            double T = U_edge[i][1+d_NEQU*k]/U_edge[i][0+d_NEQU*k];
                            mu_edge[i][k] = (1.0+cc)/(T+cc)*powf64(T,1.5)*sqrtf64(Gamma)*fs_Ma/fs_Re;
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
                            mu_edge[i][k] = mu[i][lC];
                        }
                        else
                        {
                            U_edge[i][0+d_NEQU*k] = fsnd_density;
                            U_edge[i][1+d_NEQU*k] = fsnd_pressure;
                            U_edge[i][2+d_NEQU*k] = fsnd_velocity_components[0];
                            U_edge[i][3+d_NEQU*k] = fsnd_velocity_components[1];
                            mu_edge[i][k] = fsnd_pressure/fsnd_density;
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
                    if (tangenVelocityLeft.dot_product(normalVec) > 1.0e-11)
                    {
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

                    if (abs(normalVelocityComponentRight) < fsnd_soundSpeed) // Sub sonic
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