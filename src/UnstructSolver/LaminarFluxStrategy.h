#pragma once
#include "FluxStrategy.h"
#include "UnstructIntegrator.h"
#include <algorithm>
#include "math.h"
#include <fstream>
class LaminarFluxStrategy : public FluxStrategy
{
private:
    double fs_Temperature;
    double fs_Pr;
    double fs_Ma;
    double fs_Re;

    // double **Spectrum_cell_v = NULL;

    double ** mu = NULL;
    double ** mu_edge = NULL;
public:

    double **getMu()
    {
        return mu;
    }
    double **getMuEdge()
    {
        return mu_edge;
    }
    double getMuCell(int curMesh,int curCell)
    {
        return mu[curMesh][curCell];
    }
    double getMuEdge(int curMesh,int curEdge)
    {
        return mu_edge[curMesh][curEdge];

    }
    double getMuOverPrCell(int curMesh,int curCell)
    {
        return mu[curMesh][curCell]/fs_Pr;

    }
    double getMuOverPrEdge(int curMesh,int curEdge)
    {
        return mu_edge[curMesh][curEdge]/fs_Pr;

    }


    
    
public:
    LaminarFluxStrategy(UnstructTopologyHolder *hder, TopologyHolderStrategy *hder_strategy)
        : FluxStrategy(hder, hder_strategy)
    {
        mu = new double*[d_nmesh];
        mu_edge = new double*[d_nmesh];
        // Spectrum_cell_v = new double*[d_nmesh];
        for(int i = 0;i<d_nmesh;i++)
        {
            mu[i] = new double[d_hder->nCells(i)];
            mu_edge[i] = new double[d_hder->nEdges(i)];
            // Spectrum_cell_v[i] = new double[d_hder->nCells(i)];
        }
        fs_Temperature = d_hder_strategy->fs_Temperature;
        std::cout<<"LaminarFluxStrategy::fs_Temperature "<<fs_Temperature<<'\n';
        fs_Pr = d_hder_strategy->fs_Pr;
        std::cout<<"LaminarFluxStrategy::fs_Pr "<<fs_Pr<<'\n';
        fs_Ma = d_hder_strategy->fs_mach;
        std::cout<<"LaminarFluxStrategy::fs_Ma "<<fs_Ma<<'\n';
        fs_Re = d_hder_strategy->fs_Re;
        std::cout<<"LaminarFluxStrategy::fs_Re "<<fs_Re<<'\n';
    }


    void computeMu()
    {
        const double cc = 110.5/fs_Temperature;
        const double Gamma = 1.4;
        double T;
        double** U = d_hder_strategy->getU();
        double** U_edge = d_hder_strategy->getUEdge();
        for(int i = 0;i<d_nmesh;i++)
        {
            for(int k = 0;k<d_hder->nCells(i);k++)
            {
                T = U[i][1+d_NEQU*k]/U[i][0+d_NEQU*k];
                mu[i][k] = (1.0+cc)/(T+cc)*powf64(T,1.5)*sqrtf64(Gamma)*fs_Ma/fs_Re;
            }
            for(int k = 0;k<d_hder->nEdges(i);k++)
            {
                T = U_edge[i][1+d_NEQU*k]/U_edge[i][0+d_NEQU*k];
                mu_edge[i][k] = (1.0+cc)/(T+cc)*powf64(T,1.5)*sqrtf64(Gamma)*fs_Ma/fs_Re;
            }
        }
    }

    void computeSpectrum()
    {

    }


    void computeFlux()
    {
        std::cout<<"Do nothing\n";
        return;
    }

    void computeFlux(int istage)
    {
        const double Gamma = 1.4;
        //Get mu_edge from FluxFlow
        double** U = d_hder_strategy->getU();
        double** U_edge = d_hder_strategy->getUEdge();
        double** Residual = d_hder_strategy->getResidual();
        double Fv[5];
        GeomElements::vector3d<2,double>** gradU = (GeomElements::vector3d<2,double>**)d_hder_strategy->getGradientPrimitive();
        GeomElements::vector3d<2,double>* gradEdge = new GeomElements::vector3d<2,double>[d_dim+1];
        GeomElements::vector3d<2,double>** gradT = (GeomElements::vector3d<2,double>**)d_hder_strategy->getGradientT();
        for(int i = 0; i < d_nmesh; i++)
        {
            auto& curBlk = d_hder->blk2D[i];
            GeomElements::vector3d<2,double> average;
            GeomElements::vector3d<2,double> Left2Right;
            for(int k = 0; k < d_hder->nEdges(i); k++)
            {
                auto& curEdge = curBlk.d_localEdges[k];
                int lC = curEdge.lCInd();
                int rC = curEdge.rCInd();
                double mu_at_face = d_hder_strategy->getMuEdge(i,k);
                double efficient = 2.0/3.0*mu_at_face;
                //First compute Gradient at Edge.
                if(rC>=0)//Inner
                {
                    Left2Right = curBlk.d_localCells[rC].center() - curBlk.d_localCells[lC].center();
                    for(int j = 0; j<d_dim;j++)
                    {
                        average = (gradU[i][2+j+d_NEQU*lC] + gradU[i][2+j+d_NEQU*rC])*0.5;
                        gradEdge[j] = average + Left2Right*((U[i][2+j+d_NEQU*rC]-U[i][2+j+d_NEQU*lC] - average.dot_product(Left2Right))/Left2Right.L2Square());
                    }
                    average = (gradT[i][lC] + gradT[i][rC])*0.5;
                    double rightT = U[i][1+d_NEQU*rC]/U[i][0+d_NEQU*rC];
                    double leftT = U[i][1+d_NEQU*lC]/U[i][0+d_NEQU*lC];
                    gradEdge[d_dim] = average + Left2Right*((rightT - leftT - average.dot_product(Left2Right))/Left2Right.L2Square());
                }
                else
                {
                    Left2Right = curEdge.center() - curBlk.d_localCells[lC].center();
                    for(int j = 0;j <d_dim;j++)
                    {
                        gradEdge[j] = gradU[i][2+j+d_NEQU*lC];
                    }
                    if(rC==GeomElements::edge3d<2>::BoundaryType::WALL)
                    {
                        gradEdge[d_dim] = GeomElements::vector3d<2,double>(0,0);
                    }
                    else if(rC==GeomElements::edge3d<2>::BoundaryType::FARFIELD)
                    {
                        gradEdge[d_dim] = gradT[i][lC];
                    }
                }
                //With Edge Gradient, compute tao.
                double tao[2][2];
                if(d_dim==2)
                {
                    GeomElements::vector3d<2,double> tao[2];
                    GeomElements::vector3d<2,double> velocity_edge(U_edge[i][2+d_NEQU*k],U_edge[i][3+d_NEQU*k]);
                    
                    double taoxx = efficient*(2.0*gradEdge[0][0] - gradEdge[1][1]);
                    double taoyy = efficient*(2.0*gradEdge[1][1] - gradEdge[0][0]);
                    double taoxy = mu_at_face*(gradEdge[0][1] + gradEdge[1][0]);
                    
                    tao[0][0] = taoxx;
                    tao[1][1] = taoyy;
                    tao[0][1] = taoxy;
                    tao[1][0] = taoxy;

                    double qx = - Gamma/(Gamma -1) *(d_hder_strategy->getMuOverPrEdge(i,k)) * gradEdge[d_dim][0];
                    double qy = - Gamma/(Gamma -1) *(d_hder_strategy->getMuOverPrEdge(i,k)) * gradEdge[d_dim][1];
                    GeomElements::vector3d<2,double> normal = curEdge.normal_vector();
                    Fv[0] = 0.0;
                    Fv[1] = ((velocity_edge[0]*taoxx+velocity_edge[1]*taoxy-qx)*normal[0]
                    +(velocity_edge[0]*taoxy+velocity_edge[1]*taoyy-qy)*normal[1])*curEdge.area();
                    Fv[2] = (taoxx*normal[0] + taoxy*normal[1])*curEdge.area();
                    Fv[3] = (taoxy*normal[0] + taoyy*normal[1])*curEdge.area();
                    // Fv[1] = GeomElements::vector3d<2,double>(velocity_edge.dot_product(tao[0])-qx,
                    // velocity_edge.dot_product(tao[1])-qy).dot_product(curEdge.normal_vector())*curEdge.area();
                    //Fv[2] = tao[0].dot_product(curEdge.normal_vector())*curEdge.area();
                    // Fv[3] = tao[1].dot_product(curEdge.normal_vector())*curEdge.area();
                    
                    if (lC >= 0) // Left plus right minus
                    {
                        for (int j = 0; j < d_NEQU; j++)
                        {
                            Residual[i][j+d_NEQU*lC] -= Fv[j];
                        }
                    }
                    else
                    {
                        std::cout << "left cell ind can't be negative\n";
                        std::cin.get();
                    }
                    if (rC >= 0)
                    {
                        for (int j = 0; j < d_NEQU; j++)
                        {
                            Residual[i][j+d_NEQU*rC] += Fv[j];
                        }
                    }
                }
            }
        }
        // fout.close();
        delete[] gradEdge;
    }

};