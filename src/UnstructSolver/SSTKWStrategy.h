#include "TopologyHolderStrategy.h"
#include "UnstructIntegrator.h"
#include "toolBox/edge3d_int.h"
#include "toolBox/vector3d.h"
#include <math.h>
template<int ndim>
class SSTkomegaModel
{
private:
    //Constants used.
    const double C1 = 10.0;
    const double C2 = 3.0;
    const double L = 1.0;

    const double beta1 = 0.075;
    const double beta2 = 0.0828;
    const double sigma_k1 = 0.85;
    const double sigma_k2 = 1.0;
    const double sigma_w1 = 0.5;
    const double sigma_w2 = 0.856;
    const double Cw1 = 0.533;
    const double Cw2 = 0.440;

    const double a1 = 0.31;
    const double _beta = 0.09;
    const double Cepsilion = 0.61;
    const double Comega = 0.78;

    //FreeStream Variables.
    double fs_mut = 0.0;
    double fs_pressure;
    double fs_density;
    double fs_Ma;
    double Gamma = 1.4;
    double fs_Re;
    double fs_VelocityMagnitude;

    double fs_turbulanceValue[2];
    double fs_consturbulanceValue[2];

    //From TopologyStrategy
    int d_nmesh;
    int d_NEQU = ndim + 2;
    UnstructTopologyHolder* d_hder;
    TopologyHolderStrategy* d_hder_strategy;


    double** turbulenceCellValue = NULL;
    double** turbulenceConsCellValue = NULL;
    double** cellF1 = NULL;
    double** cellF2 = NULL;

    GeomElements::vector3d<ndim,double>** turbulenceCellGrad = NULL;

    double** vCellMut = NULL;
    double** vFaceMut = NULL;

    double** vTurbulanceFaceValue = NULL;
    double** vFaceConvectiveTFlux = NULL;
    double** vFaceViscousTFlux = NULL;
    double** vCellTotalTFlux = NULL;


    //Get from Laminar Strategy.
    double** vFaceMul = NULL;
    //LUSGS should be done at LUSGS_Strategy.
    // double** vTurbulanceDiagOperator = NULL;
    // double** vTurbulanceCellDTW = NULL;


    
    //Get from TIOGA;
    double** vCellSolidDistance = NULL;


public:
    void ComputeFaceTValue();
    void ComputeCellTGradient();

    void ComputeSSTf1f2();

    void ComputeTurbulenceViscousCoefficient();

    void SolveConvectiveTFlux();

    void SolveViscousTFlux();

    double rans_P_Solver(int meshIndex,int cellIndex);

    void SolveTotalTFlux();

    SSTkomegaModel();

    void InitializeTurbulanceVariable();

    void SolveTurbulenceEquation();

    void LUSGSAdvance();

    //The last 3 looks like an AMR routine to me,I will leave it here.

};


template<int ndim>
SSTkomegaModel<ndim>::SSTkomegaModel()
{
    turbulenceCellValue = new double*[d_nmesh];
    turbulenceConsCellValue = new double*[d_nmesh];
    cellF1 = new double*[d_nmesh];
    cellF2 = new double*[d_nmesh];
    turbulenceCellGrad = new GeomElements::vector3d<ndim,double>*[d_nmesh];

    vCellMut = new double*[d_nmesh];
    vFaceMut = new double*[d_nmesh];

    vTurbulanceFaceValue = new double*[d_nmesh];
    vFaceConvectiveTFlux = new double*[d_nmesh];
    vFaceViscousTFlux = new double*[d_nmesh];
    vCellTotalTFlux = new double*[d_nmesh];

    vTurbulanceDiagOperator = new double*{d_nmesh};
    vTurbulanceCellDTW = new double*[d_nmesh];
    for(int i = 0;i<d_nmesh;i++)
    {
        turbulenceCellValue[i] = new double[d_NEQU*d_hder->nCells(i)];
        turbulenceConsCellValue[i] = new double[d_NEQU*d_hder->nCells(i)];
        cellF1[i] = new double[d_hder->nCells(i)];
        cellF2[i] = new double[d_hder->nCells(i)];
        turbulenceCellGrad[i] = new GeomElements::vector3d<ndim,double>[d_NEQU*d_hder->nCells(i)];

        vFaceConvectiveTFlux[i] = new double[d_NEQU*d_hder->nCells(i)];
        vFaceViscousTFlux[i] = new double[d_NEQU*d_hder->nCells(i)];
        vTurbulanceCellDTW[i] = new double[d_NEQU*d_hder->nCells(i)];
        vTurbulanceDiagOperator[i] = new double[d_NEQU*d_hder->nCells(i)];
    }
}

template<int ndim>
void SSTkomegaModel<ndim>::InitializeTurbulanceVariable()
{
    double T = fs_pressure/fs_density;
    double InfiniteMul = ((1.0+ 110.5/288)/(T+110.5/288))*sqrtf64(T*T*T)*fs_Ma*sqrtf64(Gamma)/fs_Re;

    fs_mut = InfiniteMul / 1000.0;

    fs_turbulanceValue[1] = C1 * fs_VelocityMagnitude /L;
    fs_turbulanceValue[0] = fs_mut * fs_turbulanceValue[1] / fs_density;

    fs_consturbulanceValue[1] = fs_turbulanceValue[1]*fs_density;
    fs_consturbulanceValue[0] = fs_turbulanceValue[0]*fs_density;
    for(int i = 0;i<d_nmesh;i++)
    {
        for(int k = 0;k<d_hder->nCells(i);k++)
        {
            int offset= 2*k;
            turbulenceCellValue[i][offset] = fs_turbulanceValue[0];
            turbulenceCellValue[i][offset+1] = fs_turbulanceValue[1];
            turbulenceConsCellValue[i][offset] = fs_consturbulanceValue[0];
            turbulenceConsCellValue[i][offset+1] = fs_consturbulanceValue[1];
        }
    }
}

template<int ndim>
void SSTkomegaModel<ndim>::ComputeFaceTValue()
{
    double U_edge = d_hder_strategy->getUEdge();
    for(int i = 0;i<d_nmesh;i++)
    {
        auto& curBlk = d_hder->blk2D[i];
        for(int k =0;k<d_hder->nEdges(i);k++)
        {
            auto& curEdge = curBlk.d_localEdges[k];
            int lC = curEdge.lCInd();
            int rC = curEdge.rCInd();
            if(rC>=0)//Inner
            {
                vTurbulanceFaceValue[2*k] = 0.5*(turbulenceCellValue[2*lC] + turbulenceCellValue[2*rC]);
                vTurbulanceFaceValue[2*k+1] = 0.5*(turbulenceCellValue[2*lC+1] + turbulenceCellValue[2*rC+1]);
            }
            else if(rC == GeomElements::edge3d<2>::BoundaryType::WALL)
            {
                vTurbulanceFaceValue[2*k] = 0.0;
                double dis = vCellSolidDistance[i][lC];
                vTurbulanceFaceValue[2*k+1] = 10.0 * 6.0 * vFaceMul[i][k] / (U_edge[i][d_NEQU*k]*beta1*dis*dis);
            }
            else if(rC == GeomElements::edge3d<2>::BoundaryType::FARFIELD)
            {
                GeomElements::vector3d<ndim,double> VelocityEdge(&U_edge[i][d_NEQU*k+2]);
                GeomElements::vector3d<ndim,double> norm = curEdge.normal_vector();
                double Vn = VelocityEdge.dot_product(norm);
                if(Vn > 0.0)
                {
                    vTurbulanceFaceValue[2*k] = turbulenceCellValue[2*lC];
                    vTurbulanceFaceValue[2*k+1] = turbulenceCellValue[2*lC+1];
                }
                else
                {
                    vTurbulanceFaceValue[2*k] = fs_turbulanceValue[0];
                    vTurbulanceFaceValue[2*k+1] = fs_turbulanceValue[1];
                }
            }
        }
    }
}

template<int ndim>
void SSTkomegaModel<ndim>::ComputeCellTGradient()
{
    for(int i = 0;i<d_nmesh;i++)
    {
        auto& curBlk = d_hder->blk2D[i];
        for(int k = 0;k<d_hder->nCells(i);k++)
        {
            turbulenceCellGrad[i][2*k].reset();
            turbulenceCellGrad[i][2*k+1].reset();
        }
        for(int k = 0;k<d_hder->nEdges(i);k++)
        {
            auto& curEdge = curBlk.d_localEdges[k];
            int lC = curEdge.lCInd();
            int rC = curEdge.rCInd();
            turbulenceCellGrad[i][2*lC] = turbulenceCellGrad[i][2*lC] + curEdge.normal_vector()*vTurbulanceFaceValue[2*k];
            turbulenceCellGrad[i][2*lC+1] = turbulenceCellGrad[i][2*lC+1] + curEdge.normal_vector()*vTurbulanceFaceValue[2*k+1];
            if(rC>=0)
            {
                turbulenceCellGrad[i][2*rC] = turbulenceCellGrad[i][2*rC] - curEdge.normal_vector()*vTurbulanceFaceValue[2*k];
                turbulenceCellGrad[i][2*rC+1] = turbulenceCellGrad[i][2*rC+1] - curEdge.normal_vector()*vTurbulanceFaceValue[2*k+1];
            }
        }
        for(int k = 0;k<d_hder->nCells(i);k++)
        {
            auto& curCell = curBlk.d_localCells[k];
            turbulenceCellGrad[i][2*k] = turbulenceCellGrad[i][2*k] / curCell.volume();
            turbulenceCellGrad[i][2*k+1] = turbulenceCellGrad[i][2*k+1] / curCell.volume(); 
        }
    }


}

template<int ndim>
void SSTkomegaModel<ndim>::ComputeSSTf1f2()
{
    double** U = d_hder_strategy->getU();
    double** mu = d_hder_strategy->getMu();
    double CD_kw, arg1, arg2, f1,f2;
    for(int i = 0;i<d_nmesh;i++)
    {
        for(int k = 0;k<d_hder->nCells(i);k++)
        {
            double gradVal = turbulenceCellGrad[i][2*k].dot_product(turbulenceCellGrad[i][2*k+1]);
            double rho = U[i][d_NEQU*k];
            double dis = vCellSolidDistance[i][k];
            double K = turbulenceCellValue[i][2*k];
            double w = turbulenceCellValue[i][2*k+1];
            double muL = mu[i][k];
            CD_kw = std::max<double>(2.0*rho*sigma_w2 / w * gradVal,1.0e-20);
            arg1 = std::max<double>(sqrtf64(K)/(0.09*w*dis),500.0*muL/(rho*w*dis*dis));
            arg1 = std::min<double>(arg1,4.0*rho*sigma_w2*K/(CD_kw*dis*dis));
            arg2 = std::max(2.0*sqrtf64(K)/(0.09*w*dis),500.0*muL/(rho*w*dis*dis));
            f1 = tanhf64(arg1*arg1*arg1*arg1);
            f2 = tanhf64(arg2*arg2);
            cellF1[i][k] = f1;
            cellF2[i][k] = f2;
        }
    }

}

template<int ndim>
void SSTkomegaModel<ndim>::ComputeTurbulenceViscousCoefficient()
{
    //First compute curlV
    double curl;
    GeomElements::vector3d<ndim,T>** gradU = (GeomElements::vector3d<ndim,T>**)d_hder_strategy->getGradientPrimitive();
    double** U_edge = d_hder_strategy->getUEdge();
    if(ndim==2)
    {
        for(int i = 0;i<d_nmesh;i++)
        {
            for(int k = 0;k<d_hder->nCells(i);k++)
            {
                curl = gradU[i][3+d_NEQU*k][0] - gradU[i][2+d_NEQU*k][1];
                vCellMut[i] = a1 * turbulenceConsCellValue[i][2*k] / std::max<double>(a1 * turbulenceCellValue[i][2*k+1], cellF2[i][k]* curl);
            }
        }
    }
    else if(ndim==3)
    {
        throw std::runtime_error("Not implemented yet\n");
    }
                
    for(int i = 0;i<d_nmesh;i++)
    {
        auto& curBlk = d_hder->blk2D[i];
        for(int k = 0;k<d_hder->nEdges(i);k++)
        {
            auto& curEdge =curBlk.d_localEdges[k];
            int lC = curEdge.lCInd();
            int rC = curEdge.rCInd();
            if(rC>=0)
            {
                vFaceMut[i][k] = 0.5*(vCellMut[lC] + vCellMut[rC])

            }
            else if(rC==GeomElements::edge3d<2>::BoundaryType::WALL)
            {
                vFaceMut[i][k] = 0.0;

            }
            else if(rC==GeomElements::edge3d<2>::BoundaryType::FARFIELD)
            {
                GeomElements::vector3d<ndim,double> velo(&(U_edge[i][2+d_NEQU*k]));
                GeomElements::vector3d<ndim,double> norm = curEdge.normal_vector();
                double Vn = velo.dot_product(norm);
                vFaceMut[i][k]  = Vn > 0.0 ? vCellMut[i][lC] : fs_mut;
            }
            else
            {
                throw std::runtime_error("Undefined boundary type\n");
            }

        }
    }
    
}

template<int ndim>
void SSTkomegaModel<ndim>::SolveConvectiveTFlux()
{
    double** U =d_hder_strategy->getU();
    double** U_edge = d_hder_strategy->getUEdge();
    for(int i = 0;i<d_nmesh;i++)
    {
        auto& curBlk = d_hder->blk2D[i];
        for(int k = 0;k<d_hder->nEdges(i);k++)
        {
            auto& curEdge = curBlk.d_localEdges[k];
            int lC = curEdge.lCInd();
            int rC = curEdge.rCInd();
            GeomElements::vector3d<ndim,double> velo(&(U_edge[i][2+d_NEQU*k]));
            GeomElements::vector3d<ndim,double> norm = curEdge.normal_vector();
            double Vn = velo.dot_product(norm);
            double VolumePassing = (Vn * curEdge.area());
            if(rC>=0)
            {

                double massPassing = VolumePassing * (Vn>0.0 ? U[i][d_NEQU*lC]: U[i][d_NEQU*rC]);
                vFaceConvectiveTFlux[i][2*k] = turbulenceCellValue[i][2*k] * massPassing;
                vFaceConvectiveTFlux[i][2*k+1] = turbulenceCellValue[i][2*k+1] * massPassing; 
            }
            else
            {
                double massPassing = VolumePassing * U_edge[i][d_NEQU*k];
                vFaceConvectiveTFlux[i][2*k] = vTurbulanceFaceValue[i][2*k] * massPassing;
                vFaceConvectiveTFlux[i][2*k+1] = vTurbulanceFaceValue[i][2*k+1] * massPassing;
            }
        }
    }

}

template<int ndim>
void SSTkomegaModel<ndim>::SolveViscousTFlux()
{
    double txx_k, tyy_k, txx_w, tyy_w;
    double tzz_k, tzz_w;
    double left_f1, right_f1, face_f1, sigma_k,sigma_w, faceMul, faceMut;
    double left_Tvalue[2], right_Tvalue[2];
    GeomElements::vector3d<ndim,double> faceTGrad[2];
    GeomElements::vector3d<ndim,double> left2Right;
    for(int i = 0;i<d_nmesh;i++)
    {
        auto& curBlk = d_hder->blk2D[i];
        for(int k = 0;k<d_hder->nEdges(i);k++)
        {
            auto& curEdge = curBlk.d_localEdges[k];
            auto norm = curEdge.normal_vector();
            int lC = curEdge.lCInd();
            int rC = curEdge.rCInd();
            if(rC>=0)
            {
                left2Right = curBlk.d_localCells[rC].center() - curBlk.d_localCells[lC].center();
                left_f1 = cellF1[i][lC];
                left_TValue[0] = turbulenceCellValue[i][2*lC];
                left_Tvalue[1] = turbulenceCellValue[i][2*lC+1];

                right_f1 = cellF1[i][rC];
                right_Tvalue[0] = turbulenceCellValue[i][2*rc];
                right_Tvalue[1] = turbulenceCellValue[i][2*rC+1];

                face_f1 = 0.5*(left_f1+right_f1);
                faceTGrad[0] = (turbulenceCellGrad[i][2*lC] + turbulenceCellGrad[i][2*rC])*0.5;
                faceTGrad[1] = (turbulenceCellGrad[i][2*lC+1] + turbulenceCellGrad[i][2*rC+1])*0.5;
            }
            else
            {
                left2Right = curEdge.center() - curBlk.d_localCells[lC].center();
                face_f1 = cellF1[i][lC];
                faceTGrad[0] = turbulenceCellGrad[i][2*lC];
                faceTGrad[1] = turbulenceCellGrad[i][2*lC+1];
                
                left_Tvalue[0] = turbulenceCellValue[i][2*lC];
                left_Tvalue[1] = turbulenceCellValue[i][2*lC+1];
                right_Tvalue[0] = vTurbulanceFaceValue[i][2*k];
                right_Tvalue[1] = vTurbulanceFaceValue[i][2*k+1];

            }
            faceTGrad[0] = faceTGrad[0] + left2Right*((right_Tvalue[0] - left_Tvalue[0]-faceTGrad[0].dot_product(left2Right))/left2Right.L2Square());
            faceTGrad[1] = faceTGrad[1] + left2Right*((right_Tvalue[1] - left_Tvalue[1]-faceTGrad[1].dot_product(left2Right))/left2Right.L2Square());

            sigma_k = face_f1 * sigma_k1 + (1-face_f1)*sigma_k2;
            sigma_w = face_f1 * sigma_w1 + (1-face_f1)*sigma_w2;

            faceMul = vFaceMul[i][k];
            faceMut = vFaceMut[i][k];

            txx_k = (faceMul + sigma_k * faceMut) * faceTGrad[0][0];
            tyy_k = (faceMul + sigma_k * faceMut) * faceTGrad[0][1];
            txx_w = (faceMul + sigma_w * faceMut) * faceTGrad[1][0];
            tyy_w = (faceMul + sigma_w * faceMut) * faceTGrad[1][1];
            if(ndim==3)
            {
                tzz_k = (faceMul + sigma_k * faceMut) * faceTGrad[0][2];
                tzz_w = (faceMul + sigma_k * faceMut) * faceTGrad[1][2];
            }
            vFaceViscousTFlux[i][2*k] = txx_k * norm[0] + tyy_k * norm[1];
            vFaceViscousTFlux[i][2*k+1] = txx_w * norm[0] + tyy_w * norm[1];
            if(ndim==3)
            {
                vFaceViscousTFlux[i][2*k] += tzz_k * norm[2];
                vFaceViscousTFlux[i][2*k+1] += tzz_w * norm[2];
            }
            vFaceViscousTFlux[i][2*k] *= curEdge.area();
            vFaceViscousTFlux[i][2*k+1] *= curEdge.area();
        }
    }

}

template<int ndim>
double SSTkomegaModel<ndim>::rans_P_Solver(int meshind, int cellIndex)
{
    GeomElements::vector3d<ndim,double>** gradU = (GeomElements::vector3d<ndim,double>**)d_hder_strategy->getGradientPrimitive();
    double ransP = 0.0;
    auto det_ij = [=](int i ,int j)->int
    {
        if(i ==j)
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }

    double div = gradU[meshind][2+d_NEQU*cellIndex][0] + gradU[meshind][3+d_NEQU*cellIndex][1];
    if(ndim==3)
    {
        div += gradU[meshind][4+d_NEQU*cellIndex][2];
    }
    for(int i = 0;i<ndim;i++)
    {
        for(int j = 0;j<ndim;j++)
        {
            double part_P = vCellMut[meshind][cellIndex] * 
            (gradU[meshind][i+2+d_NEQU*cellIndex][j]+gradU[meshind][j+2+d_NEQU*cellIndex][i]
            - 2.0 / 3.0 *det_ij(i,j)*div) - 2.0 /3.0* det_ij(i,j)*turbulenceConsCellValue[i][2*cellIndex];
            ransP += part_P * gradU[meshind][i+2+d_NEQU*cellIndex][j];
        }
    }

    return ransP;


}

template<int ndim>
void SSTkomegaModel<ndim>::SolveTotalTFlux()
{
    double Qflux[2];
    double** U  =d_hder_strategy->getU();
    for(int i = 0;i<d_nmesh;i++)
    {
        auto& curBlk = d_hder->blk2D[i];
        for(int k = 0;k<d_hder->nCells(i);k++)
        {
            auto& curCell = curBlk.d_localCells[k];
            double f1 = cellF1[i][k];
            double f2 = cellF2[i][k];
            double rdiv = turbulenceCellGrad[i][2*k].dot_product(turbulenceCellGrad[i][2*k+1]);
            double Cw = f1 * Cw1 + (1.0 - f1)* Cw2;
            double beta = f1 * beta1 + (1.0 - f1)* beta2;
            double rho = U[i][d_NEQU*k];

            double rans_P = rans_P_Solver(i,k);

            Qflux[0] = rans_P - _beta * rho * turbulenceCellValue[i][2*k] * turbulenceCellValue[i][2*k+1];
            Qflux[1] = Cw * rho * rans_P / vCellMut[i][k] - beta * rho 
            * turbulenceCellValue[i][2*k+1] * turbulenceCellValue[i][2*k+1] + 2 *(1-f1) *rho * sigma_w2/turbulenceCellValue[i][2*k+1]*rdiv;
            vCellTotalTFlux[i][2*k] = Qflux[0]*curCell.volume();
            vCellTotalTFlux[i][2*k+1] = Qflux[1]*curCell.volume();
        }
        for(int k = 0;k<d_hder->nEdges(i);k++)
        {
            auto& curEdge = curBlk.d_localEdges[k];
            int lC = curEdge.lCInd();
            int rC = curEdge.rCInd();
            vCellTotalTFlux[i][2*lC] += vFaceConvectiveTFlux[2*k] - vFaceViscousTFlux[2*k];
            vCellTotalTFlux[i][2*lC+1] += vFaceConvectiveTFlux[2*k+1] - vFaceViscousTFlux[2*k+1];
            if(rC>=0)
            {
                vCellTotalTFlux[i][2*rC] -= vFaceConvectiveTFlux[2*k] - vFaceViscousTFlux[2*k];
                vCellTotalTFlux[i][2*rC+1] -= vFaceConvectiveTFlux[2*k+1] - vFaceViscousTFlux[2*k+1];
            }
        }
    }

}

// template<int ndim>
// void SSTkomegaModel<ndim>::LUSGSAdvance()
// {
//     double** vDiagOperator = NULL;//Get Operator
//     for(int i = 0; i< d_nmesh;i++)
//     {
//         auto& curBlk = d_hder->blk2D[i];
//         for(int k = 0;k<d_hder->nCells(i);k++)
//         {
//             auto& curCell = curBlk.d_localCells[k];
//             vTurbulanceDiagOperator[i][2*k] = vDiagOperator[i][k] 
//             +2.0* _beta * turbulenceCellValue[i][2*k+1] * curCell.volume();

//             double gradVal =  
//             vTurbulanceDiagOperator[i][2*k+1] = 
//         }
//     }

// }