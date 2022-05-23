#include "HLLCFluxStrategy.h"
#include "UnstructIntegrator.h"
void HLLCFluxStrategy::computeFlux()
{
    updateLeftRight();
    // if (cur_proc == 0)
    // std::cout<<"ComputeFlux1 being called\n";
    const double Gamma = 1.4;
    double Flux[5];
    double **Residual = d_hder_strategy->getResidual();
    
    for (int i = 0; i < d_nmesh; i++)
    {
        auto &curBlk = d_hder->blk2D[i];
        for (int k = 0; k < d_hder->nEdges(i); k++)
        {
            auto &curEdge = curBlk.d_localEdges[k];
            int lC = curEdge.lCInd();
            int rC = curEdge.rCInd();
            if (lC < 0)
            {
                throw std::runtime_error("Left cell can't be boundary\n");
            }
            if (rC == GeomElements::edge3d<2>::BoundaryType::WALL)
            {
                Flux[0] = 0 * curEdge.area();
                Flux[1] = 0 * curEdge.area();
                double PatFace = UL[i][1+d_NEQU*k];
                for (int l = 0; l < d_dim; l++)
                {
                    Flux[l + 2] = PatFace * curEdge.normal_vector()[l] * curEdge.area();
                }
                
            }
            else if (rC == GeomElements::edge3d<2>::BoundaryType::FARFIELD)
            {
                double DatFace = UL[i][0+d_NEQU*k];
                double PatFace = UL[i][1+d_NEQU*k];
                GeomElements::vector3d<2, double> Velocity;
                for (int l = 0; l < d_dim; l++)
                {
                    Velocity[l] = UL[i][l + 2+d_NEQU*k];
                }
                double rhoEatFace = PatFace / (Gamma - 1) + 0.5 * DatFace * (Velocity.L2Square());
                double normalVelocity = Velocity.dot_product(curEdge.normal_vector());
                Flux[0] = normalVelocity * DatFace * curEdge.area();
                Flux[1] = normalVelocity * (rhoEatFace + PatFace) * curEdge.area();
                for (int l = 0; l < d_dim; l++)
                {
                    Flux[l + 2] = (normalVelocity * DatFace * Velocity[l] + PatFace * curEdge.normal_vector()[l]) * curEdge.area();
                }
            }
            else // Is an inner edge, use HLLC
            {
                // First compute enthalpy
                double density_left = UL[i][0+d_NEQU*k];
                double pressure_left = UL[i][1+d_NEQU*k];
                GeomElements::vector3d<2, double> velocity_left(UL[i][2+d_NEQU*k], UL[i][3+d_NEQU*k]);
                double Vn_left = velocity_left.dot_product(curEdge.normal_vector());
                double soundspeed_left = sqrtf64(Gamma * pressure_left / density_left);
                double enthalpy_left = pressure_left * Gamma / ((Gamma - 1) * density_left) + 0.5 * velocity_left.L2Square();
                double density_righ = UR[i][0+d_NEQU*k];
                double pressure_righ = UR[i][1+d_NEQU*k];
                GeomElements::vector3d<2, double> velocity_righ(UR[i][2+d_NEQU*k], UR[i][3+d_NEQU*k]);
                double Vn_righ = velocity_righ.dot_product(curEdge.normal_vector());
                double soundspeed_righ = sqrtf64(Gamma * pressure_righ / density_righ);
                double enthalpy_righ = pressure_righ * Gamma / ((Gamma - 1) * density_righ) + 0.5 * velocity_righ.L2Square();
                double efficient_left = sqrtf64(density_left) / (sqrtf64(density_left) + sqrtf64(density_righ));
                double efficient_righ = 1 - efficient_left;
                double enthalpy_tilda = efficient_left * enthalpy_left + efficient_righ * enthalpy_righ;
                double Vn_tilda = efficient_left * Vn_left + efficient_righ * Vn_righ;
                double soundspeed_tilda = sqrtf64((Gamma - 1) * (enthalpy_tilda - 0.5 * powf64(Vn_tilda, 2)));
                double Sl = std::min(Vn_left - soundspeed_left, Vn_tilda - soundspeed_tilda);
                double Sr = std::max(Vn_tilda + soundspeed_tilda, Vn_righ + soundspeed_righ);
                double Sm = (density_righ * Vn_righ * (Sr - Vn_righ) - density_left * Vn_left * (Sl - Vn_left) + pressure_left - pressure_righ) / (density_righ * (Sr - Vn_righ) - density_left * (Sl - Vn_left));
                double p_star = density_left * (Vn_left - Sl) * (Vn_left - Sm) + pressure_left;
                double p_star_scratch = density_righ * (Vn_righ - Sr) * (Vn_righ - Sm) + pressure_righ;
                if (std::abs(p_star_scratch - p_star) < 1e-10 or k < curBlk.d_nFluxedEs)
                {
                }
                else
                {
                    std::cout << p_star_scratch << "," << p_star << '\n';
                    throw std::runtime_error("SOmething wrong at p star\n");
                }
                if (Sl > 0.0)
                {
                    double DatFace = UL[i][0+d_NEQU*k];
                    double PatFace = UL[i][1+d_NEQU*k];
                    double rhoEatFace = PatFace / (Gamma - 1) + 0.5 * DatFace * velocity_left.L2Square();
                    Flux[0] = Vn_left * DatFace * curEdge.area();
                    Flux[1] = Vn_left * (rhoEatFace + PatFace) * curEdge.area();
                    for (int l = 0; l < d_dim; l++)
                    {
                        Flux[l + 2] = (Vn_left * DatFace * velocity_left[l] + PatFace * curEdge.normal_vector()[l]) * curEdge.area();
                    }
                }
                else if (Sr < 0.0)
                {
                    double DatFace = UR[i][0+d_NEQU*k];
                    double PatFace = UR[i][1+d_NEQU*k];
                    double rhoEatFace = PatFace / (Gamma - 1) + 0.5 * DatFace * velocity_righ.L2Square();
                    Flux[0] = Vn_righ * DatFace * curEdge.area();
                    Flux[1] = Vn_righ * (rhoEatFace + PatFace) * curEdge.area();
                    for (int l = 0; l < d_dim; l++)
                    {
                        Flux[l + 2] = (Vn_righ * DatFace * velocity_righ[l] + PatFace * curEdge.normal_vector()[l]) * curEdge.area();
                    }
                }
                else // Between SM and SR
                {
                    double SRef;
                    double Vn_ref;
                    double density_ref;
                    double pressure_ref;
                    GeomElements::vector3d<2, double> velocity_ref;
                    if (Sm < 0.0)
                    {
                        SRef = Sr;
                        Vn_ref = Vn_righ;
                        density_ref = density_righ;
                        pressure_ref = pressure_righ;
                        velocity_ref = velocity_righ;
                    }
                    else
                    {
                        SRef = Sl;
                        Vn_ref = Vn_left;
                        density_ref = density_left;
                        pressure_ref = pressure_left;
                        velocity_ref = velocity_left;
                    }
                    double rhoE_ref = pressure_ref / (Gamma - 1) + 0.5 * density_ref * velocity_ref.L2Square();
                    Flux[0] = (SRef - Vn_ref) * density_ref / (SRef - Sm);
                    Flux[1] = ((SRef - Vn_ref) * rhoE_ref - pressure_ref * Vn_ref + p_star * Sm) / (SRef - Sm);
                    for (int l = 0; l < d_dim; l++)
                    {
                        Flux[l + 2] = ((SRef - Vn_ref) * density_ref * velocity_ref[l] + (p_star - pressure_ref) * curEdge.normal_vector()[l]) / (SRef - Sm);
                    }
                    Flux[0] = Flux[0] * Sm * curEdge.area();
                    Flux[1] = (Flux[1] + p_star) * Sm * curEdge.area();
                    for (int l = 0; l < d_dim; l++)
                    {
                        Flux[l + 2] = (Flux[l + 2] * Sm + p_star * curEdge.normal_vector()[l]) * curEdge.area();
                    }
                }
            }

            if (lC >= 0) // Left plus right minus
            {
                for (int j = 0; j < d_NEQU; j++)
                {
                    Residual[i][j+d_NEQU*lC] += Flux[j];
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
                    Residual[i][j+d_NEQU*rC] -= Flux[j];
                }
            }
        }
    }
}