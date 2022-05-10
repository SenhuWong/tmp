#pragma once
#include "FluxStrategy.h"
#include "Euler2D.h"

// Given 4 conservatives, calculate the pressure and sound speed
struct FlowParamHolder
{
    double density;
    double pressure;
    double soundSpeed;
    double velocity[3];
};
void computeFlowParams(double ***W, int mesh_ind, int cell_ind, int dim, FlowParamHolder *hder);
class HLLCFluxStrategy : public FluxStrategy
{
public:
    HLLCFluxStrategy(UnstructTopologyHolder *hder, Euler2D *hder_strategy)
        : FluxStrategy(hder, hder_strategy)
    {
    }
    void computeFlux1()
    {
        // if (cur_proc == 0)
        // std::cout<<"ComputeFlux1 being called\n";
        const double Gamma = 1.4;
        double ***UL = d_hder_strategy->getUL();
        double ***UR = d_hder_strategy->getUR();
        double ***Flux = d_hder_strategy->getFluxEdge();
        
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
                    Flux[i][0][k] = 0 * curEdge.area();
                    Flux[i][1][k] = 0 * curEdge.area();
                    double PatFace = UL[i][1][k];
                    for (int l = 0; l < d_dim; l++)
                    {
                        Flux[i][l + 2][k] = PatFace * curEdge.normal_vector()[l] * curEdge.area();
                    }
                }
                else if (rC == GeomElements::edge3d<2>::BoundaryType::FARFIELD)
                {
                    double DatFace = UL[i][0][k];
                    double PatFace = UL[i][1][k];
                    GeomElements::vector3d<2, double> Velocity;
                    for (int l = 0; l < d_dim; l++)
                    {
                        Velocity[l] = UL[i][l + 2][k];
                    }
                    double rhoEatFace = PatFace / (Gamma - 1) + 0.5 * DatFace * (Velocity.L2Square());

                    double normalVelocity = Velocity.dot_product(curEdge.normal_vector());
                    Flux[i][0][k] = normalVelocity * DatFace * curEdge.area();
                    Flux[i][1][k] = normalVelocity * (rhoEatFace + PatFace) * curEdge.area();
                    for (int l = 0; l < d_dim; l++)
                    {
                        Flux[i][l + 2][k] = (normalVelocity * DatFace * Velocity[l] + PatFace * curEdge.normal_vector()[l]) * curEdge.area();
                    }
                }
                else // Is an inner edge, use HLLC
                {
                    // First compute enthalpy
                    double density_left = UL[i][0][k];
                    double pressure_left = UL[i][1][k];
                    GeomElements::vector3d<2, double> velocity_left(UL[i][2][k], UL[i][3][k]);
                    double Vn_left = velocity_left.dot_product(curEdge.normal_vector());
                    double soundspeed_left = sqrtf64(Gamma * pressure_left / density_left);
                    double enthalpy_left = pressure_left * Gamma / ((Gamma - 1) * density_left) + 0.5 * velocity_left.L2Square();

                    double density_righ = UR[i][0][k];
                    double pressure_righ = UR[i][1][k];
                    GeomElements::vector3d<2, double> velocity_righ(UR[i][2][k], UR[i][3][k]);
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
                    if (std::abs(p_star_scratch - p_star) < 1e-10)
                    {
                    }
                    else
                    {
                        std::cout << p_star_scratch << "," << p_star << '\n';
                        throw std::runtime_error("SOmething wrong at p star\n");
                    }
                    if (Sl > 0.0)
                    {
                        double DatFace = UL[i][0][k];
                        double PatFace = UL[i][1][k];
                        double rhoEatFace = PatFace / (Gamma - 1) + 0.5 * DatFace * velocity_left.L2Square();

                        Flux[i][0][k] = Vn_left * DatFace * curEdge.area();
                        Flux[i][1][k] = Vn_left * (rhoEatFace + PatFace) * curEdge.area();
                        for (int l = 0; l < d_dim; l++)
                        {
                            Flux[i][l + 2][k] = (Vn_left * DatFace * velocity_left[l] + PatFace * curEdge.normal_vector()[l]) * curEdge.area();
                        }
                    }
                    else if (Sr < 0.0)
                    {
                        double DatFace = UR[i][0][k];
                        double PatFace = UR[i][1][k];
                        double rhoEatFace = PatFace / (Gamma - 1) + 0.5 * DatFace * velocity_righ.L2Square();

                        Flux[i][0][k] = Vn_righ * DatFace * curEdge.area();
                        Flux[i][1][k] = Vn_righ * (rhoEatFace + PatFace) * curEdge.area();
                        for (int l = 0; l < d_dim; l++)
                        {
                            Flux[i][l + 2][k] = (Vn_righ * DatFace * velocity_righ[l] + PatFace * curEdge.normal_vector()[l]) * curEdge.area();
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
                        Flux[i][0][k] = (SRef - Vn_ref) * density_ref / (SRef - Sm);
                        Flux[i][1][k] = ((SRef - Vn_ref) * rhoE_ref - pressure_ref * Vn_ref + p_star * Sm) / (SRef - Sm);
                        for (int l = 0; l < d_dim; l++)
                        {
                            Flux[i][l + 2][k] = ((SRef - Vn_ref) * density_ref * velocity_ref[l] + (p_star - pressure_ref) * curEdge.normal_vector()[l]) / (SRef - Sm);
                        }
                        Flux[i][0][k] = Flux[i][0][k] * Sm * curEdge.area();
                        Flux[i][1][k] = (Flux[i][1][k] + p_star) * Sm * curEdge.area();
                        for (int l = 0; l < d_dim; l++)
                        {
                            Flux[i][l + 2][k] = (Flux[i][l + 2][k] * Sm + p_star * curEdge.normal_vector()[l]) * curEdge.area();
                        }
                    }
                }
            }
        }
    }

    void computeFlux()
    {
        const double Gamma = 1.4;
        double ***UL = d_hder_strategy->getUL();
        double ***UR = d_hder_strategy->getUR();
        double ***Flux = d_hder_strategy->getFluxEdge();
        FlowParamHolder curParamHolder[2];

        for (int i = 0; i < d_nmesh; i++)
        {
            auto &curBlk = d_hder->blk2D[i];
            for (int k = 0; k < d_hder->nEdges(i); k++)
            {
                auto &curEdge = d_hder->blk2D[i].d_localEdges[k];
                // If it is a natural boundary, I just use the edge variable to compute the flux
                if (curEdge.lCInd() < 0 or curEdge.rCInd() < 0)
                {
                    // std::cout<<"Natural boundary flux computed\n";
                    int rC = curEdge.rCInd();
                    if (rC == GeomElements::edge3d<2>::BoundaryType::WALL) // Wall
                    {
                        Flux[i][0][k] = 0 * curEdge.area();
                        double DatFace = UL[i][0][k];
                        double twice_VMatFace = 0;
                        GeomElements::vector3d<2, double> Velocity;
                        for (int l = 0; l < d_dim; l++)
                        {
                            // UL is the same as W[lc];
                            Velocity[l] = UL[i][l + 1][k] / DatFace;
                            twice_VMatFace += powf64(Velocity[l], 2);
                        }
                        double PatFace = (UL[i][d_dim + 1][k] - 0.5 * DatFace * twice_VMatFace) * (Gamma - 1);
                        for (int l = 0; l < d_dim; l++)
                        {
                            Flux[i][l + 1][k] = PatFace * curEdge.normal_vector()[l] * curEdge.area();
                        }
                        Flux[i][d_dim + 1][k] = 0 * curEdge.area();
                    }
                    else if (rC == GeomElements::edge3d<2>::BoundaryType::FARFIELD) // At Far field
                    {
                        double DatFace = UL[i][0][k];
                        GeomElements::vector3d<2, double> Velocity;
                        double twice_VMatFace = 0;
                        for (int l = 0; l < d_dim; l++)
                        {
                            // UL is the same as Wi[at edge]
                            Velocity[l] = UL[i][l + 1][k] / DatFace;
                            twice_VMatFace += powf64(Velocity[l], 2);
                        }
                        double normalVelocity = Velocity.dot_product(curEdge.normal_vector());
                        double PatFace = (UL[i][d_dim + 1][k] - 0.5 * DatFace * twice_VMatFace) * (Gamma - 1);
                        Flux[i][0][k] = normalVelocity * UL[i][0][k] * curEdge.area();
                        for (int l = 0; l < d_dim; l++)
                        {
                            Flux[i][l + 1][k] = (normalVelocity * UL[i][l + 1][k] + PatFace * curEdge.normal_vector()[l]) * curEdge.area();
                        }
                        Flux[i][d_dim + 1][k] = (normalVelocity * (UL[i][d_dim + 1][k] + PatFace)) * curEdge.area();
                    }
                    else
                    {
                        throw std::runtime_error("Undefined boundary type flux handle\n");
                    }
                    continue;
                }
                FlowParamHolder &leftParams = curParamHolder[0];
                FlowParamHolder &rightParams = curParamHolder[1];
                computeFlowParams(UL, i, k, 2, &(leftParams));
                computeFlowParams(UR, i, k, 2, &(rightParams));
                // Pre Pressure Estimate
                double den_scratch = 0.5 * (leftParams.density + rightParams.density);
                double soundSpeed_scratch = 0.5 * (leftParams.soundSpeed + rightParams.soundSpeed);
                double leftNormalComponent;
                double rightNormalComponent;
                if (d_dim == 2)
                {
                    GeomElements::vector3d<2, double> curNormal = curEdge.normal_vector();
                    GeomElements::vector3d<2, double> leftVelocity(leftParams.velocity[0], leftParams.velocity[1]);
                    GeomElements::vector3d<2, double> rightVelocity(rightParams.velocity[0], rightParams.velocity[1]);
                    leftNormalComponent = curNormal.dot_product(leftVelocity);
                    rightNormalComponent = curNormal.dot_product(rightVelocity);
                }
                // Pressure Estimate
                double p_pvrs = 0.5 * (leftParams.pressure + rightParams.pressure) - 0.5 * (rightNormalComponent - leftNormalComponent) * den_scratch * soundSpeed_scratch;
                double p_star = std::max(0.0, p_pvrs);
                // Wave Speed Estimate
                double Sl, Sr;
                if (leftParams.pressure >= p_star)
                {
                    Sl = leftNormalComponent - leftParams.soundSpeed;
                }
                else
                {
                    Sl = leftNormalComponent - leftParams.soundSpeed * sqrtf64(1 + (p_star / leftParams.pressure - 1) * (Gamma + 1) / (2 * Gamma));
                }
                if (rightParams.pressure >= p_star)
                {
                    Sr = rightNormalComponent + rightParams.soundSpeed;
                }
                else
                {
                    Sr = rightNormalComponent + rightParams.soundSpeed * sqrtf64(1 + (p_star / rightParams.pressure - 1) * (Gamma + 1) / (2 * Gamma));
                }
                double S_star = (rightParams.pressure - leftParams.pressure + leftParams.density * leftNormalComponent * (Sl - leftNormalComponent) - rightParams.density * rightNormalComponent * (Sr - rightNormalComponent)) / (leftParams.density * (Sl - leftNormalComponent) - rightParams.density * (Sr - rightNormalComponent));
                // HLLC flux computation
                if (k < curBlk.d_nFluxedEs)
                {
                    if (Sl > S_star or S_star > Sr or Sl > Sr)
                    {
                        std::cout << "Something rong\n";
                        std::cin.get();
                    }
                }

                if (Sl >= 0) // Fl
                {
                    Flux[i][0][k] = (leftNormalComponent * UL[i][0][k]) * curEdge.area();
                    for (int l = 0; l < d_dim; l++)
                    {
                        Flux[i][l + 1][k] = (leftNormalComponent * UL[i][l + 1][k] + leftParams.pressure * curEdge.normal_vector()[l]) * curEdge.area();
                    }
                    Flux[i][d_dim + 1][k] = leftNormalComponent * (UL[i][d_dim + 1][k] + leftParams.pressure) * curEdge.area();
                }
                else if (S_star >= 0 and Sl < 0) // F*L = FL + SL*(U*L - UL)
                {

                    double common_Efficient = leftParams.density * (Sl - leftNormalComponent) / (Sl - S_star);
                    GeomElements::vector3d<2, double> velocity(leftParams.velocity[0], leftParams.velocity[1]);
                    GeomElements::vector3d<2, double> corrected_velocity = velocity + curEdge.normal_vector() * (S_star - velocity.dot_product(curEdge.normal_vector()));

                    Flux[i][0][k] = (leftNormalComponent * UL[i][0][k] + Sl * (common_Efficient - UL[i][0][k])) * curEdge.area();
                    for (int l = 0; l < d_dim; l++)
                    {
                        Flux[i][l + 1][k] = (leftNormalComponent * UL[i][l + 1][k] + leftParams.pressure * curEdge.normal_vector()[l] + Sl * (common_Efficient * (corrected_velocity[l]) - UL[i][l + 1][k])) * curEdge.area();
                    }
                    Flux[i][d_dim + 1][k] = (leftNormalComponent * (UL[i][d_dim + 1][k] + leftParams.pressure) + Sl * (common_Efficient * (UL[i][d_dim + 1][k] / leftParams.density + (S_star - leftNormalComponent) * (S_star + leftParams.pressure / (leftParams.density * (Sl - leftNormalComponent)))) - UL[i][d_dim + 1][k])) * curEdge.area();
                }
                else if (Sr > 0 and S_star < 0) // F*R
                {
                    double common_Efficient = rightParams.density * (Sr - rightNormalComponent) / (Sr - S_star);
                    GeomElements::vector3d<2, double> velocity(rightParams.velocity[0], rightParams.velocity[1]);
                    GeomElements::vector3d<2, double> corrected_velocity = velocity + curEdge.normal_vector() * (S_star - velocity.dot_product(curEdge.normal_vector()));
                    Flux[i][0][k] = (rightNormalComponent * UR[i][0][k] + Sr * (common_Efficient - UR[i][0][k])) * curEdge.area();
                    for (int l = 0; l < d_dim; l++)
                    {
                        Flux[i][l + 1][k] = (rightNormalComponent * UR[i][l + 1][k] + rightParams.pressure * curEdge.normal_vector()[l] + Sr * (common_Efficient * (corrected_velocity[l]) - UR[i][l + 1][k])) * curEdge.area();
                    }
                    Flux[i][d_dim + 1][k] = (rightNormalComponent * (UR[i][d_dim + 1][k] + rightParams.pressure) + Sr * (common_Efficient * (UR[i][d_dim + 1][k] / rightParams.density + (S_star - rightNormalComponent) * (S_star + rightParams.pressure / (rightParams.density * (Sr - rightNormalComponent)))) - UR[i][d_dim + 1][k])) * curEdge.area();
                }
                else if (Sr <= 0)
                {
                    Flux[i][0][k] = (rightNormalComponent * UR[i][0][k]) * curEdge.area();
                    for (int l = 0; l < d_dim; l++)
                    {
                        Flux[i][l + 1][k] = (rightNormalComponent * UR[i][l + 1][k] + rightParams.pressure * curEdge.normal_vector()[l]) * curEdge.area();
                    }
                    Flux[i][d_dim + 1][k] = rightNormalComponent * (UR[i][d_dim + 1][k] + rightParams.pressure) * curEdge.area();
                }
                else
                {
                    if (k < curBlk.d_nFluxedEs)
                    {
                        std::cout << curBlk.d_nFluxedEs << '\n';
                        std::cout << "Something wrong with wave speed estimation\n";
                        std::cout << Sl << " " << Sr << '\n';
                        std::cout << "leftNormalComponent " << leftNormalComponent << " leftParams.soundSpeed :" << leftParams.soundSpeed << ",leftParams.pressure :" << leftParams.pressure << '\n';
                        std::cout << " rightNormalComponent " << rightNormalComponent << " rightParams.soundSpeed:" << rightParams.soundSpeed << ",rightParams.pressure:" << rightParams.pressure << '\n';
                        std::cout << "On proc " << cur_proc << " 's " << i << "th mesh " << k << " th Edge" << '\n';

                        std::cin.get();
                    }
                }
            }
        }
    }
};