// #include "FluxStrategy.h"
// #include "Euler2D.h"

// struct FlowParamHolder
// {
//     double density;
//     double pressure;
//     double soundSpeed;
//     double velocity[3];
//     double enthalpy;
// };

// class newHLLCStrategy : public FluxStrategy
// {
//     newHLLCStrategy(UnstructTopologyHolder *hder, Euler2D *hder_strategy)
//         : FluxStrategy(d_hder, hder_strategy)
//     {
//     }

//     void computeFlux()
//     {
//         const double Gamma = 1.4;
//         double ***WL = d_hder_strategy->getWL();
//         double ***WR = d_hder_strategy->getWR();
//         double ***Flux = d_hder_strategy->getFi();
//         FlowParamHolder curParamHolder[2];

//         for (int i = 0; i < d_nmesh; i++)
//         {
//             auto &curBlk = d_hder->blk2D[i];
//             for (int k = 0; k < d_hder->nEdges(i); k++)
//             {
//                 auto &curEdge = curBlk.d_localEdges[k];
//                 int rC = curEdge.rCInd();
//                 int lC = curEdge.lCInd();
//                 if (rC >= 0)
//                 {
//                     FlowParamHolder &leftParams = curParamHolder[0];
//                     FlowParamHolder &rightParams = curParamHolder[1];
//                     // Compute leftParam
//                     leftParams.density = WL[i][0][k];
//                     double SquareVMatFace = 0;
//                     for (int l = 0; l < d_dim; l++)
//                     {
//                         leftParams.velocity[l] = WL[i][l + 1][k] / leftParams.density;
//                         SquareVMatFace += powf64(leftParams.velocity[l], 2);
//                     }
//                     leftParams.pressure = (Gamma - 1) * (WL[i][d_dim + 1][k] - 0.5 * leftParams.density * SquareVMatFace);
//                     leftParams.soundSpeed = sqrtf64(Gamma * leftParams.pressure / leftParams.density);
//                     leftParams.enthalpy = (WL[i][d_dim + 1][k] + leftParams.pressure) / leftParams.density;
//                     // Compute rightParam
//                     rightParams.density = WR[i][0][k];
//                     SquareVMatFace = 0;
//                     for (int l = 0; l < d_dim; l++)
//                     {
//                         rightParams.velocity[l] = WR[i][l + 1][k] / rightParams.density;
//                         SquareVMatFace += powf64(rightParams.velocity[l], 2);
//                     }
//                     rightParams.pressure = (Gamma - 1) * (WR[i][d_dim + 1][k] - 0.5 * leftParams.density * SquareVMatFace);
//                     rightParams.soundSpeed = sqrtf64(Gamma * rightParams.pressure / rightParams.density);
//                     rightParams.enthalpy = (WR[i][d_dim + 1][k] + rightParams.pressure) / rightParams.density;
//                     // Compute Roe Average
//                     double Vn_RoeAverage;
//                     double VnL, VnR;
//                     double C_RoeAverage;
//                     double H_RoeAverage;
//                     double eff_left = sqrtf64(leftParams.density) / (sqrtf64(leftParams.density) + sqrtf64(rightParams.density));
//                     double eff_right = 1 - eff_left;
//                     if (d_dim == 2)
//                     {
//                         GeomElements::vector3d<2, double> velocityLeft(leftParams.velocity[0], leftParams.velocity[1]);
//                         GeomElements::vector3d<2, double> velocityRight(rightParams.velocity[0], rightParams.velocity[1]);
//                         GeomElements::vector3d<2, double> velocityRoeAver = velocityLeft * eff_left + velocityRight * eff_right;
//                         VnL = velocityLeft.dot_product(curEdge.normal_vector());
//                         VnR = velocityRight.dot_product(curEdge.normal_vector());
//                         Vn_RoeAverage = velocityRoeAver.dot_product(curEdge.normal_vector());

//                         H_RoeAverage = leftParams.enthalpy * eff_left + rightParams.enthalpy * eff_right;
//                         C_RoeAverage = sqrtf64((Gamma - 1) * (H_RoeAverage - 0.5 * powf64(velocityRoeAver.normalize(), 2)));
//                     }
//                     double Sl = std::min(VnL - C_RoeAverage, Vn_RoeAverage - C_RoeAverage);
//                     double Sr = std::max(Vn_RoeAverage + C_RoeAverage, VnR + C_RoeAverage);
//                     double SM = (rightParams.density*VnR*(Sr-VnR) - leftParams.density*VnL*(Sl-VnL) + leftParams.pressure - rightParams.pressure)/(rightParams.density*(Sr-VnR) - leftParams.density*(Sl-VnL));
//                     if(Sl>SM or SM>Sr)
//                     {
//                         throw std::runtime_error("Error at wave speed estimation\n");
//                     }
//                     if(Sl>0.0)
//                     {
                        

//                     }
//                     else if(Sr<0.0)
//                     {

//                     }
//                     else if(SM>0.0)
//                     {

//                     }              
//                     else if(SM<=0.0)
//                     {

//                     }
//                     else
//                     {
//                         throw std::runtime_error("Error at flux computation\n");

//                     }

//                 }
//                 else if (rC == GeomElements::edge3d<2>::BoundaryType::WALL)
//                 {
//                     Flux[i][0][k] = 0 * curEdge.area();
//                     double DatFace = WL[i][0][k];
//                     double SquareVMatFace = 0;
//                     GeomElements::vector3d<2, double> Velocity;
//                     for (int l = 0; l < d_dim; l++)
//                     {
//                         Velocity[l] = WL[i][l + 1][k] / DatFace;
//                         SquareVMatFace += powf64(Velocity[l], 2);
//                     }
//                     double PatFace = (Gamma - 1) * (WL[i][d_dim + 1][k] - 0.5 * DatFace * SquareVMatFace);
//                     for (int l = 0; l < d_dim; l++)
//                     {
//                         Flux[i][l + 1][k] = PatFace * curEdge.normal_vector()[l] * curEdge.area();
//                     }
//                     Flux[i][d_dim + 1][k] = 0 * curEdge.area();
//                 }
//                 else if (rC == GeomElements::edge3d<2>::BoundaryType::FARFIELD)
//                 {
//                     double DatFace = WL[i][0][k];
//                     GeomElements::vector3d<2, double> Velocity;
//                     double SquareVMatFace = 0;
//                     for (int l = 0; l < d_dim; l++)
//                     {
//                         Velocity[l] = WL[i][l + 1][k] / DatFace;
//                         SquareVMatFace += powf64(Velocity[l], 2);
//                     }
//                     double normalVelocity = Velocity.dot_product(curEdge.normal_vector());
//                     double PatFace = (Gamma - 1) * (WL[i][d_dim + 1][k] - 0.5 * DatFace * SquareVMatFace);
//                     Flux[i][0][k] = normalVelocity * WL[i][0][k] * curEdge.area();
//                     for (int l = 0; l < d_dim; l++)
//                     {
//                         Flux[i][l + 1][k] = (normalVelocity * WL[i][l + 1][k] + PatFace * curEdge.normal_vector()[l]) * curEdge.area();
//                     }
//                     Flux[i][d_dim + 1][k] = (normalVelocity * (WL[i][d_dim + 1][k] + PatFace)) * curEdge.area();
//                 }
//                 else
//                 {
//                     throw std::runtime_error("Undefined Boundary met at HLLCFlux\n");
//                 }
//             }
//         }
//     }
// };