#include "math.h"


struct FlowParamHolder
{
    double density;
    double pressure;
    double soundSpeed;
    double velocity[3];
};
void computeFlowParams(double ***W, int mesh_ind, int cell_ind, int dim, FlowParamHolder *hder)
{
    const double Gamma = 1.4;
    double density = W[mesh_ind][0][cell_ind];
    double velocity[3] = {0, 0, 0};
    double twice_kineticEnergy = 0;
    for (int i = 0; i < dim; i++)
    {
        velocity[i] = W[mesh_ind][i + 1][cell_ind] / density;
        twice_kineticEnergy += W[mesh_ind][i + 1][cell_ind] * velocity[i];
    }
    double pressure = (Gamma - 1) * (W[mesh_ind][dim + 1][cell_ind] - 0.5*twice_kineticEnergy);
    double soundSpeed = sqrtf64(Gamma * pressure / density);
    hder->density = density;
    hder->soundSpeed = soundSpeed;
    hder->pressure = pressure;
    for (int i = 0; i < dim; i++)
    {
        hder->velocity[i] = velocity[i];
    }
}