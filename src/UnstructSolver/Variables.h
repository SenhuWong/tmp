#include "vector3d.h"
const double Gamma = 1.4;

template<int ndim>
class ConservativeVariable;

template <int ndim>
class BasicVariable
{
private:
    double d_rho = 0;
    double d_pressure = 0;
    GeomElements::vector3d<ndim, double> d_velocity;

public:
    double density()
    {
        return d_rho;
    }
    double pressure()
    {
        return d_pressure;
    }
    double velocity(int ind)
    {
        return d_velocity[ind];
    }
    GeomElements::vector3d<ndim, double> velocity()
    {
        return d_velocity;
    }
    BasicVariable(double rho, double p, const GeomElements::vector3d<ndim, double> &velocity)
        : d_rho(rho), d_pressure(p), d_velocity(velocity) {}

    BasicVariable &operator=(const BasicVariable<ndim> &another)
    {
        d_rho = another.d_rho;
        d_pressure = another.d_pressure;
        d_velocity = another.d_velocity;
        return *this;
    }

    ConservativeVariable<ndim> getConservativeVariable()
    {
        double VMSquare = 0;
        for (int i = 0; i < ndim; i++)
        {
            VMSquare += d_velocity[i] * d_velocity[i];
        }
        double rhoE = d_pressure / (Gamma - 1) + 0.5 * d_rho * VMSquare;
        return ConservativeVariable<ndim>(d_rho, rhoE, d_velocity * d_rho);
    }
};

template <int ndim>
class ConservativeVariable
{
private:
    double d_rho;
    double d_rhoE;
    GeomElements::vector3d<ndim, double> d_moment;

public:
    double density()
    {
        return d_rho;
    }
    double totalEnergy()
    {
        return d_rhoE;
    }
    double moment(int ind)
    {
        return d_moment[ind];
    }
    GeomElements::vector3d<ndim, double> moment()
    {
        return d_moment;
    }
    ConservativeVariable(double rho, double rho_E, const GeomElements::vector3d<ndim, double> &moment)
        : d_rho(rho), d_rhoE(rho_E), d_moment(moment) {}

    ConservativeVariable &
    operator=(const ConservativeVariable &another)
    {
        d_rho = another.d_rho;
        d_rhoE = another.d_rhoE;
        d_moment = another.d_moment;
        return *this;
    }
    BasicVariable<ndim> getBasicVariable()
    {
        GeomElements::vector3d<ndim, double> velocity = d_moment / d_rho;
        std::cout<<velocity.L2Square()<<'\n';
        double pressure = (Gamma - 1) * (d_rhoE - 0.5 * d_rho * velocity.L2Square());
        return BasicVariable<ndim>(d_rho,pressure,velocity);
    }
};