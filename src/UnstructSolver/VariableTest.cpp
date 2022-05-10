#include "Variables.h"

int main()
{
    // Conservative Variables' function
    //      Constructor
    ConservativeVariable<2> *cons = new ConservativeVariable<2>(0.8, 3.6113, GeomElements::vector3d<2, double>(0.298476, 0.898283));
    std::cout << cons->density() << "," << cons->totalEnergy() << "," << cons->moment()[0] << "," << cons->moment()[1] << '\n';
    //      Copy Operator
    ConservativeVariable<2> copy_cons = *cons;
    std::cout << copy_cons.density() << "," << copy_cons.totalEnergy() << "," << copy_cons.moment()[0] << "," << copy_cons.moment()[1] << '\n';
    //      getBasic
    BasicVariable<2> bsic = copy_cons.getBasicVariable();
    std::cout << bsic.density() << "," << bsic.pressure() << "," << bsic.velocity()[0] << "," << bsic.velocity()[1] << '\n';
    // Basic Variables' function
    //      Constructor
    BasicVariable<2> construct_bsic = BasicVariable<2>(bsic.density(), bsic.pressure(), bsic.velocity());
    std::cout << construct_bsic.density() << "," << construct_bsic.pressure() << "," << construct_bsic.velocity()[0] << "," << construct_bsic.velocity()[1] << '\n';
    //      Copy Operator
    BasicVariable<2> copy_bsic = bsic;
    std::cout << copy_bsic.density() << "," << copy_bsic.pressure() << "," << copy_bsic.velocity()[0] << "," << copy_bsic.velocity()[1] << '\n';
    //      getConservative
    ConservativeVariable<2> getCons = copy_bsic.getConservativeVariable();
    std::cout<< getCons.density()<<","<<getCons.totalEnergy()<<","<<getCons.moment()[0]<<getCons.moment()[1]<<'\n';
}
