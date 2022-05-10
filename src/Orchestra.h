#include "UniTioga/tioga.h"
#include "reader/CobaltReader.h"
#include "samrai/PureGeoIntegrator.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "samrai/SamraiWrapper.h"
#include <memory>
#include <numeric>

class Orchestra
{
private:
    CobaltReader *rd = NULL;
    
    TIOGA::tioga *tg = NULL;

    WrapSAMRAI* samrai = NULL;
    std::shared_ptr<PureGeometricIntegrator> samrai_integrator = NULL;
    int cur_proc = -1;
    int num_proc = -1;
    //MPI_Comm scomm;
    
public:
    //In the constructor, reader, samrai and tioga's init is done and
    Orchestra(int argc, char *argv[]);
    //In the destructor, samrai reader and tioga's destructor is called
    ~Orchestra();
    //Something COnstructor wouldn't do should be put here.
    void init();
    //Time advancing and Regridding if necessary
    void singleStep(int &cur_step, double &cur_time);
    void run(int step, double time);
    //Something Destructor wouldn't do should be put here.
    void finalize();
    
    void performGeometricRefineOverset();

    void dumpAll(double current_time);
    //For now only dump the flag.
    void dumpUnstruct(const std::string& filename,double current_time);

    void dumpCartesian(const std::string& filename,double current_time);
};
