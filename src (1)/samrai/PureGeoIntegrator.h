#pragma once

#include "OversetStrategy.h"
#include "PureGeoStrategy.h"
#include "InterpStrategy.h"

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Serializable.h"
//#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
//#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/BaseGridGeometry.h"
//#include "SAMRAI/hier/PatchHierarchy.h"
//#include "SAMRAI/hier/ComponentSelector.h"
//#include "SAMRAI/xfer/RefineAlgorithm.h"
//#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
//#include "SAMRAI/xfer/CoarsenSchedule.h"
//#include "SAMRAI/hier/VariableContext.h"
//#include "SAMRAI/hier/Variable.h"

#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <memory>

//Integrator plays the key role of how a Patch initialize its data when it was created

class PureGeometricIntegrator : public SAMRAI::tbox::Serializable,
                                public SAMRAI::mesh::StandardTagAndInitStrategy
{
    //*In this integrator,we specify:
    //* Variable and Context:
    //* a context and several flag variables for postProcess in Tecplot
    //* a current context for all the other variables
    //* a scratch context (not now because I will test if SAMRAI maintains a "scratch" inside)
    //*
    //* Strategy for geometrical refine(with wall data and over data)
    //* Strategy for gradient detecting
    //*
    //* We Such an integrator should have a interface for TIOGA(We might take it in as a strategy)to read its AMR Patches.
    //*
public:
    enum VariableContextRegistry
    {
        CURRENT_ONLY,
        CURRENT_SCRATCH,
        FLAG_ONLY
    };
    //Constructor and Destructor
    PureGeometricIntegrator(
        const std::string &object_name,
        const SAMRAI::tbox::Dimension &dim,
        const std::shared_ptr<SAMRAI::tbox::Database> &input_db);

    PureGeometricIntegrator(
        const std::string &object_name,
        const SAMRAI::tbox::Dimension &dim,
        const std::shared_ptr<SAMRAI::tbox::Database> &input_db,
        PureGeoStrategy *patch_strategy,
        OversetStrategy *overset_strategy);
    PureGeometricIntegrator(
        const std::string &object_name,
        const SAMRAI::tbox::Dimension &dim,
        const std::shared_ptr<SAMRAI::tbox::Database> &input_db,
        PureGeoStrategy *patch_strategy,
        OversetStrategy *overset_strategy,
        InterpStrategy *interp_strategy);
    virtual ~PureGeometricIntegrator();

    //Interface from Serializable
    void putToRestart(const std::shared_ptr<SAMRAI::tbox::Database> &restart_db) const;

    /*initialize and restart stuff*/
    void getFromInput(const std::shared_ptr<SAMRAI::tbox::Database>, bool is_from_restart);

    void getFromRestart();

    //Interface from StandardTagAndInitStrategy

    void initializeLevelData(
        const std::shared_ptr<SAMRAI::hier::PatchHierarchy> &hierarchy,
        const int level_number,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time,
        const std::shared_ptr<SAMRAI::hier::PatchLevel> &old_level = std::shared_ptr<SAMRAI::hier::PatchLevel>(),
        const bool allocate_data = true);

    void resetHierarchyConfiguration(
        const std::shared_ptr<SAMRAI::hier::PatchHierarchy> &hierarchy,
        const int coarsest_level,
        const int finest_level);

    //End of Interface from StandardTagAndInitStrategy

    //Virtual methods from StangdardTagAndInitStrategy and need implementation
    void applyGradientDetector(
        const std::shared_ptr<SAMRAI::hier::PatchHierarchy> &hierarchy,
        const int ln,
        const double time,
        const int tag_index,
        const bool initial_time,
        const bool uses_richardson_extrapolation_too);
    //Method called by Strategy to addVariable to Patches
    void registerVariable(
        const std::shared_ptr<SAMRAI::hier::Variable> &variable,
        const SAMRAI::hier::IntVector &ghosts,
        const std::shared_ptr<SAMRAI::hier::BaseGridGeometry> transfer_geom,
        const std::string &coarsen_name,
        const std::string &refine_name,
        VariableContextRegistry var_type);

    //Do some initialize here
    void initializeIntegrator(const std::shared_ptr<SAMRAI::mesh::GriddingAlgorithm> gridding_alg);

    //Used for extracting all the possible patches for interpolation

    void findOversetCartBlocks(
        const std::shared_ptr<SAMRAI::hier::PatchHierarchy> &hierarchy,
        const int ln,
        const double time);

    void extractCartGrid(
        int* nf, int* qstride, double** qnodein,
         int** idata, double** rdata,
          int* ngridsin, int* qnodesize,
           int***ibl, double ***q,int**block_global_id, int* nblockin);
    //void extractGrid(int* nf, int* qstride, double** qnodein, int** idata, double** rdata, int* ngridsin, int* qnodesize);
    //void fillCartGridData(int *idata, double *rdata, int ngrids);
private:
    std::string d_object_name;
    SAMRAI::tbox::Dimension d_dim;
    PureGeoStrategy *d_patch_strategy;
    OversetStrategy *d_overset_strategy;
    InterpStrategy *d_interp_strategy;
    //Other stuff may be added later
    SAMRAI::hier::ComponentSelector d_boring_one;
    SAMRAI::hier::ComponentSelector d_scratch_sp;
    SAMRAI::hier::ComponentSelector d_flag_inside;

    //Looks like the only refine and Coarse algorithm I need is one
    //from current to current

    //Algorithm for transferring data from coarse level to finer after a regrid
    std::shared_ptr<SAMRAI::xfer::RefineAlgorithm> d_fill_after_regrid;
    //std::shared_ptr<SAMRAI::xfer::RefineAlgorithm> d_bdry_fill_before_advance;
    //std::shared_ptr<SAMRAI::xfer::RefineAlgorithm> d_fill_before_tagging;

    std::shared_ptr<SAMRAI::xfer::CoarsenAlgorithm> d_coarsen_algorithm;
    //std::vector<std::shared_ptr<SAMRAI::xfer::RefineSchedule>> d_bdry_sched_advance;//This is used for advancing solution on all levels
    //std::vector<std::shared_ptr<SAMRAI::xfer::CoarsenSchedule>> d_coarsen_schedule;

    std::shared_ptr<SAMRAI::hier::VariableContext> d_current;
    std::shared_ptr<SAMRAI::hier::VariableContext> d_scratch;
    std::shared_ptr<SAMRAI::hier::VariableContext> d_flag;

    //std::shared_ptr<SAMRAI::hier::VariableCOntext> d_scratched;

    std::list<std::shared_ptr<SAMRAI::hier::Variable>> d_field_variable;

    int max_level = -1;
    friend class InterpStrategy;

public:
    void set_max_level(int imax)
    {
        max_level = imax;
    }

    void checking(
        const std::shared_ptr<SAMRAI::hier::PatchHierarchy> &hierarchy)
    {
        for (int ln = 0; ln < hierarchy->getMaxNumberOfLevels() - 1; ln++)
        {
            std::shared_ptr<SAMRAI::hier::PatchLevel> level(
                hierarchy->getPatchLevel(ln));
            for (SAMRAI::hier::PatchLevel::iterator ip(level->begin()); ip != level->end(); ip++)
            {
                const std::shared_ptr<SAMRAI::hier::Patch> &p = *ip;
                d_overset_strategy->checking(*p, ln, max_level);
            }
        }
    }
};
