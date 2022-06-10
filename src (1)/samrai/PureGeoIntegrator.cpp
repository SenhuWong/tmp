#include "PureGeoIntegrator.h"
#include "SAMRAI/tbox/RestartManager.h"
//#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/hier/PatchDataRestartManager.h"

PureGeometricIntegrator::PureGeometricIntegrator(
    const std::string &object_name,
    const SAMRAI::tbox::Dimension &dim,
    const std::shared_ptr<SAMRAI::tbox::Database> &input_db)
    : d_object_name(object_name),
      d_current(SAMRAI::hier::VariableDatabase::getDatabase()->getContext("CURRENT")),
      d_scratch(SAMRAI::hier::VariableDatabase::getDatabase()->getContext("SCRATCH")),
      d_flag(SAMRAI::hier::VariableDatabase::getDatabase()->getContext("FLAG")),
      d_dim(dim)
{
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(patch_strategy != 0);

    //Set Context need to be done upon strategy is attached.
    //d_patch_strategy->set_Context(d_current, d_scratch);
    //d_overset_strategy->set_Context(d_flag);

    SAMRAI::tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);
    bool is_from_restart = SAMRAI::tbox::RestartManager::getManager()->isFromRestart();
    if (is_from_restart)
    {
        getFromRestart();
    }
    getFromInput(input_db, is_from_restart);
}

PureGeometricIntegrator::PureGeometricIntegrator(
    const std::string &object_name,
    const SAMRAI::tbox::Dimension &dim,
    const std::shared_ptr<SAMRAI::tbox::Database> &input_db,
    PureGeoStrategy *patch_strategy,
    OversetStrategy *overset_strategy) : d_object_name(object_name),
                                         d_patch_strategy(patch_strategy),
                                         d_overset_strategy(overset_strategy),
                                         d_current(SAMRAI::hier::VariableDatabase::getDatabase()->getContext("CURRENT")),
                                         d_scratch(SAMRAI::hier::VariableDatabase::getDatabase()->getContext("SCRATCH")),
                                         d_flag(SAMRAI::hier::VariableDatabase::getDatabase()->getContext("FLAG")),
                                         d_dim(dim)
{
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(patch_strategy != 0);
    d_patch_strategy->set_Context(d_current, d_scratch);
    d_overset_strategy->set_Context(d_flag);

    SAMRAI::tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);
    bool is_from_restart = SAMRAI::tbox::RestartManager::getManager()->isFromRestart();
    if (is_from_restart)
    {
        getFromRestart();
    }
    getFromInput(input_db, is_from_restart);
}

PureGeometricIntegrator::PureGeometricIntegrator(
    const std::string &object_name,
    const SAMRAI::tbox::Dimension &dim,
    const std::shared_ptr<SAMRAI::tbox::Database> &input_db,
    PureGeoStrategy *patch_strategy,
    OversetStrategy *overset_strategy,
    InterpStrategy *interp_strategy) : d_object_name(object_name),
                                       d_patch_strategy(patch_strategy),
                                       d_overset_strategy(overset_strategy),
                                       d_interp_strategy(interp_strategy),
                                       d_current(SAMRAI::hier::VariableDatabase::getDatabase()->getContext("CURRENT")),
                                       d_scratch(SAMRAI::hier::VariableDatabase::getDatabase()->getContext("SCRATCH")),
                                       d_flag(SAMRAI::hier::VariableDatabase::getDatabase()->getContext("FLAG")),
                                       d_dim(dim)
{
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(patch_strategy != 0);
    d_patch_strategy->set_Context(d_current, d_scratch);
    d_overset_strategy->set_Context(d_flag);
    d_interp_strategy->set_Context(d_current);

    SAMRAI::tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);
    bool is_from_restart = SAMRAI::tbox::RestartManager::getManager()->isFromRestart();
    if (is_from_restart)
    {
        getFromRestart();
    }
    getFromInput(input_db, is_from_restart);
}

void PureGeometricIntegrator::getFromInput(const std::shared_ptr<SAMRAI::tbox::Database> input_db, bool is_from_restart)
{
    //We don't have any efficient to get
}
void PureGeometricIntegrator::getFromRestart()
{
    //We do not have anything to restart
}

void PureGeometricIntegrator::registerVariable(
    const std::shared_ptr<SAMRAI::hier::Variable> &variable,
    const SAMRAI::hier::IntVector &ghosts,
    const std::shared_ptr<SAMRAI::hier::BaseGridGeometry> transfer_geom,
    const std::string &coarsen_name,
    const std::string &refine_name,
    VariableContextRegistry var_type)
{
    std::cout << "PureGeometricIntegrator::Registering Variable";
    //Usually we would have known which variable needs which contexts
    //For now let's assume every variable except flag is registered with current and scratch
    //Later we might need a ENUM to refer to different registeration??
    TBOX_ASSERT(variable);
    TBOX_ASSERT(transfer_geom);
    TBOX_ASSERT_OBJDIM_EQUALITY2(*variable, ghosts);
    SAMRAI::tbox::Dimension dim(ghosts.getDim());
    SAMRAI::hier::VariableDatabase *variable_db = SAMRAI::hier::VariableDatabase::getDatabase();
    const SAMRAI::hier::IntVector no_ghost(dim, 0);
    if (var_type == VariableContextRegistry::FLAG_ONLY)
    {
        //No Operation is registered
        std::cout << "registering\n";
        const int flag = variable_db->registerVariableAndContext(
            variable,
            d_flag,
            no_ghost);
        std::cout << "register successed\n";
        d_flag_inside.setFlag(flag);
        std::cout << "setFlag successed\n";
        SAMRAI::hier::PatchDataRestartManager::getManager()->registerPatchDataForRestart(flag);
    }
    else if (var_type == VariableContextRegistry::CURRENT_ONLY) //For now a CUURENT_ONLY variable is temporary flag but has some real meaning
    {
        //I don't actually know what is expected from this.
        if (!d_fill_after_regrid)
        {
            d_fill_after_regrid.reset(new SAMRAI::xfer::RefineAlgorithm());
            d_coarsen_algorithm.reset(new SAMRAI::xfer::CoarsenAlgorithm(dim));
        }

        const int current = variable_db->registerVariableAndContext(
            variable,
            d_current,
            no_ghost);
    }
    else if (var_type == VariableContextRegistry::CURRENT_SCRATCH)
    {
        if (!d_fill_after_regrid)
        {
            d_fill_after_regrid.reset(new SAMRAI::xfer::RefineAlgorithm());
            d_coarsen_algorithm.reset(new SAMRAI::xfer::CoarsenAlgorithm(dim));
        }
        //PushBack the Variable list
        d_field_variable.push_back(variable);
        //FInd the refine operator and its ghost cell width requirement
        std::shared_ptr<SAMRAI::hier::RefineOperator> refine_operator(
            transfer_geom->lookupRefineOperator(variable, refine_name));
        std::shared_ptr<SAMRAI::hier::CoarsenOperator> coarsen_operator(
            transfer_geom->lookupCoarsenOperator(variable, coarsen_name));
        //Get the unique index from variable_database
        const SAMRAI::hier::IntVector what_i_want = refine_operator->getStencilWidth(dim);
        const int current = variable_db->registerVariableAndContext(
            variable,
            d_current,
            no_ghost);
        const int scratch = variable_db->registerVariableAndContext(
            variable,
            d_scratch,
            what_i_want);

        //Set the component selector.
        d_boring_one.setFlag(current);
        d_scratch_sp.setFlag(scratch);
        //Register for restart
        SAMRAI::hier::PatchDataRestartManager::getManager()->registerPatchDataForRestart(current);
        //Register for Refine and Coarsen
        d_fill_after_regrid->registerRefine(current, current, scratch, refine_operator);
        d_coarsen_algorithm->registerCoarsen(scratch, scratch, coarsen_operator);
    }
}

PureGeometricIntegrator::~PureGeometricIntegrator()
{
}

//Inherited from Serializable
void PureGeometricIntegrator::putToRestart(const std::shared_ptr<SAMRAI::tbox::Database> &restart_db) const
{
    //DO nothing
}
//Interface for TagAndInitStrategy

void PureGeometricIntegrator::initializeLevelData(
    const std::shared_ptr<SAMRAI::hier::PatchHierarchy> &hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    const std::shared_ptr<SAMRAI::hier::PatchLevel> &old_level,
    const bool allocate_data)
{
    std::cout << "PureGeometricIntegrator::InitialzieLevel Data being called\n";
    NULL_USE(can_be_refined);
    NULL_USE(allocate_data);

    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT(hierarchy->getPatchLevel(level_number));
    TBOX_ASSERT(level_number >= 0);

    std::shared_ptr<SAMRAI::hier::PatchLevel> level(
        hierarchy->getPatchLevel(level_number));
    //Now there is only a boring selector which register everything
    level->allocatePatchData(d_boring_one, init_data_time);
    level->allocatePatchData(d_scratch_sp, init_data_time);
    level->allocatePatchData(d_flag_inside, init_data_time);
    //%TODO::Here we might need to specify the type of data transfer then implement
    //Since this is only pureGeometric,we want only to
    if (level_number > 0 or old_level)
    {
        d_fill_after_regrid->createSchedule(
                               level,            //dst_level
                               old_level,        //src_level
                               level_number - 1, //next_coarser level
                               hierarchy,
                               d_patch_strategy)
            ->fillData(init_data_time);
    }
    level->deallocatePatchData(d_scratch_sp);

    //If it is initial time ,then do initialize ,else the fillData above is enough
    for (SAMRAI::hier::PatchLevel::iterator p(level->begin()); p != level->end(); p++)
    {
        const std::shared_ptr<SAMRAI::hier::Patch> &patch = *p;
        d_patch_strategy->initializeDataOnPatch(*patch, init_data_time, initial_time);
        d_overset_strategy->tagOverBoundaryCrossCells(*patch, init_data_time);
    }
}

void PureGeometricIntegrator::initializeIntegrator(const std::shared_ptr<SAMRAI::mesh::GriddingAlgorithm> gridding_alg)
{

    NULL_USE(gridding_alg);
    TBOX_ASSERT(gridding_alg);
    d_patch_strategy->registerModelVariables(this);
    d_overset_strategy->registerModelVariables(this);
    d_interp_strategy->registerModelVariables(this);

    //Do something with d_patch_strategy
}

void PureGeometricIntegrator::resetHierarchyConfiguration(
    const std::shared_ptr<SAMRAI::hier::PatchHierarchy> &hierarchy,
    const int coarsest_level,
    const int finest_level)
{
    //DO Something to the Schedule
    int finest_hiera_level = hierarchy->getFinestLevelNumber();
    //If the level number is changed, then the schedule should have a resize

    //BUild Coarsen or refine communication schedule
}

void PureGeometricIntegrator::applyGradientDetector(
    const std::shared_ptr<SAMRAI::hier::PatchHierarchy> &hierarchy,
    const int ln,
    const double time,
    const int tag_index, //maintained by GriddingALgorithm using StandardTagANdINitialize
    const bool initial_time,
    const bool uses_richardson_extrapolation_too)
{
    std::cout << "applyGradientDetector being called\n";
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT(hierarchy->getPatchLevel(ln));
    //std::cout << "level" << ln << " " << is_finest << '\n';

    std::shared_ptr<SAMRAI::hier::PatchLevel> level(
        hierarchy->getPatchLevel(ln));
    //level->allocatePatchData(d_flag_inside,time);

    for (SAMRAI::hier::PatchLevel::iterator ip(level->begin()); ip != level->end(); ip++)
    {
        const std::shared_ptr<SAMRAI::hier::Patch> &patch = *ip;
        //A standard strategy that tags Gradient
        d_patch_strategy->tagGradientDetectorCells(*patch,
                                                   time,
                                                   initial_time,
                                                   tag_index,
                                                   uses_richardson_extrapolation_too);
        d_overset_strategy->applyOversetBoundaryDetector(*patch,
                                                         time,
                                                         initial_time,
                                                         tag_index);
    }
}

void PureGeometricIntegrator::findOversetCartBlocks(
    const std::shared_ptr<SAMRAI::hier::PatchHierarchy> &hierarchy,
    const int ln,
    const double time)
{
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT(hierarchy->getPatchLevel(ln));
    //std::cout << "level" << ln << " " << is_finest << '\n';

    std::shared_ptr<SAMRAI::hier::PatchLevel> level(
        hierarchy->getPatchLevel(ln));
    //level->allocatePatchData(d_flag_inside,time);

    for (SAMRAI::hier::PatchLevel::iterator ip(level->begin()); ip != level->end(); ip++)
    {
        const std::shared_ptr<SAMRAI::hier::Patch> &patch = *ip;
        //Call OversetStrategy for checking

        //If checking passes,then formCartBlock and allocate iblank for it.
        d_interp_strategy->formCartBlock(*patch, this, time);
    }
}

void PureGeometricIntegrator::extractCartGrid(
    int *nf, int *qstride, double **qnodein,
    int **idata, double **rdata,
    int *ngridsin, int *qnodesize,
    int ***ibl, double ***q,int** block_global_id, int *nblockin)
{
    d_interp_strategy->takeCartBlock(q, ibl, nblockin);
    *block_global_id = new int[*nblockin];
    
    d_interp_strategy->formCartGrid(
        nf, qstride, qnodein,
        ngridsin,qnodesize, idata, rdata, block_global_id);
}
