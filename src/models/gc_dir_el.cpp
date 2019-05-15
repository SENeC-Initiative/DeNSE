#include "gc_dir_el.hpp"
#include "config_impl.hpp"
#include <cassert>
#include <math.h>
#include <memory>


namespace growth
{
template <class ElongationModel, class DirectionModel>
GrowthCone_Elongation_Direction<
    ElongationModel, DirectionModel>::GrowthCone_Elongation_Direction()
    : GrowthCone()
    , DirectionModel()
    , ElongationModel()

{
}


template <class ElongationModel, class DirectionModel>
GrowthCone_Elongation_Direction<ElongationModel, DirectionModel>::
    GrowthCone_Elongation_Direction(
        const GrowthCone_Elongation_Direction<ElongationModel, DirectionModel>
            &copy)
    : //@TODO!!
      // here the right copy constructor to call is
      // ElongationModel(copy)
      // btw it's not working properly: it won't break, but neither will execute
      // the  lines into the ElongationModel copy-constructor  I think the issue is
      // related to this:
      // https://stackoverflow.com/questions/19167201/copy-constructor-of-template-class
      // but I cannot find a proper solution.
    GrowthCone(copy)
    , GrowthCone_RandomWalk(copy)
    , GrowthCone_Critical(copy)

{
}

template <class ElongationModel, class DirectionModel>
GCPtr GrowthCone_Elongation_Direction<ElongationModel, DirectionModel>::clone(
    BaseWeakNodePtr parent, NeuritePtr neurite, double distanceToParent,
    std::string binaryID, const Point &position, double angle)
{
#ifndef NDEBUG
    printf(" It's calling RandomWalk_Critical->clone! with direction %f\n",
           angle);
#endif
    auto newCone = std::make_shared<
        GrowthCone_Elongation_Direction<ElongationModel, DirectionModel>>(
        *this);
    // newCone= std::dynamic_pointer_cast<GrowthCone>(newCone);
    newCone->update_topology(parent, neurite, distanceToParent, binaryID,
                             position, angle);
    // assert(newCone->ElongationModel::get_CR_demand() ==
    // ElongationModel::get_CR_demand());
    // printf("newcone  %s critical_resource is %f \n", get_treeID().c_str(),
    // get_critical_resource_demand());
    // printf("newcone  %s critical_resource is %f \n",
    // newCone->get_treeID().c_str(),
    // newCone->get_critical_resource_demand());
    return newCone;
}


template <class ElongationModel, class DirectionModel>
void GrowthCone_Elongation_Direction<
    ElongationModel, DirectionModel>::set_status(const statusMap &status)
{
    ElongationModel::set_status(status);
    DirectionModel::average_speed_ = ElongationModel::critical_.speed_factor;
    DirectionModel::set_status(status);
}


template <class ElongationModel, class DirectionModel>
void GrowthCone_Elongation_Direction<ElongationModel,
                                     DirectionModel>::prepare_for_split()
{
    DirectionModel::prepare_for_split();
    ElongationModel::prepare_for_split();
}


template <class ElongationModel, class DirectionModel>
void GrowthCone_Elongation_Direction<ElongationModel,
                                     DirectionModel>::after_split()
{
    DirectionModel::after_split();
    ElongationModel::after_split();
}


template <class ElongationModel, class DirectionModel>
void GrowthCone_Elongation_Direction<
    ElongationModel, DirectionModel>::get_status(statusMap &status) const
{
    DirectionModel::get_status(status);
    ElongationModel::get_status(status);
}


template <class ElongationModel, class DirectionModel>
void GrowthCone_Elongation_Direction<
    ElongationModel, DirectionModel>::compute_speed(mtPtr rnd_engine)
{
    ElongationModel::compute_speed(rnd_engine);
}

// getter functions here
}
