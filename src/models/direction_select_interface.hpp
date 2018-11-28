#ifndef DIRSELMODEL_H
#define DIRSELMODEL_H

#include <vector>
#include <string>

// lib includes
#include "cttrie.hpp"
#include "elements_types.hpp"

// elements include
#include "GrowthCone.hpp"
#include "Neurite.hpp"


namespace growth
{

/**
 * @brief Direction selection interface, from which all submodels must inherit.
 *
 * This class serves as a guide for developers wanting to create a new
 * direction selection model:
 *
 * - Default (copy) constructors must not be used and are therefore deleted.
 * - Correct constructor must get a valid instance of
 *   :cpp:class:`growth::GrowthCone` and :cpp:class:`growth::Neurite`
 *   shared pointer to store in `gc_weakptr_` and `neurite_ptr_`.
 * - If necessary, the correct constructor should update the `observables_` and
 *   override the :cpp:func:`growth::DirectionSelectModel::get_state` function.
 * - There is only one method that has to be implemented for a direction
 *   model to be valid, which is
 *   :cpp:func:`growth::DirectionSelectModel::compute_target_angle`.
 *
 * See also
 * --------
 * :cpp:class:`growth::NMDirectionSelector` or
 * :cpp:class:`growth::RTDirectionSelector` for complete model
 * declarations.
 */
class DirectionSelectModel;
typedef std::shared_ptr<DirectionSelectModel> DSPtr;


class DirectionSelectModel
{
  protected:
    // GrowthCone must be destructed upon release by Neurite (don't keep it
    // alive) and we don't need to modify it (done through template)
    std::weak_ptr<GrowthCone> gc_weakptr_;
    // We must be able to modify Neurite (branching)
    std::shared_ptr<Neurite> neurite_ptr_;

    std::vector<std::string> observables_;

  public:
    DirectionSelectModel() = delete;
    DirectionSelectModel(const DirectionSelectModel& copy) = delete;

    DirectionSelectModel(GCPtr gc, NeuritePtr neurite)
      : gc_weakptr_(GCPtr(gc)), neurite_ptr_(neurite) {};

    DirectionSelectModel(const DirectionSelectModel& copy, GCPtr gc, NeuritePtr neurite)
      : gc_weakptr_(gc), neurite_ptr_(neurite), observables_(copy.observables_) {};

    
    virtual void
    compute_target_angle(
        const std::vector<double> &directions_weights, const Filopodia &filo,
        mtPtr rnd_engine, double total_proba, bool interacting,
        double old_angle, double &substep, double &step_length,
        double &new_angle, bool &stopped) = 0;

    void get_observables(std::vector<std::string> &obs) const
    {
        obs.insert(obs.end(), observables_.begin(), observables_.end());
    }

    virtual void prepare_for_split() {};
    virtual void after_split() {};
    virtual void kernel_updated() {};
    virtual double get_state(const char *observable) const { return std::nan(""); };
    virtual void set_status(const statusMap &status) = 0;
    virtual void get_status(statusMap &status) const = 0;
};

} // namespace growth

#endif
