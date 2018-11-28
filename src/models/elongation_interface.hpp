#ifndef ELMODEL_H
#define ELMODEL_H

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

class ElongationModel;
typedef std::shared_ptr<ElongationModel> EMPtr;


/**
 * @brief Elongation interface, from which all elongation models must inherit.
 *
 * This class serves as a guide for developers wanting to create a new
 * elongation model:
 *
 * - Default (copy) constructors must not be used and are therefore deleted.
 * - Correct constructor must get a valid instance of
 *   :cpp:class:`growth::GrowthCone` and :cpp:class:`growth::Neurite`
 *   shared pointer to store in `gc_weakptr_` and `neurite_ptr_`.
 * - If necessary, the correct constructor should update the `observables_` and
 *   override the :cpp:func:`growth::ElongationModel::get_state` function.
 * - There is only one method that has to be implemented for an elongation
 *   model to be valid, which is
 *   :cpp:func:`growth::ElongationModel::compute_speed`.
 *
 * See also
 * --------
 * :cpp:class:`growth::CstElongationModel` or
 * :cpp:class:`growth::ResourceBasedElongationModel` for complete model
 * declarations.
 */
class ElongationModel
{
  protected:
    // GrowthCone must be destructed upon release by Neurite (don't keep it
    // alive) and we don't need to modify it (done through template)
    std::weak_ptr<GrowthCone> gc_weakptr_;
    // We must be able to modify Neurite (branching)
    std::shared_ptr<Neurite> neurite_ptr_;
    std::vector<std::string> observables_;

  public:
    ElongationModel() = delete;
    ElongationModel(const ElongationModel& copy) = delete;

    ElongationModel(GCPtr gc, NeuritePtr neurite)
      : gc_weakptr_(GCPtr(gc)), neurite_ptr_(neurite) {};

    ElongationModel(const ElongationModel& copy, GCPtr gc, NeuritePtr neurite)
      : gc_weakptr_(gc), neurite_ptr_(neurite), observables_(copy.observables_) {};

    virtual double compute_speed(mtPtr rnd_engine, double substep) = 0;

    void get_observables(std::vector<std::string> &obs) const
    {
        obs.insert(obs.end(), observables_.begin(), observables_.end());
    }

    virtual void prepare_for_split() {};
    virtual void after_split() {};
    virtual double get_state(const char *observable) const { return std::nan(""); };
    virtual void set_status(const statusMap &status) = 0;
    virtual void get_status(statusMap &status) const = 0;
};

} // namespace growth

#endif
