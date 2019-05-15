#ifndef GC_DIR_EL_H

#define GC_DIR_EL_H

#include "gc_critical.hpp"
#include "gc_random_walk.hpp"
#include "growth_names.hpp"

namespace growth
{
template <class ElongationModel, class DirectionModel>
class GrowthCone_Elongation_Direction : public DirectionModel,
                                        public ElongationModel

{

  public:
    GrowthCone_Elongation_Direction();

    GrowthCone_Elongation_Direction(
        const GrowthCone_Elongation_Direction<ElongationModel, DirectionModel>
            &);

    //~GrowthCone_Elongation_Direction();

    virtual GCPtr clone(BaseWeakNodePtr parent, NeuritePtr neurite,
                        double distanceToParent, std::string binaryID,
                        const Point &position, double angle) override;

    void compute_speed(mtPtr rnd_engine) override final;
    void prepare_for_split() override;
    void after_split() override;

    void set_status(const statusMap &) override;
    void get_status(statusMap &) const override;
};

// GrowthCone_Elongation_Direction<GrowthCone_Critical_Langevin>;
}
#endif
