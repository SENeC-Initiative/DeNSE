#include "extension_gfluct.hpp"

#include "config.hpp"


namespace growth
{

GFluctExtensionModel::GFluctExtensionModel(GCPtr gc, NeuritePtr neurite)
  : ExtensionModel(gc, neurite)
  , speed_gc_avg_(SPEED_GROWTH_CONE)
  , speed_gc_std_(SPEED_GROWTH_CONE)  // std = mean
{
    normal_ = std::normal_distribution<double>(0, 1);
}


GFluctExtensionModel::GFluctExtensionModel(const GFluctExtensionModel &copy, GCPtr gc, NeuritePtr neurite)
  : ExtensionModel(copy, gc, neurite)
  , speed_gc_avg_(copy.speed_gc_avg_)
  , speed_gc_std_(copy.speed_gc_std_)
{
    normal_ = std::normal_distribution<double>(0, 1);
}


double GFluctExtensionModel::compute_speed(mtPtr rnd_engine, double substep)
{
    return speed_gc_avg_ +
           speed_gc_std_*sqrt(substep)*normal_(*(rnd_engine).get());
}


void GFluctExtensionModel::set_status(const statusMap &status)
{
    double std;
    bool b;

    get_param(status, names::speed_growth_cone, speed_gc_avg_);

    b = get_param(status, names::speed_variance, std);
    if (b)
    {
        if (std < 0)
        {
            throw std::invalid_argument(
                "`" + names::speed_variance + "` must be positive.");
        }

        speed_gc_std_ = std;
    }
}


void GFluctExtensionModel::get_status(statusMap &status) const
{
    set_param(status, names::speed_growth_cone, speed_gc_avg_, "micrometer/minute");
    set_param(status, names::speed_variance, speed_gc_std_, "micrometer/minute");
}

}
