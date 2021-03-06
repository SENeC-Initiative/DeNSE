/*
 * extension_gfluct.cpp
 *
 * This file is part of DeNSE.
 *
 * Copyright (C) 2019 SeNEC Initiative
 *
 * DeNSE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * DeNSE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with DeNSE. If not, see <http://www.gnu.org/licenses/>.
 */

#include "extension_gfluct.hpp"

#include "config.hpp"


namespace growth
{

GFluctExtensionModel::GFluctExtensionModel(GCPtr gc, NeuritePtr neurite)
    : ExtensionModel(gc, neurite)
    , speed_gc_avg_(SPEED_GROWTH_CONE)
    , local_speed_(SPEED_GROWTH_CONE)
    , speed_gc_std_(SPEED_GROWTH_CONE) // std = mean
    , local_std_(SPEED_GROWTH_CONE)
{
    normal_ = std::normal_distribution<double>(0, 1);
}


GFluctExtensionModel::GFluctExtensionModel(const GFluctExtensionModel &copy,
                                           GCPtr gc, NeuritePtr neurite)
    : ExtensionModel(copy, gc, neurite)
    , speed_gc_avg_(copy.speed_gc_avg_)
    , local_speed_(copy.local_speed_)
    , speed_gc_std_(copy.speed_gc_std_)
    , local_std_(copy.local_std_)
{
    normal_ = std::normal_distribution<double>(0, 1);
}


double GFluctExtensionModel::compute_speed(mtPtr rnd_engine, double substep)
{
    return local_speed_ +
           local_std_ * sqrt(substep) * normal_(*(rnd_engine).get());
}


/**
 * @brief Update growth cone average speed.
 * This function updates the average speed, e.g. because the number
 * of growth cones sustained by the neurite changed.
 * Since the average speed changed, the local speed must also be
 * updated.
 */
void GFluctExtensionModel::update_speed(double speed_factor)
{
    double area_factor = local_speed_ / speed_gc_avg_;

    speed_gc_avg_ *= speed_factor;
    local_speed_   = speed_gc_avg_ * area_factor;
}


void GFluctExtensionModel::update_local_speed(double area_factor)
{
    local_speed_ = speed_gc_avg_ * area_factor;
    local_std_   = speed_gc_std_ / area_factor;
}


double GFluctExtensionModel::get_max_speed() const
{
    return local_speed_ + 5. * local_std_;
}


bool GFluctExtensionModel::set_status(const statusMap &status)
{
    double std;
    bool bspeed, bstd;

    bspeed = get_param(status, names::speed_growth_cone, speed_gc_avg_);

    bstd = get_param(status, names::speed_variance, std);

    if (bstd)
    {
        if (std < 0)
        {
            throw std::invalid_argument("`" + names::speed_variance +
                                        "` must be positive.");
        }

        speed_gc_std_ = std;
    }

    return bspeed + bstd;
}


void GFluctExtensionModel::get_status(statusMap &status) const
{
    set_param(status, names::speed_growth_cone, speed_gc_avg_,
              "micrometer/minute");
    set_param(status, names::speed_variance, speed_gc_std_,
              "micrometer/minute");
}

} // namespace growth
