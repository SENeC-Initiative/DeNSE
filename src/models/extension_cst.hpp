/*
 * extension_cst.hpp
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

#ifndef CST_EL_H
#define CST_EL_H

// elements
#include "extension_interface.hpp"

// libgrowth
#include "config_impl.hpp"


namespace growth
{

class CstExtensionModel : public virtual ExtensionModel
{
  protected:
    double speed_growth_cone_;
    double local_avg_speed_;

  public:
    CstExtensionModel(GCPtr gc, NeuritePtr neurite)
        : ExtensionModel(gc, neurite)
        , speed_growth_cone_(SPEED_GROWTH_CONE)
        , local_avg_speed_(SPEED_GROWTH_CONE){};

    CstExtensionModel(const CstExtensionModel &copy) = delete;

    CstExtensionModel(const CstExtensionModel &copy, GCPtr gc,
                      NeuritePtr neurite)
        : ExtensionModel(copy, gc, neurite)
        , speed_growth_cone_(copy.speed_growth_cone_)
        , local_avg_speed_(copy.local_avg_speed_){};

    double compute_speed(mtPtr rnd_engine, double substep) override final
    {
        return local_avg_speed_;
    };

    bool set_status(const statusMap &status) override final
    {
        return get_param(status, names::speed_growth_cone,
                         speed_growth_cone_);
    };

    void get_status(statusMap &status) const override final
    {
        set_param(status, names::speed_growth_cone, speed_growth_cone_,
                  "micrometer / minute");
    };

    void update_speed(double speed_factor) override final
    {
        double area_factor  = local_avg_speed_ / speed_growth_cone_;

        speed_growth_cone_ *= speed_factor;
        local_avg_speed_    = speed_growth_cone_ * area_factor;
    };

    void update_local_speed(double area_factor) override final
    {
        local_avg_speed_ = speed_growth_cone_ * area_factor;
    };

    double get_max_speed() const override final
    {
        return local_avg_speed_;
    };
};

} // namespace growth

#endif
