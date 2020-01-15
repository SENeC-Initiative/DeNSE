/*
 * extension_gfluct.hpp
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

#ifndef GFLUCT_EL_H
#define GFLUCT_EL_H

// elements
#include "extension_interface.hpp"


namespace growth
{

class GFluctExtensionModel : public virtual ExtensionModel
{
  protected:
    double speed_gc_avg_;
    double speed_gc_std_;

    std::normal_distribution<double> normal_;

  public:
    GFluctExtensionModel(GCPtr gc, NeuritePtr neurite);
    GFluctExtensionModel(const GFluctExtensionModel &copy) = delete;
    GFluctExtensionModel(const GFluctExtensionModel &copy, GCPtr gc,
                         NeuritePtr neurite);

    double compute_speed(mtPtr rnd_engine, double substep) override final;

    virtual void set_status(const statusMap &status) override final;

    virtual void get_status(statusMap &status) const override final;
};

} // namespace growth

#endif
