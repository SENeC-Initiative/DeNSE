/*
 * record_manager.cpp
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

#include "record_manager.hpp"

// include from elements
#include "GrowthCone.hpp"
#include "elements_types.hpp"

// include from lib
#include "config.hpp"

// include from kernel
#include "kernel_manager.hpp"


namespace growth
{

/**
 * @brief RecordManager class
 * The RecordManager saves relevant information on the dynamics of the
 * simulation. The information are stored during the execution, it's called from
 * the SimulationManager, or from the Neurite during branch events. The desired
 * information can be set throw the Kernel_set_status method.
 */
RecordManager::RecordManager()
    : last_omp_id_(0)
{
}


/**
 * @brief initialize the record file.
 * The name corresponds to the kernelID chosen in set_status
 */
void RecordManager::initialize()
{
    stype num_omp = kernel().parallelism_manager.get_num_local_threads();
    num_threads_  = num_omp;

    for (stype i = 0; i < num_omp; i++)
    {
        omp_id_crec_.push_back(std::vector<stype>());
    }
}


/**
 * @brief finalize the record file.
 * The name corresponds to the kernelID chosen in set_status
 */
void RecordManager::finalize()
{
    omp_id_crec_.clear();
    omp_id_drec_.clear();
    c_recorders_.clear();
    d_recorders_.clear();
    neuron_to_d_recorder_.clear();
    neuron_to_c_recorder_.clear();
    events_.clear();
}


/**
 * @brief Create a new recorder with custom parameters.
 *
 * @return The number of recorders created.
 */
stype RecordManager::create_recorder(const std::vector<statusMap> &obj_params)
{
    // get new gid value from kernel
    stype first_gid   = kernel().get_num_created_objects();
    stype num_created = 0;
    stype gid;

    std::string level, event_type;

    for (const auto &status : obj_params)
    {
        // get the targets
        std::vector<stype> tgts;
        get_param(status, names::targets, tgts);

        // one recorder only takes care of neurons that are in the same OpenMP
        // thread as him, so we loop on the threads and create one recorder per
        // thread where there are target neurons
        for (int i = 0; i < num_threads_; i++)
        {
            gid = first_gid + num_created;

            // find the targets that are on the same thread
            std::vector<stype> local_targets;
            int neuron_thread;

            for (stype n : tgts)
            {
                neuron_thread = kernel().neuron_manager.get_neuron_thread(n);
                if (neuron_thread == i)
                {
                    local_targets.push_back(n);
                }
            }

            // create a local statusMap
            statusMap local_status;
            for (auto &status_it : status)
            {
                if (status_it.first != names::targets)
                {
                    local_status[status_it.first] = Property(status_it.second);
                }
            }
            local_status[names::targets] = Property(local_targets, "");

            // get recorder type
            get_param(local_status, names::level, level);
            get_param(local_status, names::event_type, event_type);

            // create recorder and set its status
            if (level == "neuron")
            {
                if (event_type == "continuous")
                {
                    c_recorders_[gid] =
                        std::make_shared<NeuronContinuousRecorder>();
                }
                else if (event_type == "discrete")
                {
                    d_recorders_[gid] =
                        std::make_shared<NeuronDiscreteRecorder>();
                }
                else
                {
                    throw BadPropertyType(
                        "event_type", "continuous' or 'discrete", event_type,
                        __FUNCTION__, __FILE__, __LINE__);
                }
            }
            else if (level == "neurite")
            {
                if (event_type == "continuous")
                {
                    c_recorders_[gid] =
                        std::make_shared<NeuriteContinuousRecorder>();
                }
                else if (event_type == "discrete")
                {
                    d_recorders_[gid] =
                        std::make_shared<NeuriteDiscreteRecorder>();
                }
                else
                {
                    throw BadPropertyType(
                        "event_type", "continuous' or 'discrete", event_type,
                        __FUNCTION__, __FILE__, __LINE__);
                }
            }
            else if (level == "growth_cone")
            {
                if (event_type == "continuous")
                {
                    c_recorders_[gid] =
                        std::make_shared<GrowthConeContinuousRecorder>();
                }
                else if (event_type == "discrete")
                {
                    d_recorders_[gid] =
                        std::make_shared<GrowthConeDiscreteRecorder>();
                }
                else
                {
                    throw BadPropertyType(
                        "event_type", "continuous' or 'discrete", event_type,
                        __FUNCTION__, __FILE__, __LINE__);
                }
            }
            else
            {
                throw BadPropertyType("level",
                                      "neuron', 'neurite' or 'growth_cone",
                                      level, __FUNCTION__, __FILE__, __LINE__);
            }

            // set status and, if discrete, map targets to the recorder
            if (event_type == "discrete")
            {
                // set recorder status
                d_recorders_[gid]->set_status(local_status);

                for (stype n : local_targets)
                {
                    auto it = neuron_to_d_recorder_.find(n);
                    if (it == neuron_to_d_recorder_.end())
                    {
                        neuron_to_d_recorder_[n] = std::vector<stype>();
                    }
                    neuron_to_d_recorder_[n].push_back(gid);
                }

                // assign recorder to a specific thread
                omp_id_drec_[last_omp_id_].push_back(gid);
            }
            else
            {
                // set recorder status
                c_recorders_[gid]->set_status(local_status);

                for (stype n : local_targets)
                {
                    auto it = neuron_to_c_recorder_.find(n);
                    if (it == neuron_to_c_recorder_.end())
                    {
                        neuron_to_c_recorder_[n] = std::vector<stype>();
                    }
                    neuron_to_c_recorder_[n].push_back(gid);
                }

                // assign recorder to a specific thread
                omp_id_crec_[last_omp_id_].push_back(gid);
            }

            // update last_omp_id_
            last_omp_id_++;

            if (last_omp_id_ == num_threads_)
            {
                last_omp_id_ = 0;
            }

            // update number of objects
            num_created++;
        }
    }

    // update kernel count
    kernel().update_num_objects(num_created);

    return num_created;
}


/**
 * @brief Make all records for the current step
 *
 * Loop over the recorders and make them record their observables.
 *
 * @warning this function must only be called inside the OpenMP loop of the
 * SimulationManager::simulate funtion!
 */
void RecordManager::record(int omp_id)
{
    for (stype gid : omp_id_crec_[omp_id])
    {
        c_recorders_[gid]->record();
    }

    for (auto event : events_)
    {
        stype neuron = std::get<edata::NEURON>(event);
        for (stype rec : omp_id_drec_[omp_id])
        {
            d_recorders_[rec]->record(event);
        }
    }
}


void RecordManager::finalize_simulation(stype steps)
{
    stype final_step = (steps > 0) ? steps - 1 : 0;

    for (auto &recorder : c_recorders_)
    {
        recorder.second->final_timestep(final_step);
    }
}


bool RecordManager::is_recorder(stype gid) const
{
    auto c_it = c_recorders_.find(gid);
    auto d_it = d_recorders_.find(gid);

    if (c_it != c_recorders_.end())
    {
        return true;
    }
    else if (d_it != d_recorders_.end())
    {
        return true;
    }

    return false;
}


void RecordManager::get_status(statusMap &status) const {}


void RecordManager::get_defaults(statusMap &status) const
{
    auto rec = std::make_shared<NeuronContinuousRecorder>();
    rec->get_status(status);
}


stype RecordManager::num_recorders() const
{
    return c_recorders_.size() + d_recorders_.size();
}


void RecordManager::num_threads_changed(int num_omp)
{
    omp_id_crec_.clear();
    omp_id_drec_.clear();
    num_threads_ = num_omp;

    for (stype i = 0; i < num_omp; i++)
    {
        omp_id_crec_.push_back(std::vector<stype>());
        omp_id_drec_.push_back(std::vector<stype>());
    }
}


void RecordManager::new_branching_event(const Event &ev)
{
    if (!d_recorders_.empty())
    {
        stype neuron_gid = std::get<edata::NEURON>(ev);
        auto it          = neuron_to_d_recorder_.find(neuron_gid);

        if (it != neuron_to_d_recorder_.end())
        {
            for (auto rec_gid : it->second)
            {
                d_recorders_[rec_gid]->record(ev);
            }
        }
    }

    if (!c_recorders_.empty())
    {
        stype neuron_gid = std::get<edata::NEURON>(ev);
        auto it          = neuron_to_c_recorder_.find(neuron_gid);

        if (it != neuron_to_c_recorder_.end())
        {
            for (auto rec_gid : it->second)
            {
                c_recorders_[rec_gid]->record(ev);
            }
        }
    }
}


void RecordManager::neurons_deleted(const std::vector<stype> &gids)
{
    for (stype neuron : gids)
    {
        auto v_crec = neuron_to_c_recorder_.find(neuron);

        if (v_crec != neuron_to_c_recorder_.end())
        {
            for (stype rec_id : v_crec->second)
            {
                c_recorders_[rec_id]->neuron_deleted(neuron);
                neuron_to_c_recorder_.erase(neuron);
            }
        }

        auto v_drec = neuron_to_d_recorder_.find(neuron);

        if (v_drec != neuron_to_d_recorder_.end())
        {
            for (stype rec_id : v_drec->second)
            {
                d_recorders_[rec_id]->neuron_deleted(neuron);
                neuron_to_d_recorder_.erase(neuron);
            }
        }
    }
}


void RecordManager::new_neurite(stype neuron, const std::string &neurite)
{
    auto v_crec = neuron_to_c_recorder_.find(neuron);

    if (v_crec != neuron_to_c_recorder_.end())
    {
        for (stype rec_id : v_crec->second)
        {
            c_recorders_[rec_id]->new_neurite(neuron, neurite);
        }
    }

    auto v_drec = neuron_to_d_recorder_.find(neuron);

    if (v_drec != neuron_to_d_recorder_.end())
    {
        for (stype rec_id : v_drec->second)
        {
            d_recorders_[rec_id]->new_neurite(neuron, neurite);
        }
    }
}


void RecordManager::gc_died(stype neuron, const std::string &neurite,
                            stype gc_id)
{
    // this information is only relevant for continuous recorders to know what
    // is the last time for the gc records.
    // for discrete recorders there is no issue since no new event will arrive.
    auto v_crec = neuron_to_c_recorder_.find(neuron);

    if (v_crec != neuron_to_c_recorder_.end())
    {
        for (stype rec_id : v_crec->second)
        {
            c_recorders_[rec_id]->gc_died(neuron, neurite, gc_id);
        }
    }
}


statusMap RecordManager::get_recorder_status(stype gid) const
{
    statusMap status;
    std::shared_ptr<BaseRecorder> rec;
    std::string event_type;

    auto c_it = c_recorders_.find(gid);
    auto d_it = d_recorders_.find(gid);

    if (c_it != c_recorders_.end())
    {
        rec        = c_recorders_.at(gid);
        event_type = "continuous";
    }
    else if (d_it != d_recorders_.end())
    {
        rec        = d_recorders_.at(gid);
        event_type = "discrete";
    }
    else
    {
        throw InvalidArg("Gid " + std::to_string(gid) + " is not a recorder.",
                         __FUNCTION__, __FILE__, __LINE__);
    }

    rec->get_status(status);

    // set status with level and event type
    set_param(status, names::event_type, event_type, "");

    unsigned int level = rec->get_level();
    std::string s;

    switch (level)
    {
    case 0:
        s = "neuron";
        set_param(status, names::level, s, "");
        break;
    case 1:
        s = "neurite";
        set_param(status, names::level, s, "");
        break;
    case 2:
        s = "growth_cone";
        set_param(status, names::level, s, "");
        break;
    default:
        throw InvalidParameter("Invalid level '" + std::to_string(level) + "'.",
                               __FUNCTION__, __FILE__, __LINE__);
        break;
    }

    return status;
}


void RecordManager::set_status(const statusMap &status) {}


void RecordManager::get_recorder_type(stype gid, std::string &level,
                                      std::string &event_type) const
{
    if (is_recorder(gid))
    {
        auto c_it = c_recorders_.find(gid);
        auto d_it = d_recorders_.find(gid);

        if (c_it != c_recorders_.end())
        {
            level      = c_recorders_.at(gid)->get_level();
            event_type = "continuous";
        }
        else if (d_it != d_recorders_.end())
        {
            level      = d_recorders_.at(gid)->get_level();
            event_type = "discrete";
        }
    }
    else
    {
        throw std::runtime_error("Object is not a recorder.");
    }
}


bool RecordManager::get_next_recording(stype gid, std::vector<Property> &ids,
                                       std::vector<double> &values)
{
    auto c_it = c_recorders_.find(gid);
    auto d_it = d_recorders_.find(gid);

    if (c_it != c_recorders_.end())
    {
        return c_recorders_[gid]->get_next_recording(ids, values);
    }
    else if (d_it != d_recorders_.end())
    {
        return d_recorders_[gid]->get_next_recording(ids, values);
    }
    else
    {
        throw InvalidArg("Gid " + std::to_string(gid) + " is not a recorder.",
                         __FUNCTION__, __FILE__, __LINE__);
    }

    return false;
}


bool RecordManager::get_next_time(stype gid, std::vector<Property> &ids,
                                  std::vector<double> &values,
                                  const std::string &time_units)
{
    auto c_it = c_recorders_.find(gid);
    auto d_it = d_recorders_.find(gid);

    if (c_it != c_recorders_.end())
    {
        return c_recorders_[gid]->get_next_time(ids, values, time_units);
    }
    else if (d_it != d_recorders_.end())
    {
        return d_recorders_[gid]->get_next_time(ids, values, time_units);
    }
    else
    {
        throw InvalidArg("Gid " + std::to_string(gid) + " is not a recorder.",
                         __FUNCTION__, __FILE__, __LINE__);
    }

    return false;
}


} // namespace growth
