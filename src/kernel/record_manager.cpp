#include "record_manager.hpp"

// include from elements
#include "GrowthCone.hpp"
#include "elements_types.hpp"

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
    size_t num_omp = kernel().parallelism_manager.get_num_local_threads();
    num_threads_   = num_omp;

    for (size_t i = 0; i < num_omp; i++)
    {
        omp_id_crec_.push_back(std::vector<size_t>());
    }
}


/**
 * @brief finalize the record file.
 * The name corresponds to the kernelID chosen in set_status
 */
void RecordManager::finalize()
{
    omp_id_crec_.clear();
    c_recorders_.clear();
    d_recorders_.clear();
    neuron_to_d_recorder_.clear();
}


/**
 * @brief Create a new recorder with custom parameters.
 *
 * @return The number of recorders created.
 */
size_t RecordManager::create_recorder(const std::vector<statusMap> &obj_params)
{
    // get new gid value from kernel
    size_t first_gid   = kernel().get_num_objects();
    size_t num_created = 0;

    size_t gid;
    std::string level, event_type;

    for (const auto &status : obj_params)
    {
        // get the targets
        std::vector<size_t> tgts;
        get_param(status, names::targets, tgts);

        // one recorder only takes care of neurons that are in the same OpenMP
        // thread as him, so we loop on the threads and create one recorder per
        // thread where there are target neurons
        for (int i = 0; i < num_threads_; i++)
        {
            gid = first_gid + num_created;

            // find the targets that are on the same thread
            std::vector<size_t> local_targets;
            int neuron_thread;

            for (size_t n : tgts)
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
            local_status[names::targets] = Property(local_targets);

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

                for (size_t n : local_targets)
                {
                    auto it = neuron_to_d_recorder_.find(n);
                    if (it == neuron_to_d_recorder_.end())
                    {
                        neuron_to_d_recorder_[n] = std::vector<size_t>();
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

                for (size_t n : local_targets)
                {
                    auto it = neuron_to_c_recorder_.find(n);
                    if (it == neuron_to_c_recorder_.end())
                    {
                        neuron_to_c_recorder_[n] = std::vector<size_t>();
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
    kernel().update_num_objects();

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
    for (size_t gid : omp_id_crec_[omp_id])
    {
        c_recorders_[gid]->record();
    }

    for (auto event : events_)
    {
        double resol = kernel().simulation_manager.get_resolution();

        size_t neuron = std::get<2>(event);
        for (size_t rec : omp_id_drec_[omp_id])
        {
            d_recorders_[rec]->record(event);
        }
    }
}


void RecordManager::finalize_simulation(size_t steps)
{
    size_t final_step = (steps > 0) ? steps - 1 : 0;

    for (auto &recorder : c_recorders_)
    {
        recorder.second->final_timestep(final_step);
    }
}


bool RecordManager::is_recorder(size_t gid) const
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


size_t RecordManager::num_recorders() const
{
    return c_recorders_.size() + d_recorders_.size();
}


void RecordManager::num_threads_changed(int num_omp)
{
    omp_id_crec_.clear();
    omp_id_drec_.clear();
    num_threads_ = num_omp;

    for (size_t i = 0; i < num_omp; i++)
    {
        omp_id_crec_.push_back(std::vector<size_t>());
        omp_id_drec_.push_back(std::vector<size_t>());
    }
}


void RecordManager::new_branching_event(const Event &ev)
{
    if (!d_recorders_.empty())
    {
        size_t neuron_gid = std::get<2>(ev);
        auto it           = neuron_to_d_recorder_.find(neuron_gid);

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
        size_t neuron_gid = std::get<2>(ev);
        auto it           = neuron_to_c_recorder_.find(neuron_gid);

        if (it != neuron_to_c_recorder_.end())
        {
            for (auto rec_gid : it->second)
            {
                c_recorders_[rec_gid]->record(ev);
            }
        }
    }
}


statusMap RecordManager::get_recorder_status(size_t gid) const
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
    set_param(status, names::event_type, event_type);

    unsigned int level = rec->get_level();
    std::string s;

    switch (level)
    {
    case 0:
        s = "neuron";
        set_param(status, names::level, s);
        break;
    case 1:
        s = "neurite";
        set_param(status, names::level, s);
        break;
    case 2:
        s = "growth_cone";
        set_param(status, names::level, s);
        break;
    default:
        throw InvalidParameter("Invalid level '" + std::to_string(level) + "'.",
                               __FUNCTION__, __FILE__, __LINE__);
        break;
    }

    return status;
}


void RecordManager::set_status(const statusMap &status) {}


void RecordManager::get_recorder_type(size_t gid, std::string &level,
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


bool RecordManager::get_next_recording(size_t gid, std::vector<Property> &ids,
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


bool RecordManager::get_next_time(size_t gid, std::vector<Property> &ids,
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


} /* namespace */
