#include "recorders.hpp"

// element include
#include "Neuron.hpp"
#include "Neurite.hpp"
#include "GrowthCone.hpp"

// libgrowth include
#include "config_impl.hpp"
#include "growth_names.hpp"

// kernel includes
#include "kernel_manager.hpp"


namespace growth
{

/**
 * Constructor for BaseRecorder
 */
BaseRecorder::BaseRecorder()
  : record_to_("memory")
  , record_file_()
  , interval_(-1)
  , restrict_to_("")
  , v_iterating_(false)
  , t_iterating_(false)
  , cstr_obs_("")
{
}


void BaseRecorder::get_status(statusMap &status) const
{
    set_param(status, names::record_to, record_to_);
    set_param(status, names::interval, interval_);
    set_param(status, names::restrict_to, restrict_to_);
    set_param(status, names::observable, observable_);

    // get the neuron ids
    std::vector<size_t> gids;
    for (auto const& neuron : targets_)
    {
        gids.push_back(neuron.first);
    }
    set_param(status, names::targets, gids);

    // @todo: find a way to bring the data back to Python
}


void BaseRecorder::set_status(const statusMap &status)
{
    get_param(status, names::record_to, record_to_);
    get_param(status, names::interval, interval_);

    // for the last parameters, check that simulation was not started
    // or that nothing was previously recorded
    if (targets_.empty() || kernel().simulation_manager.get_current_step() == 0)
    {
        get_param(status, names::restrict_to, restrict_to_);
        get_param(status, names::observable, observable_);
        cstr_obs_ = observable_.c_str();
    }

    // @todo if record_to_ changed, do the necessary updates
}


void BaseRecorder::record() {}


/**
 * @brief Set the final times for the continuous recorders
 */
void BaseRecorder::final_timestep(size_t step) {}


unsigned int BaseRecorder::get_level() const {
    return 3;
}


unsigned int BaseRecorder::get_event_type() const {
    return 2;
}


bool BaseRecorder::get_next_recording(std::vector<Property>& ids,
                                      std::vector<double>& values)
{
    return false;
}


bool BaseRecorder::get_next_time(std::vector<Property>& ids,
                                 std::vector<double>& values)
{
    return false;
}


void BaseRecorder::reset_iterations()
{
    v_iterating_ = false;
    t_iterating_ = false;
}


/**
 * Constructor for NeuronContinuousRecorder
 */
NeuronContinuousRecorder::NeuronContinuousRecorder()
{
    times_ = std::array<double, 3>({0., 0., 0.});
}


/**
 * @brief Set the status of "targets"
 *
 * Call default BaseRecorder::set_status, then do a custom update for targets
 */
void NeuronContinuousRecorder::set_status(const statusMap &status)
{
    // call default parent set_status
    BaseRecorder::set_status(status);

    // set targets and recording_
    if (targets_.empty() || kernel().simulation_manager.get_current_step() == 0)
    {
        // get the targets
        std::vector<size_t> gids;
        bool targets_changed = get_param(status, names::targets, gids);

        if (targets_changed)
        {
            recording_.clear();
            times_ = std::array<double, 3>({0., 0., 0.});
        }

        // initialize the recording container
        for (size_t gid : gids)
        {
            NeuronPtr n = kernel().neuron_manager.get_neuron(gid);
            targets_.insert({gid, n});

            recording_.insert({gid, std::vector<double>()});
        }
    }
}


void NeuronContinuousRecorder::record()
{
    if (record_to_ == "memory")
    {
        for (const auto& neuron : targets_)
        {
            recording_[neuron.first].push_back(
                neuron.second->get_state(cstr_obs_));
        }
        times_[2] += 1;
    }
    else
    {
        // @todo: see with Alessio
    }
}


/**
 * @brief Set the final times for the continuous records
 */
void NeuronContinuousRecorder::final_timestep(size_t step)
{
    times_[1] = step;
}


/**
 * @brief Return the level at which the recorder works
 *
 * Returns 0 if recording at neuron level, 1 if recording at neurite level and
 * 2 for growth cone level.
 * This function returns 0.
 */
unsigned int NeuronContinuousRecorder::get_level() const
{
    return 0;
}


/**
 * @brief Return the type of event recorded by the recorder
 *
 * Returns 0/1 if continuous/discrete events are recorded.
 * This function returns 0.
 */
unsigned int NeuronContinuousRecorder::get_event_type() const
{
    return 0;
}


/**
 * @brief Access the data of the next object recorded, return false when done.
 *
 * Data about the next object is given through a list properties (gid of the
 * neuron, name of the neurite, id of the growth cone) and a list of values.
 * The function returns true if there is still data left, otherwise it returns
 * false. Once all data has been looped on, the iterators are reset to their
 * initial values.
 */
bool NeuronContinuousRecorder::get_next_recording(std::vector<Property>& ids,
                                                  std::vector<double>& values)
{
    if (!v_iterating_)
    {
        v_iterating_ = true;
        neuron_it_ = recording_.cbegin();
    }

    if (neuron_it_ != recording_.cend())
    {
        Property p = Property(neuron_it_->first);
        ids.push_back(p);
        values.reserve(values.size() + neuron_it_->second.size());
        values.insert(
            values.end(), neuron_it_->second.cbegin(),
            neuron_it_->second.cend());

        neuron_it_++;

        if (neuron_it_ != recording_.cend())
        {
            return true;
        }
        else
        {
            v_iterating_ = false;
            return false;
        }
    }
    else
    {
        v_iterating_ = false;
        return false;
    }
}


/**
 * @brief Access the times for the next object recorded, return false when done.
 *
 * Times about the next object is given through a list properties (gid of the
 * neuron, name of the neurite, id of the growth cone) and a list of values.
 * The function returns true if there is still data left, otherwise it returns
 * false. Once all data has been looped on, the iterators are reset to their
 * initial values.
 */
bool NeuronContinuousRecorder::get_next_time(std::vector<Property>& ids,
                                             std::vector<double>& values)
{
    values.insert(values.end(), times_.begin(), times_.end());
    return false;
}


/**
 * Constructor for NeuronDiscreteRecorder
 */
NeuronDiscreteRecorder::NeuronDiscreteRecorder()
{}


/**
 * @brief Set the status of "targets"
 *
 * Call default BaseRecorder::set_status, then do a custom update for targets
 */
void NeuronDiscreteRecorder::set_status(const statusMap &status)
{
    // call default parent set_status
    BaseRecorder::set_status(status);

    // set targets and recording_
    if (targets_.empty() || kernel().simulation_manager.get_current_step() == 0)
    {
        // get the targets
        std::vector<size_t> gids;
        bool targets_changed = get_param(status, names::targets, gids);

        if (targets_changed)
        {
            recording_.clear();
            times_.clear();
        }

        // initialize the recording container
        for (size_t gid : gids)
        {
            NeuronPtr n = kernel().neuron_manager.get_neuron(gid);
            targets_.insert({gid, n});

            recording_.insert({gid, std::vector<double>()});
            times_.insert({gid, std::vector<double>()});
        }
    }
}


void NeuronDiscreteRecorder::record()
{}


unsigned int NeuronDiscreteRecorder::get_level() const
{
    return 0;
}


unsigned int NeuronDiscreteRecorder::get_event_type() const
{
    return 1;
}


bool NeuronDiscreteRecorder::get_next_recording(std::vector<Property>& ids,
                                                  std::vector<double>& values)
{
    if (!v_iterating_)
    {
        v_iterating_ = true;
        neuron_it_ = recording_.cbegin();
    }

    if (neuron_it_ != recording_.cend())
    {
        // set the values
        Property p = Property(neuron_it_->first);
        ids.push_back(p);

        // set the values
        values.reserve(values.size() + neuron_it_->second.size());
        values.insert(
            values.end(), neuron_it_->second.cbegin(),
            neuron_it_->second.cend());

        neuron_it_++;

        if (neuron_it_ != recording_.cend())
        {
            return true;
        }
        else
        {
            v_iterating_ = false;
            return false;
        }
    }
    else
    {
        v_iterating_ = false;
        return false;
    }
}


bool NeuronDiscreteRecorder::get_next_time(std::vector<Property>& ids,
                                             std::vector<double>& values)
{
    if (!t_iterating_)
    {
        t_iterating_ = true;
        time_it_ = times_.cbegin();
    }

    if (time_it_ != times_.cend())
    {
        // set the Property vector
        Property p = Property(time_it_->first);
        ids.push_back(p);

        // set the values
        values.reserve(values.size() + time_it_->second.size());
        values.insert(
            values.end(), time_it_->second.cbegin(), time_it_->second.cend());

        time_it_++;

        if (time_it_ != times_.cend())
        {
            return true;
        }
        else
        {
            t_iterating_ = false;
            return false;
        }
    }
    else
    {
        t_iterating_ = false;
        return false;
    }
}


/**
 * Constructor for NeuriteContinuousRecorder
 */
NeuriteContinuousRecorder::NeuriteContinuousRecorder()
{
    times_ = std::array<double, 3>({0., 0., 0.});
}


/**
 * @brief Set the status of "targets"
 *
 * Call default BaseRecorder::set_status, then do a custom update for targets
 */
void NeuriteContinuousRecorder::set_status(const statusMap &status)
{
    // call default parent set_status
    BaseRecorder::set_status(status);

    // set targets and recording_
    if (targets_.empty() || kernel().simulation_manager.get_current_step() == 0)
    {
        // get the targets
        std::vector<size_t> gids;
        bool targets_changed = get_param(status, names::targets, gids);

        if (targets_changed)
        {
            recording_.clear();
            times_ = std::array<double, 3>({0., 0., 0.});
        }

        // initialize the recording container
        for (size_t gid : gids)
        {
            NeuronPtr n = kernel().neuron_manager.get_neuron(gid);
            targets_.insert({gid, n});

            recording_.insert({gid, mapNameVecDouble()});

            //@todo: finish init
        }
    }
}


void NeuriteContinuousRecorder::record()
{
    if (record_to_ == "memory")
    {
        for (const auto& neuron : targets_)
        {
            for (const auto& neurite : neuron.second->neurites_)
            {
                recording_[neuron.first][neurite.first].push_back(
                    neurite.second->get_state(cstr_obs_));
            }
        }
        times_[2] += 1;
    }
    else
    {
        // @todo: see with Alessio
    }
}


/**
 * @brief Set the final times for the continuous records
 */
void NeuriteContinuousRecorder::final_timestep(size_t step)
{
    times_[1] = step;
}


unsigned int NeuriteContinuousRecorder::get_level() const
{
    return 1;
}


unsigned int NeuriteContinuousRecorder::get_event_type() const
{
    return 0;
}


bool NeuriteContinuousRecorder::get_next_recording(std::vector<Property>& ids,
                                                   std::vector<double>& values)
{
    if (!v_iterating_)
    {
        v_iterating_ = true;
        neuron_it_ = recording_.cbegin();
        if (neuron_it_ != recording_.cend())
        {
            neurite_it_ = neuron_it_->second.cbegin();
            neurite_endit_ = neuron_it_->second.cend();
        }
    }

    if (neuron_it_ != recording_.cend())
    {
        // check for neuron without neurites
        if (neurite_it_ != neurite_endit_)
        {
            // make the Property vector
            Property p_neuron = Property(neuron_it_->first);
            ids.push_back(p_neuron);
            Property p_neurite = Property(neurite_it_->first);
            ids.push_back(p_neurite);

            // set the values
            values.reserve(values.size() + neurite_it_->second.size());
            values.insert(
                values.end(), neurite_it_->second.cbegin(),
                neurite_it_->second.cend());

            // increment iterators
            neurite_it_++;
        }
        // increment all iterators
        if (neurite_it_ == neurite_endit_)
        {
            neuron_it_++;
            if (neuron_it_ != recording_.cend())
            {
                neurite_it_ = neuron_it_->second.cbegin();
                neurite_endit_ = neuron_it_->second.cend();
            }
            else
            {
                v_iterating_ = false;
                return false;
            }
        }
    }
    else
    {
        v_iterating_ = false;
        return false;
    }

    return true;
}


bool NeuriteContinuousRecorder::get_next_time(std::vector<Property>& ids,
                                              std::vector<double>& values)
{
    values.insert(values.end(), times_.begin(), times_.end());
    return false;
}


/**
 * Constructor for NeuriteContinuousRecorder
 */
NeuriteDiscreteRecorder::NeuriteDiscreteRecorder()
{}


void NeuriteDiscreteRecorder::record()
{}


unsigned int NeuriteDiscreteRecorder::get_level() const
{
    return 1;
}


unsigned int NeuriteDiscreteRecorder::get_event_type() const
{
    return 1;
}


/**
 * @brief Set the status of "targets"
 *
 * Call default BaseRecorder::set_status, then do a custom update for targets
 */
void NeuriteDiscreteRecorder::set_status(const statusMap &status)
{
    // call default parent set_status
    BaseRecorder::set_status(status);

    // set targets and recording_
    if (targets_.empty() || kernel().simulation_manager.get_current_step() == 0)
    {
        // get the targets
        std::vector<size_t> gids;
        bool targets_changed = get_param(status, names::targets, gids);

        if (targets_changed)
        {
            recording_.clear();
            times_.clear();
        }

        // initialize the recording container
        for (size_t gid : gids)
        {
            NeuronPtr n = kernel().neuron_manager.get_neuron(gid);
            //@todo: finish init
        }
    }
}


bool NeuriteDiscreteRecorder::get_next_recording(std::vector<Property>& ids,
                                                 std::vector<double>& values)
{
    if (!v_iterating_)
    {
        v_iterating_ = true;
        v_neuron_it_ = recording_.cbegin();
        if (v_neuron_it_ != recording_.cend())
        {
            v_neurite_it_ = v_neuron_it_->second.cbegin();
            v_neurite_endit_ = v_neuron_it_->second.cend();
        }
    }

    if (v_neuron_it_ != recording_.cend())
    {
        // check for neuron without neurites
        if (v_neurite_it_ != v_neurite_endit_)
        {
            // make the Property vector
            Property p_neuron = Property(v_neuron_it_->first);
            ids.push_back(p_neuron);
            Property p_neurite = Property(v_neurite_it_->first);
            ids.push_back(p_neurite);

            // set the values
            values.reserve(values.size() + v_neurite_it_->second.size());
            values.insert(
                values.end(), v_neurite_it_->second.cbegin(),
                v_neurite_it_->second.cend());

            // increment iterator
            v_neurite_it_++;
        }
        // increment all iterators
        if (v_neurite_it_ == v_neurite_endit_)
        {
            v_neuron_it_++;
            if (v_neuron_it_ != recording_.cend())
            {
                v_neurite_it_ = v_neuron_it_->second.cbegin();
                v_neurite_endit_ = v_neuron_it_->second.cend();
            }
            else
            {
                v_iterating_ = false;
                return false;
            }
        }
    }
    else
    {
        v_iterating_ = false;
        return false;
    }

    return true;
}


bool NeuriteDiscreteRecorder::get_next_time(std::vector<Property>& ids,
                                            std::vector<double>& values)
{
    if (!t_iterating_)
    {
        t_iterating_ = true;
        t_neuron_it_ = times_.cbegin();
        if (t_neuron_it_ != times_.cend())
        {
            t_neurite_it_ = t_neuron_it_->second.cbegin();
            t_neurite_endit_ = t_neuron_it_->second.cend();
        }
    }


    if (t_neuron_it_ != times_.cend())
    {
        // check for neurons without neurites
        if (t_neurite_it_!= t_neurite_endit_)
        {
            // make the Property vector
            Property p_neuron = Property(t_neuron_it_->first);
            ids.push_back(p_neuron);
            Property p_neurite = Property(t_neurite_it_->first);
            ids.push_back(p_neurite);

            // set the values
            values.reserve(values.size() + t_neurite_it_->second.size());
            values.insert(
                values.end(), t_neurite_it_->second.cbegin(),
                t_neurite_it_->second.cend());

            // increment iterator
            t_neurite_it_++;
        }
        // increment all iterators
        if (t_neurite_it_ == t_neurite_endit_)
        {
            t_neuron_it_++;
            if (t_neuron_it_ != times_.cend())
            {
                t_neurite_it_ = t_neuron_it_->second.cbegin();
                t_neurite_endit_ = t_neuron_it_->second.cend();
            }
            else
            {
                t_iterating_ = false;
                return false;
            }
        }
    }
    else
    {
        t_iterating_ = false;
        return false;
    }

    return true;
}


/**
 * Constructor for GrowthConeContinuousRecorder
 */
GrowthConeContinuousRecorder::GrowthConeContinuousRecorder()
  : v_gc_pos_(0)
  , v_gc_endpos_(0)
  , t_gc_pos_(0)
  , t_gc_endpos_(0)
{}


void GrowthConeContinuousRecorder::record()
{
    double step = kernel().simulation_manager.get_current_step();
    if (record_to_ == "memory")
    {
        for (const auto& neuron : targets_)
        {
            for (const auto& neurite : neuron.second->neurites_)
            {
                std::vector< std::vector<double> >& gc_values =
                    recording_[neuron.first][neurite.first];
                std::vector< std::array<double, 3> >& gc_times =
                    times_[neuron.first][neurite.first];

                for (const auto& gc : neurite.second->growth_cones_)
                {
                    if (gc_values.size() <= gc.first)
                    {
                        gc_values.push_back(std::vector<double>(
                            {gc.second->get_state(cstr_obs_)}));
                        gc_times.push_back(
                            std::array<double, 3>({step, step, 1.}));
                    }
                    else
                    {
                        gc_values[gc.first].push_back(
                            gc.second->get_state(cstr_obs_));
                        gc_times[gc.first][2] += 1;
                    }
                }
            }
        }
    }
    else
    {
        // @todo: see with Alessio
    }
}


/**
 * @brief Set the final times for the continuous records
 */
void GrowthConeContinuousRecorder::final_timestep(size_t step)
{
    for (auto& neuron_it : times_)
    {
        for (auto& neurite_it : neuron_it.second)
        {
            for (auto ttt : neurite_it.second)
            {
                ttt[1] = step;
            }
        }
    }
}


unsigned int GrowthConeContinuousRecorder::get_level() const
{
    return 2;
}


unsigned int GrowthConeContinuousRecorder::get_event_type() const
{
    return 0;
}


/**
 * @brief Set the status of "targets"
 *
 * Call default BaseRecorder::set_status, then do a custom update for targets
 */
void GrowthConeContinuousRecorder::set_status(const statusMap &status)
{
    // call default parent set_status
    BaseRecorder::set_status(status);
    printf("Set BaseRecorder status\n");

    size_t step = kernel().simulation_manager.get_current_step();
    double dstep = static_cast<double>(step);

    // set targets and recording_
    if (targets_.empty() || step == 0)
    {
        // get the targets
        std::vector<size_t> gids;
        bool targets_changed = get_param(status, names::targets, gids);

        if (targets_changed)
        {
            recording_.clear();
            times_.clear();
        }
        printf("Got targets and cleared\n");

        // initialize the recording container
        for (size_t gid : gids)
        {
            NeuronPtr n = kernel().neuron_manager.get_neuron(gid);

            recording_[gid] = std::unordered_map<std::string, vVecDouble>();
            times_[gid]     = std::unordered_map<std::string, vArray3Double>();
            for (const auto& neurite : n->neurites_)
            {
                printf("Creating neurite %s level containers\n", neurite.first.c_str());
                size_t size = neurite.second->growth_cones_.size();

                recording_[gid][neurite.first] = vVecDouble(size);
                times_[gid][neurite.first]     = vArray3Double(
                    size, {dstep, dstep, 0.});
            }
        }
    }
}


bool GrowthConeContinuousRecorder::get_next_recording(
    std::vector<Property>& ids, std::vector<double>& values)
{
    if (!v_iterating_)
    {
        v_iterating_ = true;
        v_neuron_it_ = recording_.cbegin();
        if (v_neuron_it_ != recording_.cend())
        {
            v_neurite_it_ = v_neuron_it_->second.cbegin();
            v_neurite_endit_ = v_neuron_it_->second.cend();
            if (v_neurite_it_ != v_neurite_endit_)
            {
                v_gc_pos_ = 0;
                v_gc_endpos_ = v_neurite_it_->second.size();
            }
        }
    }

    if (v_neuron_it_ != recording_.cend())
    {
        // check for neuron without neurites
        if (v_neurite_it_ != v_neurite_endit_)
        {
            // check for dead neurites
            if (v_gc_pos_ != v_gc_endpos_)
            {
                // make the Property vector
                Property p_neuron = Property(v_neuron_it_->first);
                ids.push_back(p_neuron);
                Property p_neurite = Property(v_neurite_it_->first);
                ids.push_back(p_neurite);
                Property p_gc = Property(v_gc_pos_);
                ids.push_back(p_gc);

                // set the values
                values.reserve(
                    values.size() + v_neurite_it_->second[v_gc_pos_].size());
                values.insert(
                    values.end(), v_neurite_it_->second[v_gc_pos_].cbegin(),
                    v_neurite_it_->second[v_gc_pos_].cend());

                // increment gc
                v_gc_pos_++;
            }

            // increment neurite after gc
            if (v_gc_pos_ == v_gc_endpos_)
            {
                v_neurite_it_++;
                if (v_neurite_it_ != v_neurite_endit_)
                {
                    v_gc_pos_ = 0;
                    v_gc_endpos_ = v_neurite_it_->second.size();
                }
            }
        }

        // increment neuron after neurite
        if (v_neurite_it_ == v_neurite_endit_)
        {
            v_neuron_it_++;
            if (v_neuron_it_ != recording_.cend())
            {
                v_neurite_it_ = v_neuron_it_->second.cbegin();
                v_neurite_endit_ = v_neuron_it_->second.cend();
                if (v_neurite_it_ != v_neurite_endit_)
                {
                    v_gc_pos_ = 0;
                    v_gc_endpos_ = v_neurite_it_->second.size();
                }
            }
            else
            {
                v_iterating_ = false;
                return false;
            }
        }
    }
    else
    {
        v_iterating_ = false;
        return false;
    }

    return true;
}


bool GrowthConeContinuousRecorder::get_next_time(
    std::vector<Property>& ids, std::vector<double>& values)
{
    if (!t_iterating_)
    {
        t_iterating_ = true;
        t_neuron_it_ = times_.cbegin();
        if (t_neuron_it_ != times_.cend())
        {
            t_neurite_it_ = t_neuron_it_->second.cbegin();
            t_neurite_endit_ = t_neuron_it_->second.cend();
            if (t_neurite_it_ != t_neurite_endit_)
            {
                t_gc_pos_ = 0;
                t_gc_endpos_ = t_neurite_it_->second.size();
            }
        }
    }

    if (t_neuron_it_ != times_.cend())
    {
        // check for neuron without neurites
        if (t_neurite_it_ != t_neurite_endit_)
        {
            // check for dead neurites
            if (t_gc_pos_ != t_gc_endpos_)
            {
                // make the Property vector
                Property p_neuron = Property(t_neuron_it_->first);
                ids.push_back(p_neuron);
                Property p_neurite = Property(t_neurite_it_->first);
                ids.push_back(p_neurite);
                Property p_gc = Property(t_gc_pos_);
                ids.push_back(p_gc);

                // set the values
                values.insert(values.end(),
                              t_neurite_it_->second[t_gc_pos_].begin(),
                              t_neurite_it_->second[t_gc_pos_].end());

                // increment gc
                t_gc_pos_++;
            }

            // increment neurite after gc
            if (t_gc_pos_ == t_gc_endpos_)
            {
                t_neurite_it_++;
                if (t_neurite_it_ != t_neurite_endit_)
                {
                    t_gc_pos_ = 0;
                    t_gc_endpos_ = t_neurite_it_->second.size();
                }
            }
        }

        // increment neuron after neurite
        if (t_neurite_it_ == t_neurite_endit_)
        {
            t_neuron_it_++;
            if (t_neuron_it_ != times_.cend())
            {
                t_neurite_it_ = t_neuron_it_->second.cbegin();
                t_neurite_endit_ = t_neuron_it_->second.cend();
                if (t_neurite_it_ != t_neurite_endit_)
                {
                    t_gc_pos_ = 0;
                    t_gc_endpos_ = t_neurite_it_->second.size();
                }
            }
            else
            {
                t_iterating_ = false;
                return false;
            }
        }
    }
    else
    {
        t_iterating_ = false;
        return false;
    }

    return true;
}


/**
 * Constructor for GrowthConeDiscreteRecorder
 */
GrowthConeDiscreteRecorder::GrowthConeDiscreteRecorder()
  : v_gc_pos_(0)
  , v_gc_endpos_(0)
  , t_gc_pos_(0)
  , t_gc_endpos_(0)
{}


void GrowthConeDiscreteRecorder::record()
{}


unsigned int GrowthConeDiscreteRecorder::get_level() const
{
    return 2;
}


unsigned int GrowthConeDiscreteRecorder::get_event_type() const
{
    return 1;
}


/**
 * @brief Set the status of "targets"
 *
 * Call default BaseRecorder::set_status, then do a custom update for targets
 */
void GrowthConeDiscreteRecorder::set_status(const statusMap &status)
{
    // call default parent set_status
    BaseRecorder::set_status(status);

    // set targets and recording_
    if (targets_.empty() || kernel().simulation_manager.get_current_step() == 0)
    {
        // get the targets
        std::vector<size_t> gids;
        bool targets_changed = get_param(status, names::targets, gids);

        if (targets_changed)
        {
            recording_.clear();
            times_.clear();
        }

        // initialize the recording container
        for (size_t gid : gids)
        {
            NeuronPtr n = kernel().neuron_manager.get_neuron(gid);
            //@todo: finish init
        }
    }
}


bool GrowthConeDiscreteRecorder::get_next_recording(
    std::vector<Property>& ids, std::vector<double>& values)
{
    if (!v_iterating_)
    {
        v_iterating_ = true;
        v_neuron_it_ = recording_.cbegin();
        if (v_neuron_it_ != recording_.cend())
        {
            v_neurite_it_ = v_neuron_it_->second.cbegin();
            v_neurite_endit_ = v_neuron_it_->second.cend();
            if (v_neurite_it_ != v_neurite_endit_)
            {
                v_gc_pos_ = 0;
                v_gc_endpos_ = v_neurite_it_->second.size();
            }
        }
    }

    if (v_neuron_it_ != recording_.cend())
    {
        // check for neuron without neurites
        if (v_neurite_it_ != v_neurite_endit_)
        {
            // check for dead neurites
            if (v_gc_pos_ != v_gc_endpos_)
            {
                // make the Property vector
                Property p_neuron = Property(v_neuron_it_->first);
                ids.push_back(p_neuron);
                Property p_neurite = Property(v_neurite_it_->first);
                ids.push_back(p_neurite);
                Property p_gc = Property(v_gc_pos_);
                ids.push_back(p_gc);

                // set the values
                values.reserve(
                    values.size() + v_neurite_it_->second[v_gc_pos_].size());
                values.insert(
                    values.end(), v_neurite_it_->second[v_gc_pos_].cbegin(),
                    v_neurite_it_->second[v_gc_pos_].cend());

                // increment gc
                v_gc_pos_++;
            }

            // increment neurite after gc
            if (v_gc_pos_ == v_gc_endpos_)
            {
                v_neurite_it_++;
                if (v_neurite_it_ != v_neurite_endit_)
                {
                    v_gc_pos_ = 0;
                    v_gc_endpos_ = v_neurite_it_->second.size();
                }
            }
        }

        // increment neuron after neurite
        if (v_neurite_it_ == v_neurite_endit_)
        {
            v_neuron_it_++;
            if (v_neuron_it_ != recording_.cend())
            {
                v_neurite_it_ = v_neuron_it_->second.cbegin();
                v_neurite_endit_ = v_neuron_it_->second.cend();
                if (v_neurite_it_ != v_neurite_endit_)
                {
                    v_gc_pos_ = 0;
                    v_gc_endpos_ = v_neurite_it_->second.size();
                }
            }
            else
            {
                v_iterating_ = false;
                return false;
            }
        }
    }
    else
    {
        v_iterating_ = false;
        return false;
    }

    return true;
}


bool GrowthConeDiscreteRecorder::get_next_time(
    std::vector<Property>& ids, std::vector<double>& values)
{
    if (!t_iterating_)
    {
        t_iterating_ = true;
        t_neuron_it_ = times_.cbegin();
        if (t_neuron_it_ != times_.cend())
        {
            t_neurite_it_ = t_neuron_it_->second.cbegin();
            t_neurite_endit_ = t_neuron_it_->second.cend();
            if (t_neurite_it_ != t_neurite_endit_)
            {
                t_gc_pos_ = 0;
                t_gc_endpos_ = t_neurite_it_->second.size();
            }
        }
    }

    if (t_neuron_it_ != times_.cend())
    {
        // check for neuron without neurites
        if (t_neurite_it_ != t_neurite_endit_)
        {
            // check for dead neurites
            if (t_gc_pos_ != t_gc_endpos_)
            {
                // make the Property vector
                Property p_neuron = Property(t_neuron_it_->first);
                ids.push_back(p_neuron);
                Property p_neurite = Property(t_neurite_it_->first);
                ids.push_back(p_neurite);
                Property p_gc = Property(t_gc_pos_);
                ids.push_back(p_gc);

                // set the values
                values.reserve(
                    values.size() + t_neurite_it_->second[t_gc_pos_].size());
                values.insert(
                    values.end(), t_neurite_it_->second[t_gc_pos_].cbegin(),
                    t_neurite_it_->second[t_gc_pos_].cend());

                // increment gc
                t_gc_pos_++;
            }

            // increment neurite after gc
            if (t_gc_pos_ == t_gc_endpos_)
            {
                t_neurite_it_++;
                if (t_neurite_it_ != t_neurite_endit_)
                {
                    t_gc_pos_ = 0;
                    t_gc_endpos_ = t_neurite_it_->second.size();
                }
            }
        }

        // increment neuron after neurite
        if (t_neurite_it_ == t_neurite_endit_)
        {
            t_neuron_it_++;
            if (t_neuron_it_ != times_.cend())
            {
                t_neurite_it_ = t_neuron_it_->second.cbegin();
                t_neurite_endit_ = t_neuron_it_->second.cend();
                if (t_neurite_it_ != t_neurite_endit_)
                {
                    t_gc_pos_ = 0;
                    t_gc_endpos_ = t_neurite_it_->second.size();
                }
            }
            else
            {
                t_iterating_ = false;
                return false;
            }
        }
    }
    else
    {
        t_iterating_ = false;
        return false;
    }

    return true;
}


/*
 * TOOLS
 */


} /* namespace */
