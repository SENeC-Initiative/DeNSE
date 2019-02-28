#include "recorders.hpp"

// element include
#include "GrowthCone.hpp"
#include "Neurite.hpp"
#include "Neuron.hpp"

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
    set_param(status, names::record_to, record_to_, "");
    set_param(status, names::interval, interval_, "minute");
    set_param(status, names::restrict_to, restrict_to_, "");
    set_param(status, names::observable, observable_, "");

    // get the neuron ids
    std::vector<size_t> gids;
    for (auto const &neuron : targets_)
    {
        gids.push_back(neuron.first);
    }
    set_param(status, names::targets, gids, "");

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


void BaseRecorder::record(const Event &ev) {}


/**
 * @brief Set the final times for the continuous recorders
 */
void BaseRecorder::final_timestep(size_t step) {}


unsigned int BaseRecorder::get_level() const { return 3; }


unsigned int BaseRecorder::get_event_type() const { return 2; }


bool BaseRecorder::get_next_recording(std::vector<Property> &ids,
                                      std::vector<double> &values)
{
    return false;
}


bool BaseRecorder::get_next_time(std::vector<Property> &ids,
                                 std::vector<double> &values,
                                 const std::string &time_units)
{
    return false;
}


void BaseRecorder::reset_iterations()
{
    v_iterating_ = false;
    t_iterating_ = false;
}


void BaseRecorder::neuron_deleted(size_t gid)
{
    auto it = targets_.find(gid);

    if (it != targets_.end())
    {
        targets_.erase(it);
    }
}


/**
 * Constructor for NeuronContinuousRecorder
 */
NeuronContinuousRecorder::NeuronContinuousRecorder()
    : num_times_(0)
{
    Time t0 = kernel().simulation_manager.get_time();
    times_  = std::array<Time, 2>({t0, t0});
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
            Time t0 = kernel().simulation_manager.get_time();
            times_  = std::array<Time, 2>({t0, t0});
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
        for (const auto &neuron : targets_)
        {
            recording_[neuron.first].push_back(
                neuron.second->get_state(cstr_obs_));
        }
        num_times_ += 1;
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
    times_[1] = kernel().simulation_manager.get_time();
}


/**
 * @brief Return the level at which the recorder works
 *
 * Returns 0 if recording at neuron level, 1 if recording at neurite level and
 * 2 for growth cone level.
 * This function returns 0.
 */
unsigned int NeuronContinuousRecorder::get_level() const { return 0; }


/**
 * @brief Return the type of event recorded by the recorder
 *
 * Returns 0/1 if continuous/discrete events are recorded.
 * This function returns 0.
 */
unsigned int NeuronContinuousRecorder::get_event_type() const { return 0; }


/**
 * @brief Access the data of the next object recorded, return false when done.
 *
 * Data about the next object is given through a list properties (gid of the
 * neuron, name of the neurite, id of the growth cone) and a list of values.
 * The function returns true if there is still data left, otherwise it returns
 * false. Once all data has been looped on, the iterators are reset to their
 * initial values.
 */
bool NeuronContinuousRecorder::get_next_recording(std::vector<Property> &ids,
                                                  std::vector<double> &values)
{
    if (!v_iterating_)
    {
        v_iterating_ = true;
        neuron_it_   = recording_.cbegin();
    }

    if (neuron_it_ != recording_.cend())
    {
        Property p = Property(neuron_it_->first, "");
        ids.push_back(p);
        values.reserve(values.size() + neuron_it_->second.size());
        values.insert(values.end(), neuron_it_->second.cbegin(),
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
bool NeuronContinuousRecorder::get_next_time(std::vector<Property> &ids,
                                             std::vector<double> &values,
                                             const std::string &time_units)
{
    double t0, tf;
    const char *ctu = time_units.c_str();

    TRIE(ctu)
    CASE("seconds")
    t0 = times_[0].get_total_seconds();
    tf = times_[1].get_total_seconds();
    CASE("minutes")
    t0 = times_[0].get_total_minutes();
    tf = times_[1].get_total_minutes();
    CASE("hours")
    t0 = times_[0].get_total_hours();
    tf = times_[1].get_total_hours();
    CASE("days")
    t0 = times_[0].get_total_days();
    tf = times_[1].get_total_days();
    ENDTRIE;

    values.push_back(t0);
    values.push_back(tf);
    values.push_back(num_times_);

    return false;
}


/**
 * Constructor for NeuronDiscreteRecorder
 */
NeuronDiscreteRecorder::NeuronDiscreteRecorder() {}


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

            Time t0      = kernel().simulation_manager.get_time();
            double n_gc0 = n->get_num_neurites();

            recording_.insert({gid, std::vector<double>({n_gc0})});
            times_.insert({gid, std::vector<Time>({t0})});
        }
    }
}


void NeuronDiscreteRecorder::record(const Event &ev)
{
    Time event_time     = std::get<edata::TIME>(ev);
    size_t neuron       = std::get<edata::NEURON>(ev);
    signed char ev_type = std::get<edata::EV_TYPE>(ev);

    // test which data is recorded

    bool branching_event =
        (ev_type == names::lateral_branching || ev_type == names::gc_splitting);

    // test which data is recorded

    TRIE(cstr_obs_)
    CASE("num_growth_cones")
    if (branching_event)
    {
        recording_[neuron].push_back(recording_[neuron].back() + 1);
    }
    else if (ev_type == names::gc_deletion)
    {
        recording_[neuron].push_back(recording_[neuron].back() - 1);
    }
    times_[neuron].push_back(event_time);
    ENDTRIE;
}


unsigned int NeuronDiscreteRecorder::get_level() const { return 0; }


unsigned int NeuronDiscreteRecorder::get_event_type() const { return 1; }


bool NeuronDiscreteRecorder::get_next_recording(std::vector<Property> &ids,
                                                std::vector<double> &values)
{
    if (!v_iterating_)
    {
        v_iterating_ = true;
        neuron_it_   = recording_.cbegin();
    }

    if (neuron_it_ != recording_.cend())
    {
        // set the values
        Property p = Property(neuron_it_->first, "");
        ids.push_back(p);

        // set the values
        values.reserve(values.size() + neuron_it_->second.size());
        values.insert(values.end(), neuron_it_->second.cbegin(),
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


bool NeuronDiscreteRecorder::get_next_time(std::vector<Property> &ids,
                                           std::vector<double> &values,
                                           const std::string &time_units)
{
    const char *ctu = time_units.c_str();

    if (!t_iterating_)
    {
        t_iterating_ = true;
        time_it_     = times_.cbegin();
    }

    if (time_it_ != times_.cend())
    {
        // set the Property vector
        Property p = Property(time_it_->first, "");
        ids.push_back(p);

        // set the values
        values.reserve(values.size() + time_it_->second.size());

        TRIE(ctu)
        CASE("seconds")
        for (auto t : time_it_->second)
        {
            values.push_back(t.get_total_seconds());
        }
        CASE("minutes")
        for (auto t : time_it_->second)
        {
            values.push_back(t.get_total_minutes());
        }
        CASE("hours")
        for (auto t : time_it_->second)
        {
            values.push_back(t.get_total_hours());
        }
        CASE("days")
        for (auto t : time_it_->second)
        {
            values.push_back(t.get_total_days());
        }
        ENDTRIE;

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
    : num_times_(0)
{
    Time t0 = kernel().simulation_manager.get_time();
    times_  = std::array<Time, 2>({t0, t0});
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
            Time t0 = kernel().simulation_manager.get_time();
            times_  = std::array<Time, 2>({t0, t0});
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
        for (const auto &neuron : targets_)
        {
            for (const auto &neurite : neuron.second->neurites_)
            {
                recording_[neuron.first][neurite.first].push_back(
                    neurite.second->get_state(cstr_obs_));
            }
        }
        num_times_ += 1;
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
    times_[1] = kernel().simulation_manager.get_time();
}


void NeuriteContinuousRecorder::new_neurite(size_t neuron,
                                            const std::string& neurite)
{
    // initialize the recording container
    recording_[neuron].insert({neurite, std::vector<double>()});
}


unsigned int NeuriteContinuousRecorder::get_level() const { return 1; }


unsigned int NeuriteContinuousRecorder::get_event_type() const { return 0; }


bool NeuriteContinuousRecorder::get_next_recording(std::vector<Property> &ids,
                                                   std::vector<double> &values)
{
    if (!v_iterating_)
    {
        v_iterating_ = true;
        neuron_it_   = recording_.cbegin();
        if (neuron_it_ != recording_.cend())
        {
            neurite_it_    = neuron_it_->second.cbegin();
            neurite_endit_ = neuron_it_->second.cend();
        }
    }

    if (neuron_it_ != recording_.cend())
    {
        // check for neuron without neurites
        if (neurite_it_ != neurite_endit_)
        {
            // make the Property vector
            Property p_neuron = Property(neuron_it_->first, "");
            ids.push_back(p_neuron);
            Property p_neurite = Property(neurite_it_->first, "");
            ids.push_back(p_neurite);

            // set the values
            values.reserve(values.size() + neurite_it_->second.size());
            values.insert(values.end(), neurite_it_->second.cbegin(),
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
                neurite_it_    = neuron_it_->second.cbegin();
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


bool NeuriteContinuousRecorder::get_next_time(std::vector<Property> &ids,
                                              std::vector<double> &values,
                                              const std::string &time_units)
{
    double t0, tf;
    const char *ctu = time_units.c_str();

    TRIE(ctu)
    CASE("seconds")
    t0 = times_[0].get_total_seconds();
    tf = times_[1].get_total_seconds();
    CASE("minutes")
    t0 = times_[0].get_total_minutes();
    tf = times_[1].get_total_minutes();
    CASE("hours")
    t0 = times_[0].get_total_hours();
    tf = times_[1].get_total_hours();
    CASE("days")
    t0 = times_[0].get_total_days();
    tf = times_[1].get_total_days();
    ENDTRIE;

    values.push_back(t0);
    values.push_back(tf);
    values.push_back(num_times_);

    return false;
}


/**
 * Constructor for NeuriteContinuousRecorder
 */
NeuriteDiscreteRecorder::NeuriteDiscreteRecorder() {}


void NeuriteDiscreteRecorder::record(const Event &ev)
{
    Time event_time     = std::get<edata::TIME>(ev);
    size_t neuron       = std::get<edata::NEURON>(ev);
    std::string neurite = std::get<edata::NEURITE>(ev);
    signed char ev_type = std::get<edata::EV_TYPE>(ev);

    // test which data is recorded

    bool branching_event =
        (ev_type == names::lateral_branching || ev_type == names::gc_splitting);

#ifndef NDEBUG
    printf("got branching event at %f\n", event_time.get_total_minutes());
#endif

    TRIE(cstr_obs_)
    CASE("num_growth_cones")
    double old_val = recording_[neuron][neurite].back();

    if (branching_event)
    {
        recording_[neuron][neurite].push_back(old_val + 1);
    }
    else if (ev_type == names::gc_deletion)
    {
        recording_[neuron][neurite].push_back(old_val - 1);
    }

    times_[neuron][neurite].push_back(event_time);
    ENDTRIE;
}


unsigned int NeuriteDiscreteRecorder::get_level() const { return 1; }


unsigned int NeuriteDiscreteRecorder::get_event_type() const { return 1; }


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
            targets_.insert({gid, n});

            recording_.insert({gid, mapNameVecDouble()});
            times_.insert({gid, mapNameVecTime()});

            if (observable_ == names::num_growth_cones)
            {
                // initialize all num_growth_cones at 1 at initial_time
                auto it_neurite = n->neurite_cbegin();
                while (it_neurite != n->neurite_cend())
                {
                    Time t0 = kernel().simulation_manager.get_time();
                    recording_[gid][it_neurite->first] =
                        std::vector<double>({1});
                    times_[gid][it_neurite->first] = std::vector<Time>({t0});
                    it_neurite++;
                }
            }
        }
    }
}


void NeuriteDiscreteRecorder::new_neurite(size_t neuron,
                                          const std::string& neurite)
{
    // initialize the recording container
    if (observable_ == names::num_growth_cones)
    {
        Time t0 = kernel().simulation_manager.get_time();
        times_[neuron].insert({neurite, std::vector<Time>({t0})});
        recording_[neuron].insert({neurite, std::vector<double>({1.})});
    }
    else
    {
        times_[neuron].insert({neurite, std::vector<Time>()});
        recording_[neuron].insert({neurite, std::vector<double>()});
    }
}


bool NeuriteDiscreteRecorder::get_next_recording(std::vector<Property> &ids,
                                                 std::vector<double> &values)
{
    if (!v_iterating_)
    {
        v_iterating_ = true;
        v_neuron_it_ = recording_.cbegin();
        if (v_neuron_it_ != recording_.cend())
        {
            v_neurite_it_    = v_neuron_it_->second.cbegin();
            v_neurite_endit_ = v_neuron_it_->second.cend();
        }
    }

    if (v_neuron_it_ != recording_.cend())
    {
        // check for neuron without neurites
        if (v_neurite_it_ != v_neurite_endit_)
        {
            // make the Property vector
            Property p_neuron = Property(v_neuron_it_->first, "");
            ids.push_back(p_neuron);
            Property p_neurite = Property(v_neurite_it_->first, "");
            ids.push_back(p_neurite);

            // set the values
            values.reserve(values.size() + v_neurite_it_->second.size());
            values.insert(values.end(), v_neurite_it_->second.cbegin(),
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
                v_neurite_it_    = v_neuron_it_->second.cbegin();
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


bool NeuriteDiscreteRecorder::get_next_time(std::vector<Property> &ids,
                                            std::vector<double> &values,
                                            const std::string &time_units)
{
    const char *ctu = time_units.c_str();

    if (!t_iterating_)
    {
        t_iterating_ = true;
        t_neuron_it_ = times_.cbegin();
        if (t_neuron_it_ != times_.cend())
        {
            t_neurite_it_    = t_neuron_it_->second.cbegin();
            t_neurite_endit_ = t_neuron_it_->second.cend();
        }
    }


    if (t_neuron_it_ != times_.cend())
    {
        // check for neurons without neurites
        if (t_neurite_it_ != t_neurite_endit_)
        {
            // make the Property vector
            Property p_neuron = Property(t_neuron_it_->first, "");
            ids.push_back(p_neuron);
            Property p_neurite = Property(t_neurite_it_->first, "");
            ids.push_back(p_neurite);

            // set the values
            values.reserve(values.size() + t_neurite_it_->second.size());

            TRIE(ctu)
            CASE("seconds")
            for (auto t : t_neurite_it_->second)
            {
                values.push_back(t.get_total_seconds());
            }
            CASE("minutes")
            for (auto t : t_neurite_it_->second)
            {
                values.push_back(t.get_total_minutes());
            }
            CASE("hours")
            for (auto t : t_neurite_it_->second)
            {
                values.push_back(t.get_total_hours());
            }
            CASE("days")
            for (auto t : t_neurite_it_->second)
            {
                values.push_back(t.get_total_days());
            }
            ENDTRIE;

            // increment iterator
            t_neurite_it_++;
        }
        // increment all iterators
        if (t_neurite_it_ == t_neurite_endit_)
        {
            t_neuron_it_++;
            if (t_neuron_it_ != times_.cend())
            {
                t_neurite_it_    = t_neuron_it_->second.cbegin();
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
GrowthConeContinuousRecorder::GrowthConeContinuousRecorder() {}


void GrowthConeContinuousRecorder::record()
{
    if (record_to_ == "memory")
    {
        for (const auto &neuron : targets_)
        {
            for (const auto &neurite : neuron.second->neurites_)
            {
                auto &gc_values    = recording_[neuron.first][neurite.first];
                auto &gc_times     = times_[neuron.first][neurite.first];
                auto &gc_num_times = num_times_[neuron.first][neurite.first];

                for (const auto &gc : neurite.second->growth_cones_)
                {
                    auto it = gc_values.find(gc.first);
                    if (it == gc_values.end())
                    {
                        Time t0 = kernel().simulation_manager.get_time();
                        gc_times[gc.first]     = std::array<Time, 2>({t0, t0});
                        gc_num_times[gc.first] = 1;
                        gc_values[gc.first]    = std::vector<double>(
                            {gc.second->get_state(cstr_obs_)});
                    }
                    else
                    {
                        gc_values[gc.first].push_back(
                            gc.second->get_state(cstr_obs_));
                        gc_num_times[gc.first] += 1;
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
 * Take branching event into account
 */
void GrowthConeContinuousRecorder::record(const Event &ev)
{
    // branching event occured on a neuron
    Time event_time     = std::get<edata::TIME>(ev);
    size_t neuron       = std::get<edata::NEURON>(ev);
    std::string neurite = std::get<edata::NEURITE>(ev);

    std::unordered_map<size_t, std::vector<double>> &gc_values =
        recording_[neuron][neurite];
    std::unordered_map<size_t, std::array<Time, 2>> &gc_times =
        times_[neuron][neurite];
    std::unordered_map<size_t, size_t> &gc_num_times =
        num_times_[neuron][neurite];

    // get the neuron and neurite objects
    NeuronPtr n     = kernel().neuron_manager.get_neuron(neuron);
    auto neurite_it = n->neurites_.find(neurite);

    // add the new growth cone that just appeared
    for (const auto &gc : neurite_it->second->growth_cones_)
    {
        auto it = gc_values.find(gc.first);
        if (it == gc_values.end())
        {
            Time t0 = kernel().simulation_manager.get_time();
            gc_values[gc.first] =
                std::vector<double>({gc.second->get_state(cstr_obs_)});
            gc_times[gc.first] = std::array<Time, 2>({event_time, event_time});
            gc_num_times[gc.first] = 1;
        }
    }
}


/**
 * @brief Set the final times for the continuous records
 */
void GrowthConeContinuousRecorder::final_timestep(size_t step)
{
    for (auto &neuron_it : times_)
    {
        for (auto &neurite_it : neuron_it.second)
        {
            auto dead_set = dead_cones_[neuron_it.first][neurite_it.first];

            for (auto &ttt : neurite_it.second)
            {
                if (dead_set.count(ttt.first) == 0)
                {
                    ttt.second[1] = kernel().simulation_manager.get_time();
                }
            }
        }
    }
}


void GrowthConeContinuousRecorder::new_neurite(size_t neuron,
                                               const std::string& neurite)
{
    times_[neuron].insert({neurite, mapNumArrayTime()});
    recording_[neuron].insert({neurite, mapNumVecDouble()});
}


void GrowthConeContinuousRecorder::gc_died(size_t neuron,
                                           const std::string& neurite,
                                           size_t gc_id)
{
    // add cone to dead cones
    dead_cones_[neuron][neurite].insert(gc_id);

    // set last time value
    times_[neuron][neurite][gc_id][1] = kernel().simulation_manager.get_time();
}


unsigned int GrowthConeContinuousRecorder::get_level() const { return 2; }


unsigned int GrowthConeContinuousRecorder::get_event_type() const { return 0; }


/**
 * @brief Set the status of "targets"
 *
 * Call default BaseRecorder::set_status, then do a custom update for targets
 */
void GrowthConeContinuousRecorder::set_status(const statusMap &status)
{
    // call default parent set_status
    BaseRecorder::set_status(status);

    Time t = kernel().simulation_manager.get_time();

    // set targets and recording_
    if (targets_.empty() || t.get_total_seconds() == 0.)
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

            recording_[gid] =
                std::unordered_map<std::string, mapNumVecDouble>();
            times_[gid] = std::unordered_map<std::string, mapNumArrayTime>();
            num_times_[gid] =
                std::unordered_map<std::string,
                                   std::unordered_map<size_t, size_t>>();
            dead_cones_[gid] =
                std::unordered_map< std::string, std::unordered_set<size_t> >();

            for (const auto &neurite : n->neurites_)
            {
                size_t size = neurite.second->growth_cones_.size();

                recording_[gid][neurite.first]  = mapNumVecDouble();
                times_[gid][neurite.first]      = mapNumArrayTime();
                dead_cones_[gid][neurite.first] = std::unordered_set<size_t>();
                num_times_[gid][neurite.first]  =
                    std::unordered_map<size_t, size_t>();
            }
        }
    }
}


bool GrowthConeContinuousRecorder::get_next_recording(
    std::vector<Property> &ids, std::vector<double> &values)
{
    if (!v_iterating_)
    {
        v_iterating_ = true;
        v_neuron_it_ = recording_.cbegin();
        if (v_neuron_it_ != recording_.cend())
        {
            v_neurite_it_    = v_neuron_it_->second.cbegin();
            v_neurite_endit_ = v_neuron_it_->second.cend();
            if (v_neurite_it_ != v_neurite_endit_)
            {
                v_gc_pos_    = v_neurite_it_->second.cbegin();
                v_gc_endpos_ = v_neurite_it_->second.cend();
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
                Property p_neuron = Property(v_neuron_it_->first, "");
                ids.push_back(p_neuron);
                Property p_neurite = Property(v_neurite_it_->first, "");
                ids.push_back(p_neurite);
                Property p_gc = Property(v_gc_pos_->first, "");
                ids.push_back(p_gc);

                // set the values
                values.reserve(
                    values.size() +
                    v_neurite_it_->second.at(v_gc_pos_->first).size());
                values.insert(
                    values.end(),
                    v_neurite_it_->second.at(v_gc_pos_->first).cbegin(),
                    v_neurite_it_->second.at(v_gc_pos_->first).cend());

                // increment gc
                v_gc_pos_++;
            }

            // increment neurite after gc
            if (v_gc_pos_ == v_gc_endpos_)
            {
                v_neurite_it_++;
                if (v_neurite_it_ != v_neurite_endit_)
                {
                    v_gc_pos_    = v_neurite_it_->second.cbegin();
                    v_gc_endpos_ = v_neurite_it_->second.cend();
                }
            }
        }

        // increment neuron after neurite
        if (v_neurite_it_ == v_neurite_endit_)
        {
            v_neuron_it_++;
            if (v_neuron_it_ != recording_.cend())
            {
                v_neurite_it_    = v_neuron_it_->second.cbegin();
                v_neurite_endit_ = v_neuron_it_->second.cend();
                if (v_neurite_it_ != v_neurite_endit_)
                {
                    v_gc_pos_    = v_neurite_it_->second.cbegin();
                    v_gc_endpos_ = v_neurite_it_->second.cend();
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


bool GrowthConeContinuousRecorder::get_next_time(std::vector<Property> &ids,
                                                 std::vector<double> &values,
                                                 const std::string &time_units)
{
    const char *ctu = time_units.c_str();

    if (!t_iterating_)
    {
        t_iterating_ = true;
        t_neuron_it_ = times_.cbegin();
        if (t_neuron_it_ != times_.cend())
        {
            t_neurite_it_    = t_neuron_it_->second.cbegin();
            t_neurite_endit_ = t_neuron_it_->second.cend();
            if (t_neurite_it_ != t_neurite_endit_)
            {
                t_gc_pos_    = t_neurite_it_->second.cbegin();
                t_gc_endpos_ = t_neurite_it_->second.cend();
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
                Property p_neuron = Property(t_neuron_it_->first, "");
                ids.push_back(p_neuron);
                Property p_neurite = Property(t_neurite_it_->first, "");
                ids.push_back(p_neurite);
                Property p_gc = Property(t_gc_pos_->first, "");
                ids.push_back(p_gc);

                // set the values
                TRIE(ctu)
                CASE("seconds")
                for (auto t : t_neurite_it_->second.at(t_gc_pos_->first))
                {
                    values.push_back(t.get_total_seconds());
                }
                CASE("minutes")
                for (auto t : t_neurite_it_->second.at(t_gc_pos_->first))
                {
                    values.push_back(t.get_total_minutes());
                }
                CASE("hours")
                for (auto t : t_neurite_it_->second.at(t_gc_pos_->first))
                {
                    values.push_back(t.get_total_hours());
                }
                CASE("days")
                for (auto t : t_neurite_it_->second.at(t_gc_pos_->first))
                {
                    values.push_back(t.get_total_days());
                }
                ENDTRIE;
                // increment gc
                t_gc_pos_++;
            }

            // increment neurite after gc
            if (t_gc_pos_ == t_gc_endpos_)
            {
                t_neurite_it_++;
                if (t_neurite_it_ != t_neurite_endit_)
                {
                    t_gc_pos_    = t_neurite_it_->second.cbegin();
                    t_gc_endpos_ = t_neurite_it_->second.cend();
                }
            }
        }

        // increment neuron after neurite
        if (t_neurite_it_ == t_neurite_endit_)
        {
            t_neuron_it_++;
            if (t_neuron_it_ != times_.cend())
            {
                t_neurite_it_    = t_neuron_it_->second.cbegin();
                t_neurite_endit_ = t_neuron_it_->second.cend();
                if (t_neurite_it_ != t_neurite_endit_)
                {
                    t_gc_pos_    = t_neurite_it_->second.cbegin();
                    t_gc_endpos_ = t_neurite_it_->second.cend();
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
GrowthConeDiscreteRecorder::GrowthConeDiscreteRecorder() {}


void GrowthConeDiscreteRecorder::record(const Event &event)
{
    //~ size_t step         = std::get<0>(event);
    //~ double substep      = std::get<1>(event);
    //~ size_t neuron       = std::get<2>(event);
    //~ signed char ev_type = std::get<4>(event);

    //~ Time event_time = Time::from_steps(step, substep);

    //~ bool branching_event = (ev_type == names::lateral_branching
    //~ || ev_type == names::gc_splitting)

    //~ // test which data is recorded

    //~ TRIE(observable_)
    //~ CASE(names::num_growth_cones)
    //~ if (branching_event)
    //~ {
    //~ recording_.push_back();
    //~ }
    //~ else if (ev_type == names::gc_deletion)
    //~ {
    //~ recording_.push_back();
    //~ }
    //~ times_.push_back();
    //~ ENDTRIE;
}


unsigned int GrowthConeDiscreteRecorder::get_level() const { return 2; }


unsigned int GrowthConeDiscreteRecorder::get_event_type() const { return 1; }


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


bool GrowthConeDiscreteRecorder::get_next_recording(std::vector<Property> &ids,
                                                    std::vector<double> &values)
{
    if (!v_iterating_)
    {
        v_iterating_ = true;
        v_neuron_it_ = recording_.cbegin();
        if (v_neuron_it_ != recording_.cend())
        {
            v_neurite_it_    = v_neuron_it_->second.cbegin();
            v_neurite_endit_ = v_neuron_it_->second.cend();
            if (v_neurite_it_ != v_neurite_endit_)
            {
                v_gc_pos_    = v_neurite_it_->second.cbegin();
                v_gc_endpos_ = v_neurite_it_->second.cend();
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
                Property p_neuron = Property(v_neuron_it_->first, "");
                ids.push_back(p_neuron);
                Property p_neurite = Property(v_neurite_it_->first, "");
                ids.push_back(p_neurite);
                Property p_gc = Property(v_gc_pos_->first, "");
                ids.push_back(p_gc);

                // set the values
                values.reserve(
                    values.size() +
                    v_neurite_it_->second.at(v_gc_pos_->first).size());
                values.insert(
                    values.end(),
                    v_neurite_it_->second.at(v_gc_pos_->first).cbegin(),
                    v_neurite_it_->second.at(v_gc_pos_->first).cend());

                // increment gc
                v_gc_pos_++;
            }

            // increment neurite after gc
            if (v_gc_pos_ == v_gc_endpos_)
            {
                v_neurite_it_++;
                if (v_neurite_it_ != v_neurite_endit_)
                {
                    v_gc_pos_    = v_neurite_it_->second.cbegin();
                    v_gc_endpos_ = v_neurite_it_->second.cend();
                }
            }
        }

        // increment neuron after neurite
        if (v_neurite_it_ == v_neurite_endit_)
        {
            v_neuron_it_++;
            if (v_neuron_it_ != recording_.cend())
            {
                v_neurite_it_    = v_neuron_it_->second.cbegin();
                v_neurite_endit_ = v_neuron_it_->second.cend();
                if (v_neurite_it_ != v_neurite_endit_)
                {
                    v_gc_pos_    = v_neurite_it_->second.cbegin();
                    v_gc_endpos_ = v_neurite_it_->second.cend();
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


bool GrowthConeDiscreteRecorder::get_next_time(std::vector<Property> &ids,
                                               std::vector<double> &values,
                                               const std::string &time_units)
{
    const char *ctu = time_units.c_str();

    if (!t_iterating_)
    {
        t_iterating_ = true;
        t_neuron_it_ = times_.cbegin();
        if (t_neuron_it_ != times_.cend())
        {
            t_neurite_it_    = t_neuron_it_->second.cbegin();
            t_neurite_endit_ = t_neuron_it_->second.cend();
            if (t_neurite_it_ != t_neurite_endit_)
            {
                t_gc_pos_    = t_neurite_it_->second.cbegin();
                t_gc_endpos_ = t_neurite_it_->second.cend();
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
                Property p_neuron = Property(t_neuron_it_->first, "");
                ids.push_back(p_neuron);
                Property p_neurite = Property(t_neurite_it_->first, "");
                ids.push_back(p_neurite);
                Property p_gc = Property(t_gc_pos_->first, "");
                ids.push_back(p_gc);

                // set the values
                values.reserve(
                    values.size() +
                    t_neurite_it_->second.at(t_gc_pos_->first).size());

                TRIE(ctu)
                CASE("seconds")
                for (auto t : t_neurite_it_->second.at(t_gc_pos_->first))
                {
                    values.push_back(t.get_total_seconds());
                }
                CASE("minutes")
                for (auto t : t_neurite_it_->second.at(t_gc_pos_->first))
                {
                    values.push_back(t.get_total_minutes());
                }
                CASE("hours")
                for (auto t : t_neurite_it_->second.at(t_gc_pos_->first))
                {
                    values.push_back(t.get_total_hours());
                }
                CASE("days")
                for (auto t : t_neurite_it_->second.at(t_gc_pos_->first))
                {
                    values.push_back(t.get_total_days());
                }
                ENDTRIE;

                // increment gc
                t_gc_pos_++;
            }

            // increment neurite after gc
            if (t_gc_pos_ == t_gc_endpos_)
            {
                t_neurite_it_++;
                if (t_neurite_it_ != t_neurite_endit_)
                {
                    t_gc_pos_    = t_neurite_it_->second.cbegin();
                    t_gc_endpos_ = t_neurite_it_->second.cend();
                }
            }
        }

        // increment neuron after neurite
        if (t_neurite_it_ == t_neurite_endit_)
        {
            t_neuron_it_++;
            if (t_neuron_it_ != times_.cend())
            {
                t_neurite_it_    = t_neuron_it_->second.cbegin();
                t_neurite_endit_ = t_neuron_it_->second.cend();
                if (t_neurite_it_ != t_neurite_endit_)
                {
                    t_gc_pos_    = t_neurite_it_->second.cbegin();
                    t_gc_endpos_ = t_neurite_it_->second.cend();
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


} // namespace growth
