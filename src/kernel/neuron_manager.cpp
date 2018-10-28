#include "neuron_manager.hpp"

#include <cmath>
#include <cstring>
#include <limits>
#include <mutex>

#include "config_impl.hpp"
#include "kernel_manager.hpp"

#include "Neuron.hpp"


namespace growth
{

NeuronManager::NeuronManager() {}


void NeuronManager::initialize()
{
    // create default neuron with two neurites (axon + dendrite)
    statusMap empty_params;
    statusMap params({{names::num_neurites, Property(2, "")},
                      {names::growth_cone_model, Property("", "")},
                      {"x", Property(0., "micrometer")},
                      {"y", Property(0., "micrometer")}});

    // set their status to use all possible parameters to have them all
    // when using GetDefaults
    std::vector<std::string> options({"use_critical_resource", "use_van_pelt",
                                      "use_uniform_branching",
                                      "use_flpl_branching"});

    for (std::string opt : options)
    {
        Property p  = Property(true, "");
        params[opt] = p;
    }

    // get random engine
    mtPtr rnd_ptr = kernel().rng_manager.get_rng(0);

    // create default neuron
    model_neuron_ = std::make_shared<Neuron>(0);
    model_neuron_->init_status(params, empty_params, empty_params, rnd_ptr);
    // remove information about angles
    model_neuron_->neurite_angles_.clear();
}


void NeuronManager::finalize()
{
    neurons_.clear();
    neurons_on_thread_.clear();
}


/**
 * @brief Create new neurons with custom parameters.
 *
 * @return Number of objects created.
 */
size_t
NeuronManager::create_neurons(const std::vector<statusMap> &neuron_params,
                              const std::vector<statusMap> &axon_params,
                              const std::vector<statusMap> &dendrites_params)
{
    size_t first_id             = kernel().get_num_objects();
    size_t previous_num_neurons = neurons_.size();
    // put the neurons on the thread list they belong to
    size_t num_omp = kernel().parallelism_manager.get_num_local_threads();
    int omp_id     = kernel().parallelism_manager.get_thread_local_id();
    std::vector<std::vector<size_t>> thread_neurons(num_omp);
    for (size_t i = 0; i < neuron_params.size(); i++)
    {
        double x, y;
        get_param(neuron_params[i], "x", x);
        get_param(neuron_params[i], "y", y);

        if (kernel().space_manager.has_environment())
        {
            // printf("check point %f, %f in the environment \n",x,y);
            if (not kernel().space_manager.env_contains(Point(x, y), omp_id))
            {
                throw std::runtime_error(
                    " a Neuron was positioned outside the environment\n");
            }
        }

        // @TODO change temporary round-robin for neuron assignement
        unsigned int omp_id = (first_id + i) % num_omp;
        thread_neurons[omp_id].push_back(first_id + i);
        thread_of_neuron_[first_id + i] = omp_id;
    }

    // exception_capture_flag "guards" captured_exception. std::called_once()
    // guarantees that will only execute any of its Callable(s) ONCE for each
    // unique std::once_flag. See C++11 Standard Library documentation
    // (``<mutex>``). These tools together ensure that we can capture exceptions
    // from OpenMP parallel regions in a thread-safe way.
    std::once_flag exception_capture_flag;
    // pointer-like object that manages an exception captured with
    // std::capture_exception(). We use this to capture exceptions thrown from
    // the OpenMP parallel region.
    std::exception_ptr captured_exception;

// create the neurons on the respective threads
#pragma omp parallel
    {
        std::vector<NeuronPtr> local_neurons;
        int omp_id       = kernel().parallelism_manager.get_thread_local_id();
        mtPtr rnd_engine = kernel().rng_manager.get_rng(omp_id);
        std::vector<size_t> gids(thread_neurons[omp_id]);
        for (size_t gid : gids)
        {
            size_t idx       = gid - first_id;
            NeuronPtr neuron = std::make_shared<Neuron>(gid);

            try
            {
                neuron->init_status(neuron_params[idx], axon_params[idx],
                                    dendrites_params[idx], rnd_engine);
            }
            catch (const std::exception &except)
            {
                std::call_once(exception_capture_flag, [&captured_exception]() {
                    captured_exception = std::current_exception();
                });
            }

            local_neurons.push_back(neuron);
        }
#pragma omp critical
        {
            for (size_t i = 0; i < gids.size(); i++)
            {
                neurons_.insert({gids[i], local_neurons[i]});
                neurons_on_thread_[omp_id].push_back(local_neurons[i]);
            }
        }
    }

    // check if an exception was thrown there
    if (captured_exception != nullptr)
    {
        // rethrowing nullptr is illegal
        std::rethrow_exception(captured_exception);
    }

    // tell the kernel manager to update the number of objects
    kernel().update_num_objects();
    return neurons_.size() - previous_num_neurons;
}


void NeuronManager::init_neurons_on_thread(unsigned int num_local_threads)
{
    assert(neurons_.size() == 0); // no changes once neurons exist
    neurons_on_thread_ = std::vector<std::vector<NeuronPtr>>(num_local_threads);
    max_resolutions_   = std::vector<std::unordered_map<size_t, double>>(
        num_local_threads, std::unordered_map<size_t, double>());
}


void NeuronManager::update_kernel_variables()
{
    // exception_capture_flag "guards" captured_exception. std::called_once()
    // guarantees that will only execute any of its Callable(s) ONCE for each
    // unique std::once_flag. See C++11 Standard Library documentation
    // (``<mutex>``). These tools together ensure that we can capture exceptions
    // from OpenMP parallel regions in a thread-safe way.
    std::once_flag exception_capture_flag;
    // pointer-like object that manages an exception captured with
    // std::capture_exception(). We use this to capture exceptions thrown from
    // the OpenMP parallel region.
    std::exception_ptr captured_exception;

#pragma omp parallel
    {
        try
        {
            int omp_id = kernel().parallelism_manager.get_thread_local_id();

            gidNeuronMap local_neurons = get_local_neurons(omp_id);

            for (auto &neuron : local_neurons)
            {
                neuron.second->update_kernel_variables();
            }
        }
        catch (const std::exception &except)
        {
            std::call_once(exception_capture_flag, [&captured_exception]() {
                captured_exception = std::current_exception();
            });
        }
    }

    // check if an exception was thrown there
    if (captured_exception != nullptr)
    {
        // rethrowing nullptr is illegal
        std::rethrow_exception(captured_exception);
    }
}


// getters

NeuronPtr NeuronManager::get_neuron(size_t gid) { return neurons_[gid]; }


gidNeuronMap NeuronManager::get_local_neurons(int local_thread_id)
{
    gidNeuronMap local_neurons;

    for (auto &n : neurons_on_thread_[local_thread_id])
    {
        local_neurons[n->get_gid()] = n;
    }

    return local_neurons;
}


void NeuronManager::get_all_neurons(std::vector<NeuronPtr> &neuron_ptr_vec)
{
    for (const auto &neuron : neurons_)
    {
        neuron_ptr_vec.push_back(neuron.second);
    }
}


void NeuronManager::get_defaults(statusMap &status,
                                 const std::string &object,
                                 bool detailed) const
{
    if (object == "neuron")
    {
        if (detailed)
        {
            model_neuron_->get_neurite_status(status, "dendrite", "neurite");
            model_neuron_->get_neurite_status(status, "axon", "neurite");
        }
        model_neuron_->get_status(status);
    }
    else if (object == "axon")
    {
        model_neuron_->get_neurite_status(status, "axon", "neurite");
    }
    else if (object == "dendrite" || object == "neurite")
    {
        model_neuron_->get_neurite_status(status, "dendrite", "neurite");
    }
}


const statusMap NeuronManager::get_neuron_status(size_t gid) const
{
    NeuronPtr n = neurons_.at(gid);
    statusMap neuron_status;

    n->get_status(neuron_status);

    return neuron_status;
}


std::vector<size_t> NeuronManager::get_gids() const
{
    size_t n = neurons_.size();
    std::vector<size_t> gids;
    for (auto it = neurons_.begin(); it != neurons_.end(); it++)
        gids.push_back(it->first);
    return gids;
}


const statusMap
NeuronManager::get_neurite_status(size_t gid, const std::string &neurite,
                                  const std::string &level) const
{
    statusMap status;
    neurons_.at(gid)->get_neurite_status(status, neurite, level);
    return status;
}


bool NeuronManager::is_neuron(size_t gid) const
{
    auto it = neurons_.find(gid);

    if (it != neurons_.end())
        return true;
    return false;
}


gidNeuronMap::const_iterator NeuronManager::iter_neurons()
{
    return neurons_.begin();
}


size_t NeuronManager::num_neurons() const { return neurons_.size(); }


GCPtr NeuronManager::get_model(std::string model) { return model_map_[model]; }


GCPtr NeuronManager::get_default_model()
// default model is passed by model manager
{
    return model_map_["default"];
}


int NeuronManager::get_neuron_thread(size_t gid) const
{
    return thread_of_neuron_.at(gid);
}


void NeuronManager::register_model(std::string model_name, GCPtr model_ptr)
{
    model_map_[model_name] = model_ptr;
}


void NeuronManager::get_models(std::vector<std::string> &models)
{
    for (auto it : model_map_)
    {
        models.push_back(it.first);
    }
}


void NeuronManager::set_max_resol(size_t neuron, double max_resol)
{
    max_resolutions_[thread_of_neuron_[neuron]][neuron] = max_resol;
}


double NeuronManager::get_max_resol() const
{
    double max_resol = std::numeric_limits<double>::max();

    for (const auto &map : max_resolutions_)
    {
        for (const auto &pair : map)
        {
            max_resol = std::min(max_resol, pair.second);
        }
    }

    return max_resol;
}

} // namespace growth
