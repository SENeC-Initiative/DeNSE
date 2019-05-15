#ifndef NEURON_M_H
#define NEURON_M_H

#include <string>
#include <unordered_map>

#include "config.hpp"
#include "elements_types.hpp"
#include "manager_interface.hpp"

namespace growth
{

typedef std::unordered_map<size_t, NeuronPtr> NeuronMap;
typedef std::unordered_map<std::string, GCPtr> ModelMap;
typedef std::vector<std::vector<NeuronPtr>> ThreadNeurons;

class Neuron;
class GrowthCone;

class NeuronManager : public ManagerInterface
{
  public:
    NeuronManager();

    virtual void initialize();
    virtual void finalize();

    /**
     * Create neurons.
     */
    size_t create_neurons(const std::vector<statusMap> &neuron_params,
                          const std::vector<statusMap> &axon_params,
                          const std::vector<statusMap> &dendrites_params);

    NeuronPtr get_neuron(size_t gid);
    void get_all_neurons(std::vector<NeuronPtr> &);
    std::vector<size_t> get_gids() const;
    std::vector<NeuronPtr> get_local_neurons(int local_thread_id);

    void init_neurons_on_thread(unsigned int num_local_threads);
    void update_kernel_variables();

    const statusMap get_neuron_status(size_t gid) const;
    const statusMap get_neurite_status(size_t gid,
                                       const std::string &type) const;

    bool is_neuron(size_t gid) const;

    size_t num_neurons() const;

    NeuronMap::const_iterator iter_neurons();

    GCPtr get_model(std::string model_name);
    void get_models(std::vector<std::string> &models);
    GCPtr get_default_model();

    //! here we call the models file and register each model in the map
    void register_model(std::string, GCPtr);
    ModelMap model_map_;

  private:
    NeuronMap neurons_;
    ThreadNeurons neurons_on_thread_;
};
void _fill_skel(const SkelNeurite &source_container,
                SkelNeurite &target_container, bool add_nan);
}

#endif // NEURON_M_H
