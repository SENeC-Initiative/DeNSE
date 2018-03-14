#ifndef NEURON_M_H
#define NEURON_M_H

#include <string>
#include <unordered_map>

#include "config.hpp"
#include "elements_types.hpp"
#include "manager_interface.hpp"

namespace growth
{

// typedefs
typedef std::unordered_map<size_t, NeuronPtr> gidNeuronMap;
typedef std::unordered_map<size_t, int> gidThreadMap;
typedef std::unordered_map<std::string, GCPtr> modelMap;
typedef std::vector<std::vector<NeuronPtr>> threadNeurons;


// forward declarations
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
    gidNeuronMap get_local_neurons(int local_thread_id);
    int get_neuron_thread(size_t gid) const;

    void init_neurons_on_thread(unsigned int num_local_threads);
    void update_kernel_variables();

    void get_defaults(statusMap &status, const std::string &object) const;
    const statusMap get_neuron_status(size_t gid) const;
    const statusMap get_neurite_status(size_t gid,
                                       const std::string &type,
                                       const std::string& level) const;

    bool is_neuron(size_t gid) const;

    size_t num_neurons() const;

    gidNeuronMap::const_iterator iter_neurons();

    GCPtr get_model(std::string model_name);
    void get_models(std::vector<std::string> &models);
    GCPtr get_default_model();

    void set_max_resol(size_t neuron, double max_resol);
    double get_max_resol() const;

    //! here we call the models file and register each model in the map
    void register_model(std::string, GCPtr);
    modelMap model_map_;

  private:
    NeuronPtr model_neuron_;           // unused model neuron for get_defaults
    gidNeuronMap neurons_;             // get neuron from gid
    threadNeurons neurons_on_thread_;  // group neurons by thread
    gidThreadMap thread_of_neuron_;    // get thread from gid
    std::vector<std::unordered_map<size_t, double>> max_resolutions_;  // max allowed resol
};


void _fill_skel(const SkelNeurite &source_container,
                SkelNeurite &target_container, bool add_nan);
}

#endif // NEURON_M_H
