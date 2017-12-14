#include "Neuron.hpp"

// c++ includes
#include <cmath>
#include <functional>
#include <stdexcept>

// elements includes
#include "GrowthCone.hpp"
#include "Node.hpp"
#include "growth_names.hpp"

// kernel includes
#include "Skeleton.hpp"
#include "Swc.hpp"
#include "config_impl.hpp"
#include "kernel_manager.hpp"
#include "neuron_manager.hpp"


namespace growth
{

//--------------------------------------------------------------------------
// Neuron class

// Constructors and destructor

//! Default constructor
Neuron::Neuron(size_t gid)
    : gid_(gid)
    , details()
    , observables_({"length", "speed"})
    , use_actin_waves_(false)
    , aw_generation_step_(-1)
    , actin_content_(0)
    , next_actin_event_(0)
    , axon_angle_(0)
    , axon_angle_set_(false)
{
    // initialize the soma
    soma_          = std::make_shared<BaseNode>();
    soma_->somaID_ = "soma";
    uniform_       = std::uniform_real_distribution<double>(0., 1.);
    normal_        = std::normal_distribution<double>(0., 1.);
}


//###################################################
//                  Init functions
//###################################################

void Neuron::init_status(const statusMap &status, const statusMap &astatus,
                         const statusMap &dstatus, mtPtr rnd_engine)
{
    // set growth cone model and actin wave
    set_status(status);
    int num_neurites = 0;
    get_param(status, names::num_neurites, num_neurites);
    get_param(status, names::growth_cone_model, growth_cone_model_);
    GCPtr gc_model;
    // get growth cone model from the model manager;
    if (growth_cone_model_ != "")
    {
        gc_model = kernel().neuron_manager.get_model(growth_cone_model_);
#ifndef NDEBUG
        printf("model is %s :", growth_cone_model_.c_str());
#endif
    }
    else
    {
#ifndef NDEBUG
        printf("WARNING: DEFAULT MODEL IN USE--> DEFINE A MODEL! <-- \n");
#endif
        gc_model = kernel().neuron_manager.get_default_model();
    }
#ifndef NDEBUG

    // this neuron will be initialized with this property and nothing else
    printf("\n#########################\n"
           "Neuron initialization: \n"
           "soma radius: %f  \n",
           "growth cone model:  %s\n"
           "use actin waves: %d \n",
           "num_neurites: %d \n", "#########################\n",
           details.soma_radius, growth_cone_model_.c_str(), use_actin_waves_);
#endif
    // initialize the soma giving the position of neuron
    //@TODO move the neuron in a new position with all the branches
    auto pos_x = status.find("x");
    if (pos_x != status.end())
    {
        double x, y;
        get_param(status, "x", x);
        get_param(status, "y", y);
        soma_->set_position(Point(x, y));
    }
    else
    {
        throw InvalidParameter(
            "Position was not set.", __FUNCTION__, __FILE__, __LINE__);
    }

    // create neurites
    // first angle nd growth cone model are the only non default params
    // check if customs growth_cone_models are required
    GCPtr axon_gc     = gc_model;
    GCPtr dendrite_gc = gc_model;
    std::string model_name;
    auto it = astatus.find(names::growth_cone_model);
    if (it != astatus.end())
    {
        get_param(astatus, names::growth_cone_model, model_name);
        axon_gc = kernel().neuron_manager.get_model(model_name);
    }
    it = dstatus.find(names::growth_cone_model);
    if (it != dstatus.end())
    {
        get_param(astatus, names::growth_cone_model, model_name);
        dendrite_gc = kernel().neuron_manager.get_model(model_name);
    }
    // create the neurites with default param
    for (int i = 0; i < num_neurites; i++)
    {
        if (i + neurites_.size() == 0)
        {
            new_neurite("axon", "axon", axon_gc, rnd_engine);
        }
        else
        {
            std::string name = "dendrite_" + std::to_string(neurites_.size());
            new_neurite(name, "dendrite", dendrite_gc, rnd_engine);
        }
    }
}


void Neuron::initialize_next_event(mtPtr rnd_engine, double new_resolution,
                                   size_t previous_step)
{
    if (use_actin_waves_)
    {
        if (next_actin_event_ == 0)
        {
            next_actin_event(rnd_engine);
        }
        else
        {
            next_actin_event_ =
                (size_t)(next_actin_event_ - previous_step) * new_resolution;
        }
    }
    for (auto &neurite : neurites_)
    {
        neurite.second->branching_model_.initialize_next_event(
            rnd_engine, new_resolution, previous_step);
    }
}


//###################################################
//                  Growth
//###################################################

/**
 * Function that simulates the growth of the neuron.
 */
void Neuron::grow(mtPtr rnd_engine, size_t current_step, double substep)
{
    // if we use actin waves, tell each neurite to update them
    if (use_actin_waves_)
    {
        // new actin wave generation
        // @todo: update this to check substep
        if (current_step == aw_generation_step_)
        {
            next_actin_event(rnd_engine);
            // pick random neurite
            size_t idx_neurite =
                neurites_.size() * uniform_(*(rnd_engine).get());
            auto rnd_neurite = std::next(std::begin(neurites_), idx_neurite);
            rnd_neurite->second->start_actin_wave(actin_content_);
        }
        // actin wave update
        for (const auto &neurite : neurites_)
        {
            neurite.second->update_actin_waves(rnd_engine, substep);
        }
    }

    // For each neurite of the neuron, apply growth
    for (auto &neurite : neurites_)
    {
        neurite.second->grow(rnd_engine, current_step, substep);
    }
}


/**
 * @brief branch according to the event
 */
void Neuron::branch(mtPtr rnd_engine, const branchingEvent& ev)
{
    std::string name_neurite = std::get<3>(ev);
    neurites_[name_neurite]->branching_model_.branching_event(rnd_engine, ev);
}


// Neurite

/*! Function creating a new ``Neurite`` object for the ``Neuron``.
 *
 *  \param name: string - name of the neurite.
 *  \param neurite_type: string, optional (default: "dendrite") -
 *      type of neurite (either "dendrite" or "axon").
 *  \return string name of the created ``Neurite`` object.
 *  the neurite constructor will initialize the ''Growth Cone'' with
 *  the proper model.
 */
std::string Neuron::new_neurite(const std::string &name,
                                const std::string &neurite_type,
                                const GCPtr gc_model, mtPtr rnd_engine)
{
    // update the neurite in the neuron vector
    NeuronWeakPtr my_weak_ptr =
        std::dynamic_pointer_cast<Neuron>(shared_from_this());
    neurites_.insert(
        {name, std::make_shared<Neurite>(name, neurite_type, my_weak_ptr)});

    //#####################################
    // add first growth cone to the neurite
    //#####################################
    // initialize parameters
    double angle = 0;
    Point cone_start_point;

    // a minimal trophism approach: set the neurite on the other side of the
    // neuron
    // then ensure that the first cone it's not created outside the environment.
    if (name == "axon")
    {
        do
        {
            if (not axon_angle_set_)
            {
                angle       = uniform_(*(rnd_engine).get()) * 2 * M_PI;
                axon_angle_ = angle;
            }
            Point position = get_position();
            cone_start_point =
                Point(position.at(0) + details.soma_radius * cos(angle),
                      position.at(1) + details.soma_radius * sin(angle));
        } while (not kernel().space_manager.env_contains(cone_start_point));
        neurites_[name]->init_first_node(soma_, get_position(), name,
                                         details.soma_radius,
                                         details.axon_diameter);
    }
    else
    {
        int counter = 0;
        do
        {
            counter++;
            angle = fmod(axon_angle_ + M_PI, 2 * M_PI) +
                    normal_(*(rnd_engine).get()) * M_PI / 4. * counter;
            Point position = get_position();
            cone_start_point =
                Point(position.at(0) + details.soma_radius * cos(angle),
                      position.at(1) + details.soma_radius * sin(angle));

        } while (not kernel().space_manager.env_contains(cone_start_point));
        neurites_[name]->init_first_node(soma_, get_position(), name,
                                         details.soma_radius,
                                         details.dendrite_diameter);
    }

    // initialize the neurite firstNode, a copy of the soma with child!
    // @todo: clean up this mess!!

    // eventually create the growth cone, cloning the default model.
    neurites_[name]->growth_cone_model_ = growth_cone_model_;
    GCPtr first_gc = gc_model->clone(
        neurites_[name]->get_first_node(), neurites_[name],
        details.soma_radius, name + "0", cone_start_point, angle);
    neurites_[name]->growth_cones_[0] = first_gc;
    first_gc->set_cone_ID();
    first_gc->set_diameter(
        neurites_[name]->get_first_node()->biology_.diameter);

    // then reset the branch with the right initial position
    first_gc->set_first_point(cone_start_point, details.soma_radius);
    // and adjust the topology of the new growth cone
    neurites_[name]->nodes_[0]->children_.push_back(first_gc);
    neurites_[name]->update_tree_structure(neurites_[name]->get_first_node());

#ifndef NDEBUG
    printf("neurite %s is generated with angle %f\n"
           "soma position is %f, %f\n"
           "first growth cone has name: %s \n",
           name.c_str(), angle, get_position().at(0), get_position().at(1),
           neurites_[name]->growth_cones_[0]->get_treeID().c_str());
// and then add the first cone to soma-ghost node
// printf(
//"with parent I get: %s, %lu\n", neurites_[name]
//->growth_cones_[0]
//->get_parent()
//.lock()
//->get_treeID()
//.c_str(),
// neurites_[name]->growth_cones_[0]->get_parent().lock()->get_nodeID());
#endif
    return name;
}


void Neuron::next_actin_event(mtPtr rnd_engine)
{
    // aw_generation_step_ += (size_t)1. / actin_freq_;
    // printf( "next actin event at %lu", aw_generation_step_);
}


//###################################################
//                  Getter/setter functions
//###################################################

void Neuron::set_status(const statusMap &status)
{
    get_param(status, names::soma_radius, details.soma_radius);
    get_param(status, names::axon_diameter, details.axon_diameter);
    get_param(status, names::dendrite_diameter, details.dendrite_diameter);
    //@TODO change the growth cone_model during the growth
    get_param(status, names::growth_cone_model, growth_cone_model_);
    if (get_param(status, names::axon_angle, axon_angle_))
    {
        axon_angle_set_ = true;
    };

    //@TODO implement actin waves
    // auto use_actin_waves_old = use_actin_waves_;
    // get_param(status, names::use_actin_waves, use_actin_waves_);
    // if (use_actin_waves_ != use_actin_waves_old)
    //{
    // next_actin_event_ = 0;
    //}
    // generate time for first actin wave it it wasn't set before
}


void Neuron::set_neurite_status(const std::string &neurite_type,
                                const statusMap &status)
{
    if (neurite_type == "axon")
    {
        neurites_["axon"]->set_status(status);
    }
    else
    {
        for (auto &neurite : neurites_)
        {
            if (neurite.first != "axon")
            {
                neurite.second->set_status(status);
            }
        }
    }
}


/**
 * @brief Get the current value of one of the observables
 */
double Neuron::get_state(const char* observable) const
{
    double value = 0.;

    TRIE(observable)
    CASE("length")
        for (const auto& neurite : neurites_)
        {
            value += neurite.second->get_state(observable);
        }
    CASE("speed")
        for (const auto& neurite : neurites_)
        {
            value += neurite.second->get_state(observable);
        }
    ENDTRIE;

    return value;
}


void Neuron::get_status(statusMap &status) const
{
    set_param(status, names::soma_radius, details.soma_radius);
    set_param(status, names::axon_diameter, details.axon_diameter);
    set_param(status, names::dendrite_diameter, details.dendrite_diameter);
    set_param(status, names::observables, observables_);

    // set position
    Point pos = soma_->get_position();
    set_param(status, "x", pos.at(0));
    set_param(status, "y", pos.at(1));
}


void Neuron::get_neurite_status(statusMap &status, std::string neurite_type)
{

    for (const auto &neurite : neurites_)
    {
        if (neurite.second->neurite_type_ == neurite_type)
        {
            neurite.second->get_status(status);
            break;
        }
    }
}


void Neuron::update_kernel_variables()
{
    for (const auto &neurite : neurites_)
    {
        neurite.second->update_kernel_variables();
    }
}


double Neuron::get_soma_radius() const
{
    return details.soma_radius;
}


BaseNodePtr Neuron::get_soma() const
{
    return soma_;
}


int Neuron::get_num_neurites() const
{
    return neurites_.size();
}


std::string Neuron::get_gc_model() const
{
    return growth_cone_model_;
}


Point Neuron::get_position() const
{
    return soma_->get_position();
}


size_t Neuron::get_gid() const 
{
    return gid_;
}

} // namespace
