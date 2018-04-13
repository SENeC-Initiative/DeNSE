#include "Neurite.hpp"

// c++ includes
#include <algorithm>
#include <cmath>
#include <sstream>
#define _USE_MATH_DEFINES

// elements includes
#include "GrowthCone.hpp"
#include "Neuron.hpp"

// kernel includes
#include "growth_time.hpp"
#include "kernel_manager.hpp"
#include "search.hpp"

// models includes
#include "ActinWave.hpp"
#include "Node.hpp"
#include "growth_names.hpp"

// debug
#include <assert.h>
#include <typeinfo>


namespace growth
{

/**
 * Neurite constructor:
 * \param neurite_type th etype of the neurite (either "axon" or "dendrite").
 * \param p pointer to parent class (Neuron).
 *
 * The neurite has no points but branches, so the constructor will call a
 * growthConebuild which will instanciate an branch with a random angle.
 */
Neurite::Neurite(std::string name, const std::string &neurite_type,
                 NeuronWeakPtr p)
    : parent_(p)
    , branching_model_(Branching())
    , name_(name)
    , observables_({"length", "speed", "num_growth_cones", "retraction_time", "stopped"})
    , num_created_nodes_(0)
    , num_created_cones_(0)
    , growth_cone_model_("")
    , neurite_type_(neurite_type)
    // parameters for van Pelt branching
    , lateral_branching_angle_mean_(LATERAL_BRANCHING_ANGLE_MEAN)
    , lateral_branching_angle_std_(LATERAL_BRANCHING_ANGLE_STD)
    , gc_split_angle_mean_(GC_SPLIT_ANGLE_MEAN)
    , gc_split_angle_std_(GC_SPLIT_ANGLE_STD)
    , diameter_eta_exp_(DIAMETER_ETA_EXP)
    , diameter_variance_(DIAMETER_VARIANCE)
{
    uniform_ = std::uniform_real_distribution<double>(0., 1.);
    poisson_ = std::poisson_distribution<>(0);
    normal_  = std::normal_distribution<double>(0, 1);
    //~ #ifndef NDEBUG
    //~ printf(" the neurite %s model is %s \n", neurite_type_.c_str(),
    //~ parent_.lock()->get_gc_model().c_str());
    //~ #endif
}


//! default destructor
Neurite::~Neurite()
{
    // it's important to clear actinwaves and growthcones
    // to use properly the SmartPointers
    growth_cones_.clear();
    growth_cones_tmp_.clear();
    dead_cones_.clear();
    actinDeck_.clear();
}


//#######################################################
//                  Initialization
//#######################################################

/**
 * @brief Initialize the first node of the Neurite
 *
 * It's called from neuron during neurite initialization,
 * it creates a new node with same position and attribute of the soma,
 * Now the branching algorithm is standardized for all the nodes
 *
 * @param soma
 * @param pos position of the soma
 * @param neurite_name
 */
void Neurite::init_first_node(BaseWeakNodePtr soma, Point pos, std::string name,
                              double soma_radius, double neurite_diameter)
{
    auto firstNode = std::make_shared<Node>(soma, 0, pos, name);
    firstNode->set_diameter(neurite_diameter);
    firstNode->topology_.has_child         = true;
    firstNode->topology_.centrifugal_order = 0;
    firstNode->geometry_.dis_to_soma       = soma_radius;
    nodes_.insert({num_created_nodes_, firstNode});
    assert(firstNode->get_branch()->size() == 0);
    num_created_nodes_++;

    // also initialize branching model
    if (branching_model_.neurite_ == nullptr)
    {
        branching_model_ = Branching(shared_from_this());
    }
}


void Neurite::finalize()
{
    // include the growth cones that were created during the previous substep
    if (growth_cones_tmp_.size() > 0)
    {
        growth_cones_.insert(growth_cones_tmp_.begin(),
                             growth_cones_tmp_.end());
    }

    growth_cones_tmp_.clear();
}


void Neurite::set_soma_angle(const double angle)
{
    soma_angle_ = angle;
}


double Neurite::get_soma_angle() const
{
    return soma_angle_;
}


/**
 * @brief store kernel parameters locally to avoid repeated calls
 *
 * This functions stores the time resolution from the kernel, it is called
 * every time the user changes that variable in the kernel.
 * It also updates these variables in all active growth cones.
 */
void Neurite::update_kernel_variables()
{
    for (auto &gc : growth_cones_)
    {
        gc.second->update_kernel_variables();
    }
}


//#######################################################
//                  Growth
//#######################################################

/**
 * @brief Neurite grow up with branching and elongation
 *
 * All the dynamical event of the neurite are computed here,
 * at this stage it manages for:
 * - growth cone splitting
 * - lateral branching event
 * - actin wave propagation
 * This function is called once for each simulation step
 *
 * @param rnd_engine
 */
void Neurite::grow(mtPtr rnd_engine, size_t current_step, double substep)
{
    // include the growth cones that were created during the previous substep
    if (growth_cones_tmp_.size() > 0)
    {
        growth_cones_.insert(growth_cones_tmp_.begin(),
                             growth_cones_tmp_.end());
    }

    growth_cones_tmp_.clear();

    // call the branching model specific update
    branching_model_.update_growth_cones(rnd_engine);

    // grow all the growth cones
    for (auto &gc : growth_cones_)
    {
        assert(gc.second.use_count() == 2);
        gc.second->grow(rnd_engine, gc.first, substep);
    }

    // We remove the growth cones afterwards because death occurs while looping
    // on growth_cones_, hence we cannot modify it.
    std::sort(dead_cones_.begin(), dead_cones_.end(), reverse_sorting);
    for (auto &cone_n : dead_cones_)
    {
        // if dead_cones_ is not ordered, then this will fail!
        assert(growth_cones_[cone_n].use_count() == 1);
        growth_cones_.erase(cone_n);
    }

    dead_cones_.clear();
}


/**
 * Remove the parent node of two growth cone, one of which just died.
 * Transfer the remaining child to the grand-parent.
 *
 * @param parent    - the parent of the dead growth cone
 * @param child_id  - the id of the remaining, living growth cone child
 */
void Neurite::delete_parent_node(NodePtr parent, int living_child_id)
{
    auto child = parent->children_[living_child_id];
    size_t grand_parent_ID = parent->get_parent().lock()->get_nodeID();
    NodePtr grand_parent   = nodes_[grand_parent_ID];

#ifndef NDEBUG
    TNodePtr other = parent->children_[1 - living_child_id];
    printf("     size        @@@        size       \n"
           "      %lu         %s         %lu       \n"
           "                  ||                   \n "
           "    (%s)==========(%s)========> %s     \n "
           "     %i            %i           %i     \n \n",
           parent->get_branch()->size(), other->get_treeID().c_str(),
           child->get_branch()->size(), grand_parent->get_treeID().c_str(),
           parent->get_treeID().c_str(), child->get_treeID().c_str(),
           grand_parent->get_centrifugal_order(),
           parent->get_centrifugal_order(), child->get_centrifugal_order());
    printf("positions are:\n"
           "  - (%f; %f)\n  - (%f; %f)\n  - (%f; %f)\n",
           parent->get_position()[0], parent->get_position()[1],
           child->get_position()[0], child->get_position()[1],
           other->get_position()[0], other->get_position()[1]);
#endif
    // reconcile the branch
    parent->biology_.branch->append_branch(child->get_branch());
    child->biology_.branch = parent->biology_.branch;

    // change parental relations
    for (size_t i = 0; i < grand_parent->children_.size(); i++)
    {
        //~ if (grand_parent->get_child(i) == parent) // does this work?
        if (grand_parent->get_child(i)->get_nodeID() == parent->get_nodeID())
        {
            grand_parent->children_[i] = child;
            child->topology_.parent    = grand_parent;
        }
    }

    /*    printf("Old size: %lu ; new_size: %lu. Id was %lu\n", old_nodes_num,*/
    /*nodes_.size(), parent->get_nodeID());*/
    // update tree structure
    update_tree_structure(grand_parent);

#ifndef NDEBUG
    size_t old_nodes_num = nodes_.size();
    printf("           size: %lu                   \n "
           "    (%s)======================> %s     \n "
           "    A                           C      \n"
           "     %i                         %i     \n",
           child->get_branch()->size(), grand_parent->get_treeID().c_str(),
           child->get_treeID().c_str(), grand_parent->get_centrifugal_order(),
           child->get_centrifugal_order());
#endif

    nodes_.erase(parent->get_nodeID());
    assert(nodes_.size() == old_nodes_num - 1);
    //~ dead_nodes_.push_back(parent->get_nodeID());
}


/**
 * Remove a growth cone that was absorbed back inside the parent branch,
 * then call :cpp:func:`delete_parent_node` to remove the lone parent
 * which has only one child.
 */
void Neurite::delete_cone(size_t cone_n)
{
    // delete only if not last growth cone (neurites cannot die)
    if (growth_cones_.size() - dead_cones_.size() > 1)
    {
        GCPtr dead_cone = growth_cones_[cone_n];
#ifndef NDEBUG
        printf(" ############ Cone Deletion #########  \n");
        printf(" dead cone is %s \n", dead_cone->topology_.binaryID.c_str());
#endif
        dead_cone->biology_.dead = true;
#ifndef NDEBUG
        int omp_id = kernel().parallelism_manager.get_thread_local_id();
        printf("got bio on %i\n", omp_id);
        assert(dead_cone->get_branch()->size() == 0);
        printf("asserted\n");
        if (dead_cone->topology_.parent.expired())
        {
            printf("invalid pointer coming\n");
        }
#endif
        size_t parent_ID = dead_cone->get_parent().lock()->get_nodeID();
#ifndef NDEBUG
        printf("got parent ID\n");
#endif
        auto parent_node = nodes_[parent_ID];

#ifndef NDEBUG
        printf("pushing %lu to dead cones\n", cone_n);
#endif
        dead_cones_.push_back(cone_n);

        for (size_t i = 0; i < parent_node->children_.size(); i++)
        {
            if (parent_node->get_child(i)->is_dead() == false)
            {
                delete_parent_node(parent_node, i);
            }
        }

        assert(parent_node.use_count() == 1);
    }
}


//#######################################################
//                  Branching
//#######################################################

/**
 * @brief Compute the angle and the diameter of branching neurite
 * the output variables (old|new)_(diameter|angle) represent the branching cone
 * and the newly created cone
 * both are reference& to the growthcone and are actively used in the
 * computation too.
 * YOU NEED TO PASS THE ACTUAL ANGLE AND DIAMETER TO COMPUTE THE NEXT
 *
 * This function would require a clear model for the branching which relates the
 * angle to the diameter of the neurite.
 * Such a reference it's missing at this moment, excepet for this article: DOI
 * 10.1002/neu.20108
 *
 * @param rnd_engine random enginer
 * @param old_angle
 * @param new_angle
 * @param old_diameter
 * @param new_diameter
 */
void Neurite::gc_split_angles_diameter(mtPtr rnd_engine, double &old_angle,
                                       double &new_angle, double &old_diameter,
                                       double &new_diameter)
{
    // draw the branching angle from a gaussian distribution as hypothized in
    // reference article
    double branching_angle = gc_split_angle_mean_ +
                             gc_split_angle_std_ * normal_(*(rnd_engine).get());

    // ratio between the diameters of the two neurites,
    // it's a gaussian distributed value arround 1.
    double diameter_variance_ = 0.02; //@TODO here the variance was fixed,
                                      // requires to be set by the user!
    double ratio = 1 + diameter_variance_ * normal_(*(rnd_engine).get());
    // The diameters are computed on the base of: d^eta = d_1^eta + d_2^eta
    // looks weird but it's a simple optimization, next two lines compute the
    // new diameter
    // from the old diameter of the neurite and the ratio between the rising
    // cones.
    new_diameter = old_diameter * powf(1. + powf(ratio, diameter_eta_exp_),
                                       -1. / diameter_eta_exp_);
    old_diameter = powf(powf(old_diameter, diameter_eta_exp_) -
                            powf(new_diameter, diameter_eta_exp_),
                        1. / diameter_eta_exp_);

    // the diameter affect the branching angle: the largest cone goes straighter
    // than the other.
    double eps = (old_diameter - new_diameter) / (old_diameter + new_diameter);
    double half_tang = eps * tan(branching_angle / 2);
    new_angle        = old_angle;

    new_angle = -branching_angle / 2. - half_tang;
    old_angle = +branching_angle / 2. - half_tang;

#ifndef NDEBUG
    printf("@@@@@  Growth cone split @@@@@ \n"
           "gc_split angle distribution: %f, +- %f, eta_expo: %f, "
           "diam_variance: %f \n"
           "the branching angle is %f, new_angle: %f, old_angle %f \n"
           "the new diameters are: %f , %f, ratio: %f\n,",
           gc_split_angle_mean_ * 180 / M_PI, gc_split_angle_std_ * 180 / M_PI,
           diameter_eta_exp_, diameter_variance_, branching_angle * 180 / M_PI,
           new_angle * 180 / M_PI, old_angle * 180 / M_PI, new_diameter,
           old_diameter, ratio);
#endif
}


void Neurite::update_parent_nodes(NodePtr new_node, TNodePtr branching)
{
    // neurites_[name]->growth_cones_.back()->set_first_point(pos,soma_radius);
    new_node->topology_.nodeID = num_created_nodes_;

    // update parent node
    NodePtr parent_node = nodes_[new_node->get_parent().lock()->get_nodeID()];
    assert(parent_node->get_nodeID() >= 0);
    assert(parent_node->has_child() == true);

    // I would remove the "first_node" of neurite and just use the Soma,
    // with a vector of TopologicalNode as children and use
    for (size_t i = 0; i < parent_node->children_.size(); i++)
    {
        if (parent_node->get_child(i) == branching)
        {
            parent_node->children_[i] = new_node;
            break;
        }
    }
    branching->topology_.parent   = new_node;
    new_node->topology_.has_child = true;
    nodes_.insert({new_node->get_nodeID(), new_node});
    num_created_nodes_++;
}


GCPtr Neurite::create_branching_cone(const TNodePtr branching_node,
                                     NodePtr new_node, double new_length,
                                     Point xy, double new_cone_angle)
{
    // create a new growth cone
    auto sibling = growth_cones_.begin()->second->clone(
        new_node, shared_from_this(), new_length,
        branching_node->get_treeID() + "1", xy, -3.14);
    statusMap status;
    sibling->set_cone_ID();
    // Here we copy model and status from a random growth cone
    // in the neurite since all the growth cones have same status and
    // model. It's not possible to copy from the splitting cone since this
    // function is common with lateral branching
    // growth_cones_.back()->get_status(status);
    // sibling->set_status(status);
    // insert new elements in the tree
    new_node->children_.push_back(branching_node);
    new_node->children_.push_back(sibling);
    // sibling->topology_.parent       = new_node;
    sibling->biology_.branch = std::make_shared<Branch>(
        xy, new_node->get_branch()->get_last_point()[2]);
    sibling->set_first_point(xy, new_node->get_branch()->get_last_point()[2]);
    add_cone(sibling);
    return sibling;
}


/**
 * @brief Branch from a node of the neuritic tree
 *
 * The branching event can happen wherever along the branch.
 * Both the internal nodes and the leaves (the GrowthCones) can have a
 * lateral branching event.
 * This function will create a new node at the branching event xy and a new
 * growth cone cloning from the 'models_manager'.
 *
 * @param branching_node the node which is going to branch
 * @param branchpoint index of the point in the branch where the
 *        branching occurs.
 * @param new_length the length of the newborn GrowthCone
 * @param rnd_engine
 */
void Neurite::lateral_branching(TNodePtr branching_node, size_t branch_point,
                                double new_length, mtPtr rnd_engine)
{
    // Locate the event and the event parameters
    double branch_direction = 0;
    Point xy;
    locate_from_idx(xy, branch_direction, branching_node->get_branch(),
                    branch_point);
    char branching_side = 1;
    if (uniform_(*(rnd_engine).get()) < 0.5)
    {
        branching_side = -1;
    }
    double angle = branching_side * (lateral_branching_angle_mean_) +
                   lateral_branching_angle_std_ * normal_(*(rnd_engine).get());
    // create a new node and update the node counter
    NodePtr new_node = std::make_shared<Node>(*branching_node);
    update_parent_nodes(new_node, branching_node);
    auto sibling = create_branching_cone(branching_node, new_node, new_length,
                                         xy, branch_direction + angle);
    sibling->move_.angle = branch_direction + angle;

#ifndef NDEBUG
    printf("angle is %f\n", angle * 180 / M_PI);
    printf("angle is %f\n", sibling->move_.angle * 180 / M_PI);
    printf("angle is %f\n", branch_direction * 180 / M_PI);
#endif

    sibling->set_diameter(0.5 * branching_node->get_diameter());
    // lateral branching specific
    branching_node->biology_.branch =
        new_node->biology_.branch->resize_head(branch_point);
    new_node->get_branch()->resize_tail(branch_point + 1);
    // Point xz = Point(10,10);
    // printf("give position %f %f \n", xz.at(0),xz.at(1));
    new_node->set_position(xy);
    // printf("get position %f %f \n", new_node->geometry_.position.at(0),
    // new_node->geometry_.position.at(1));

    // modify sequent nodes
    update_tree_structure(new_node);
#ifndef NDEBUG
    printf("xy: %f, %f, id_x: %lu,  and direction %f \n",
           new_node->get_position().at(0), new_node->get_position().at(1),
           branch_point, branch_direction * 180 / M_PI);
    printf("parent_node has get_treeID: %lu and ID: %s\n",
           new_node->get_parent().lock()->get_nodeID(),
           new_node->get_parent().lock()->get_treeID().c_str());
    printf("biology_.branchsize newcone: %lu               \n"
           " branch size    %s         branch size   \n"
           " oldNode        @@@        newNode/Cone  \n"
           "      %lu       ||             %lu       \n"
           "               || (%.2f)%.2f              \n "
           "    ==========(%s)============>%.2f %s     \n "
           "                                       \n"
           "the centrifugal order of new node is %i \n"
           "the centrifugal order of new cone is %i \n"
           "the centrifugal order of old node is %i \n",
           sibling->get_branch()->size(), sibling->get_treeID().c_str(),
           new_node->get_branch()->size(), branching_node->get_branch()->size(),
           angle * 180 / M_PI, angle * 180 / M_PI + branch_direction * 180 / M_PI,
           new_node->get_treeID().c_str(), branch_direction * 180 / M_PI,
           branching_node->get_treeID().c_str(),
           new_node->get_centrifugal_order(), sibling->get_centrifugal_order(),
           branching_node->get_centrifugal_order());
#endif /* NDEBUG */
    assert(sibling->get_centrifugal_order() ==
           branching_node->get_centrifugal_order());
    assert(new_node->get_child(0) == branching_node);
    assert(new_node->get_child(1) == sibling);
}


/**
 * Manage the "growth cone splitting" branching event, computing the new
 * directions of the two GrowthCone objects.
 *
 * GrowthCone splitting means that an existing GrowthCone will divide into
 * two new GrowthCone objects.
 * This results in a new TopologicalNode being created at the position
 * where the split occured. This new node becomes the new parent of the
 * existing GrowthCone and of the new GrowthCone that is created.
 */
bool Neurite::growth_cone_split(GCPtr branching_cone, double new_length,
                                double new_angle, double old_angle,
                                double new_diameter, double old_diameter)
{
    if (not branching_cone->is_dead())
    {
        auto direction = branching_cone->move_.angle;
#ifndef NDEBUG
        printf("splitting cone %s \n", branching_cone->get_treeID().c_str());
        printf("angles: %f %f, \n", (direction + new_angle) * 180 / 3.13,
               (direction + old_angle) * 180 / 3.14);
#endif
        // branching_cone->split_update();
        // create new node as branching point and add to nodes_ + update counter
        branching_cone->prepare_for_split();
        NodePtr new_node = std::make_shared<Node>(*branching_cone);
        update_parent_nodes(new_node, branching_cone);
        auto sibling = create_branching_cone(
            branching_cone, new_node, new_length,
            branching_cone->geometry_.position, direction + new_angle);

        sibling->move_.angle        = direction + new_angle;
        branching_cone->move_.angle = direction + old_angle;
        // move the old growth cone --> growth cone split specific
        branching_cone->biology_.diameter = old_diameter;
        sibling->biology_.diameter        = new_diameter;
        branching_cone->after_split();
        sibling->after_split();
        branching_cone->topological_advance();
        /*        branching_cone->set_angle(direction+old_angle);*/
        /*sibling->set_angle(direction+new_angle);*/
        branching_cone->biology_.branch = std::make_shared<Branch>(
            branching_cone->get_position(),
            branching_cone->get_branch()->get_distance_to_soma());


#ifndef NDEBUG
        printf(""
               "                      +%.2f ->(%s)    \n"
               " %f                    //                \n"
               " ->(%s)    ==>  (%s)                    \n "
               "                      \\                \n"
               "                      -%.2f ->(%s)    \n"
               "the centrifugal order of new cone is %i \n"
               "the centrifugal order of old cone is %i \n",
               (sibling->move_.angle) * 180 / M_PI,
               branching_cone->get_treeID().c_str(),

               direction * 180 / M_PI, new_node->get_treeID().c_str(),
               new_node->get_treeID().c_str(),

               (branching_cone->move_.angle) * 180 / M_PI,
               sibling->topology_.binaryID.c_str(),

               sibling->get_centrifugal_order(),
               branching_cone->get_centrifugal_order());
#endif /* NDEBUG */
        assert(sibling->get_centrifugal_order() ==
               branching_cone->get_centrifugal_order());
        assert(new_node->get_child(0) == branching_cone);
        assert(nodes_[new_node->get_nodeID()]->has_child() == true);
        return true;
    }
    return false;
}


void Neurite::update_tree_structure(TNodePtr root)
{

    /*    printf("parent is %s  and his centr order is %i \n",*/
    // root->parent_.lock()->get_treeID().c_str(),
    /*root->parent_.lock()->get_centrifugal_order());*/
    std::deque<TNodePtr> nodes{root};
    root->topology_.centrifugal_order =
        root->topology_.parent.lock()->get_centrifugal_order() + 1;
    while (not nodes.empty())
    {
        TNodePtr node = nodes.front();
        nodes.pop_front();
        if (node->has_child())
        {
            NodePtr mynode = std::dynamic_pointer_cast<Node>(node);
            for (size_t i = 0; i < mynode->children_.size(); i++)
            {
                std::stringstream ss;
                ss << i;
                mynode->children_[i]->topology_.binaryID =
                    mynode->topology_.binaryID + ss.str();
                mynode->children_[i]->topology_.centrifugal_order =
                    mynode->get_centrifugal_order() + 1;
                nodes.push_back(mynode->children_[i]);
            }
        }
    }
}


const Branching *Neurite::get_branching_model() const
{
    return &branching_model_;
}


//#######################################################
//              Actin Model
//#######################################################

void Neurite::add_actin(ActinPtr actin) { actinDeck_.push_front(actin); }


void Neurite::update_actin_waves(mtPtr rnd_engine, double substep)
{
    // count the number of actin elements because we only want to loop over
    // them once (they put themselves or new AW back in the deque in step.
    unsigned long aw_count = 0;
    auto size              = actinDeck_.size();
    while (aw_count < size)
    {
        aw_count++;
        auto aw = actinDeck_.front();
        aw->step(rnd_engine, substep);
        actinDeck_.pop_front();
    }
}


void Neurite::start_actin_wave(double actin_content)
{
    actinDeck_.push_back(std::make_shared<ActinWave>(
        get_first_node(), actin_content, shared_from_this()));
}


//#######################################################
//              Utilities Functions
//#######################################################

unsigned int Neurite::num_growth_cones() const
{
    return growth_cones_.size() + growth_cones_tmp_.size();
}


/**
 * @brief Add a GrowthCone to the neurite
 *
 * @param GCPtr pointer to the GrowthCone
 */
void Neurite::add_cone(GCPtr cone)
{
    growth_cones_tmp_[num_created_cones_] = cone;
    num_created_cones_++;
}


std::unordered_map<size_t, GCPtr>::const_iterator Neurite::gc_cbegin() const
{
    return growth_cones_.cbegin();
}


std::unordered_map<size_t, GCPtr>::const_iterator Neurite::gc_cend() const
{
    return growth_cones_.cend();
}


NodePtr Neurite::get_first_node() const { return nodes_.at(0); }


NeuronWeakPtr Neurite::get_parent_neuron() const { return parent_; }


std::string Neurite::get_name() const { return name_; }


size_t Neurite::get_and_increment_gc_ID()
{
    num_created_cones_ += 1;
    return num_created_cones_;
}

//###################################################
//                  get/set status
//###################################################
//

void Neurite::set_status(const statusMap &status)
{
    //~ #ifndef NDEBUG
    //~ printf("\nNEURITE SET STATUS\n");
    //~ #endif
    get_param(status, names::diameter_eta_exp, diameter_eta_exp_);
    get_param(status, names::diameter_variance, diameter_variance_);

    bool is_rad = get_param(status, names::gc_split_angle_mean, gc_split_angle_mean_);
    bool is_rad_std = get_param(status, names::gc_split_angle_std, gc_split_angle_std_);
    //printf("angle %f, std %f \n", gc_split_angle_mean_, gc_split_angle_std_);

    if (is_rad and not kernel().angles_in_radians())
    {
        gc_split_angle_mean_ = _rad_from_deg(gc_split_angle_mean_);
    }

    if (is_rad_std and not kernel().angles_in_radians())
    {
        gc_split_angle_std_ = _rad_from_deg(gc_split_angle_std_);
    }


    is_rad_std = false;
    is_rad = false;
    is_rad = get_param(status, names::lateral_branching_angle_mean, lateral_branching_angle_mean_);
    if ( is_rad and not kernel().angles_in_radians() )
    {
        lateral_branching_angle_mean_ =
            _rad_from_deg(lateral_branching_angle_mean_);
    }

    is_rad = false;
    is_rad = get_param(status, names::lateral_branching_angle_std,lateral_branching_angle_std_);
    if ( is_rad and not kernel().angles_in_radians())
    {
        lateral_branching_angle_std_ =
            _rad_from_deg(lateral_branching_angle_std_);
    }
    is_rad = false;


    if (branching_model_.neurite_ == nullptr)
    {
        branching_model_ = Branching(shared_from_this());
    }
    bool use_cr;
    if (get_param(status, names::use_critical_resource, use_cr))
    {
        if (use_cr)
        {
            observables_.push_back("A"); // add A as observable
        }
        else
        {
            observables_.pop_back(); // delete A as observable (last)
        }
    }
    branching_model_.set_status(status);

    for (auto &gc : growth_cones_)
    {
        gc.second->set_status(status);
    }
}


//@TODO
void Neurite::get_status(statusMap &status, const std::string& level) const
{
    if (level == "neurite")
    {
        set_param(status, names::gc_split_angle_mean,
                  _deg_from_rad(gc_split_angle_mean_));
        set_param(status, names::gc_split_angle_std,
                  _deg_from_rad(gc_split_angle_std_));
        set_param(status, names::lateral_branching_angle_std,
                  _deg_from_rad(lateral_branching_angle_mean_));
        set_param(status, names::lateral_branching_angle_mean,
                  _deg_from_rad(lateral_branching_angle_std_));
        set_param(status, names::observables, observables_);
        // branching properties
        branching_model_.get_status(status);
    }

    // growth_cone properties
    set_param(status, names::growth_cone_model, growth_cone_model_);
    growth_cones_.begin()->second->get_status(status);
}


/**
 * @brief Get the current value of one of the observables
 */
double Neurite::get_state(const char *observable) const
{
    double value = 0.;

    TRIE(observable)
    // default case, just sum up
    for (const auto &gc : growth_cones_)
    {
        value += gc.second->get_state(observable);
    }
    CASE("A")
    value = branching_model_.cr_neurite_.available;
    ENDTRIE;

    return value;
}

double Neurite::get_max_resol() const
{
    // check speed limit
    double max_speed = growth_cones_.begin()->second->avg_speed_ +
                       5*growth_cones_.begin()->second->speed_variance_;
    double length    = growth_cones_.begin()->second->filopodia_.finger_length;

    // check angle limit for simple random walk (sigma depends on resolution)
    if (growth_cone_model_ == "simple_random_walk")
    {
        double sigma_sqrd   = growth_cones_.begin()->second->sensing_angle_ *
                              growth_cones_.begin()->second->sensing_angle_;
        double dt_angle_max = M_PI*M_PI/(9*sigma_sqrd);
        return std::min(length / max_speed, dt_angle_max);
    }

    return length / max_speed;
}

} // namespace

