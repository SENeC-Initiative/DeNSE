#include "Skeleton.hpp"

#include "GrowthCone.hpp"
#include "Neuron.hpp"
#include "Node.hpp"

namespace growth
{
// for plotting the neurons and axons
// TO BE EXPLAINED !
Skeleton::Skeleton() {}

Skeleton::Skeleton(const Neuron *neuron, unsigned int resolution)
{
    soma_position        = neuron->get_position();
    soma_radius          = neuron->get_soma_radius();
    axon                 = SkelNeurite();
    dendrites            = SkelNeurite();
    branching_points     = SkelNeurite();
    growth_cones         = SkelNeurite();
    SkelNeurite branches = SkelNeurite();
    //~ #ifndef NDEBUG
    //~ printf(" %lu neurites inside the neuron skeleton \n",
    //~ neuron->get_num_neurites());
    //~ #endif
    for (const auto &neurite : neuron->neurites_)
    {
        if (neurite.first == "axon")
        {
            NodePtr node = neurite.second->get_first_node();
            std::deque<TNodePtr> nodes{node->get_child(0)};
            while (not nodes.empty())
            {
                TNodePtr node = nodes.back();
                BranchPtr b   = node->get_branch();
                size_t i, last(b->size());
                BPoint p;

                for (i=0; i<last; i+=resolution)
                {
                    p = b->xy_at(i);
                    axon.first.push_back(p.x());
                    axon.second.push_back(p.y());
                }
                if (i != last-1)
                {
                    p = b->get_last_xy();
                    axon.first.push_back(p.x());
                    axon.second.push_back(p.y());
                }

                axon.first.push_back(NAN);
                axon.second.push_back(NAN);

                nodes.pop_back();

                if (node->has_child())
                {
                    branching_points.first.push_back(
                        node->get_position().x());
                    branching_points.second.push_back(
                        node->get_position().y());
                    NodePtr mynode = std::dynamic_pointer_cast<Node>(node);
                    nodes.push_front(mynode->get_child(0));
                    nodes.push_front(mynode->get_child(1));
                }
            }

            // if no dendrites, still mark it with NaNs
            if (neuron->neurites_.size() == 1)
            {
                dendrites.first.push_back(NAN);
                dendrites.second.push_back(NAN);
            }
        }
        else
        {
            NodePtr node = neurite.second->get_first_node();
            std::deque<TNodePtr> nodes{node->get_child(0)};
            while (not nodes.empty())
            {

                TNodePtr node = nodes.back();
                BranchPtr b   = node->get_branch();
                size_t i, last(b->size());
                BPoint p;

                for (i=0; i<last; i+=resolution)
                {
                    p = b->xy_at(i);
                    dendrites.first.push_back(p.x());
                    dendrites.second.push_back(p.y());
                }
                if (i != last-1)
                {
                    p = b->get_last_xy();
                    dendrites.first.push_back(p.x());
                    dendrites.second.push_back(p.y());
                }

                dendrites.first.push_back(NAN);
                dendrites.second.push_back(NAN);

                nodes.pop_back();

                if (node->has_child())
                {
                    branching_points.first.push_back(
                        node->get_position().x());
                    branching_points.second.push_back(
                        node->get_position().y());
                    NodePtr mynode = std::dynamic_pointer_cast<Node>(node);
                    nodes.push_front(mynode->get_child(0));
                    nodes.push_front(mynode->get_child(1));
                }
            }

            // if no axon, still mark it with NaNs
            if (neuron->neurites_.size() == 1)
            {
                axon.first.push_back(NAN);
                axon.second.push_back(NAN);
            }
        }

        // anyway get active growth cones
        for (auto it : neurite.second->active_gc_range())
        {
            double x = it.second->get_position().x();
            double y = it.second->get_position().y();
            growth_cones.first.push_back(x);
            growth_cones.second.push_back(y);
        }

        // add NaNs to help separate data from different neurons afterwards
        branching_points.first.push_back(NAN);
        branching_points.second.push_back(NAN);
        growth_cones.first.push_back(NAN);
        growth_cones.second.push_back(NAN);
    }
}

} // namespace growth
