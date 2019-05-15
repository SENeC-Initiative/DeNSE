#include "Swc.hpp"


#include "GrowthCone.hpp"
#include "Neuron.hpp"
#include "Node.hpp"

namespace growth
{

/**
 * @brief Create a file with neuron in SWC format
 *
 * We are using the CNIC format as described in
 * http://www.neuronland.org/NLMorphologyConverter/MorphologyFormats/SWC/Spec.html
 *
 * The constructor will create the SWC class which is used to write all neurons,
 * each neuron is written separatly.
 * @param output_file swc format file
 * @param resolution swc point density respect to netgrowth branch point
 */
Swc::Swc(std::string output_file, unsigned int resolution)
    : resolution_(resolution)
{
    swc_file_.open(output_file);
    swc_file_ << "#SWC file properties\n";
    swc_file_ << "#sample element_ID X Y Z radius parent\n";
}

void Swc::close_file() { swc_file_.close(); }
/**
 * @brief Write a neuron to swc
 *
 * @param Neuron neuron the smart pointer to the required neuron
 * @param size_t gid the identificative number of the neuron in the
 * NeuronManager
 */
void Swc::to_swc(const Neuron *neuron, size_t gid)
{
    // write the number of the neuron, it can be whichever identifier used to
    // distinguish between neurons.
    swc_file_ << "# neuron number " << gid << "\n";

    // start the swc format, the first point is the soma
    // and it s written as a single point which is meant as a
    // sphere af the diameter of the neurite diamater field
    // as described in
    // https://www.neuron.yale.edu/phpbb/viewtopic.php?f=13&t=2161/
    // the previous point is set to -1
    size_t last_sample = 0;
    size_t sample      = 1;
    Point pos          = neuron->get_position();
    swc_file_ << 1 << " " << somaID << " " << pos.at(0) << " " << pos.at(1)
              << " " << 0 << " " << neuron->get_soma_radius() << " " << -1
              << "\n";
    int neuriteID = -1;
    int ID        = -1;
    for (const auto &neurite : neuron->neurites_)
    {
        // attach each neaurite to the soma.
        // the neurites are always starting there!
        last_sample = 1;
        if (neurite.first == "axon")
        {
            neuriteID = axonID;
        }
        else
        {
            neuriteID = b_dendriteID;
        }

        // cycle: for each neurite start from the soma and procede
        // following the branch container.
        // since branching can happen save the branching point and when the
        // present branch is over start again from there.
        TNodePtr node = neurite.second->get_first_node()->get_child(0);
        std::deque<std::pair<size_t, TNodePtr>> nodes{std::make_pair(1, node)};
        while (not nodes.empty())
        {
            auto node          = nodes.back();
            size_t branch_size = node.second->get_branch()->size();
            ID                 = forkID;
            last_sample        = node.first;
            for (size_t idx = 0; idx < branch_size; idx += resolution_)
            {
                sample++;
                swc_file_ << sample << " " << ID << " "
                          << node.second->get_branch()->points[0].at(idx) << " "
                          << node.second->get_branch()->points[1].at(idx) << " "
                          << 0 << " " << node.second->get_diameter() * 0.5
                          << " " << last_sample << "\n";
                ID          = neuriteID;
                last_sample = sample;
            }
            nodes.pop_back();

            // when the neuron branch is at the end: there is no more point to
            // add then verify if it's a leaves or a an internal branch of the
            // tree
            // in that case add childs and start again from another branch.
            if (node.second->has_child())
            {
                NodePtr mynode = std::dynamic_pointer_cast<Node>(node.second);
                nodes.push_front(
                    std::make_pair(last_sample, mynode->get_child(1)));
                nodes.push_front(
                    std::make_pair(last_sample, mynode->get_child(0)));
            }
            else
            {
                // add growth cone at the end!
                sample++;
                swc_file_ << sample << " " << endID << " "
                          << node.second->get_position().at(0) << " "
                          << node.second->get_position().at(1) << " " << 0
                          << " " << 1 << " " << last_sample << "\n";
                /*                if (not nodes.empty())*/
                //{
                // last_sample   = nodes.front().first;
                // TNodePtr node = nodes.front().second;
                /*}*/
            }
        }
    }
}

Swc::~Swc() { swc_file_.close(); }
}
