/*
 * Swc.cpp
 *
 * This file is part of DeNSE.
 *
 * Copyright (C) 2019 SeNEC Initiative
 *
 * DeNSE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * DeNSE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with DeNSE. If not, see <http://www.gnu.org/licenses/>.
 */

#include "Swc.hpp"

#include "GrowthCone.hpp"
#include "Neuron.hpp"
#include "Node.hpp"


namespace growth
{

typedef std::array<std::vector<double>, 3> PointsArray;

/**
 * @brief Create a file with neuron in SWC format
 *
 * We are using the CNIC format as described in
 * http://www.neuronland.org/NLMorphologyConverter/MorphologyFormats/SWC/Spec.html
 *
 * The constructor will create the SWC class which is used to write all neurons,
 * each neuron is written separatly.
 * @param output_file swc format file
 * @param resolution swc point density respect to DeNSE branch point
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
 * @param stype gid the identificative number of the neuron in the
 * NeuronManager
 */
void Swc::to_swc(const Neuron *neuron, stype gid)
{
    // write the number of the neuron, it can be whichever identifier used to
    // distinguish between neurons.
    swc_file_ << "#start_neuron " << gid << "\n";

    // start the swc format, the first point is the soma
    // and it s written as a single point which is meant as a
    // sphere af the diameter of the neurite diamater field
    // as described in
    // https://www.neuron.yale.edu/phpbb/viewtopic.php?f=13&t=2161/
    // the previous point is set to -1
    stype last_sample = 0;
    stype sample      = 1;
    BPoint pos         = neuron->get_position();
    swc_file_ << 1 << " " << somaID << " " << pos.x() << " " << pos.y()
              << " " << 0 << " " << neuron->get_soma_radius() << " " << -1
              << "\n";
    int neuriteID = -1;
    int ID        = -1;
    stype idxp;
    double final_dts, tap_r, final_diam;
    PointArray pp;
    BranchPtr b; // points and lengthof the branch

    for (const auto &neurite : neuron->neurites_)
    {
        swc_file_ << "#start_neurite " << gid << "." << neurite.first << "\n";
        // attach each neurite to the soma.
        // the neurites are always starting there!
        last_sample = 1;
        tap_r       = neurite.second->get_taper_rate();

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
        std::deque<std::pair<stype, TNodePtr>> nodes{std::make_pair(1, node)};

        while (not nodes.empty())
        {
            idxp               = 0;
            auto node          = nodes.back();
            final_diam         = node.second->get_diameter();
            stype branch_size  = node.second->get_branch()->size();
            b                  = node.second->get_branch();
            final_dts          = b->final_distance_to_soma();
            ID                 = forkID;
            last_sample        = node.first;
            
            for (stype idx = 0; idx < branch_size; idx += resolution_)
            {
                pp = b->at(idx);

                sample++;

                swc_file_ << sample << " " << ID << " " << pp[0] << " " << pp[1]
                          << " " << 0 << " "
                          << 0.5*(final_diam + tap_r*(final_dts - pp[2]))
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
                          << node.second->get_position().x() << " "
                          << node.second->get_position().y() << " " << 0
                          << " " << 0.5*node.second->get_diameter() << " "
                          << last_sample << "\n";
                /*                if (not nodes.empty())*/
                //{
                // last_sample   = nodes.front().first;
                // TNodePtr node = nodes.front().second;
                /*}*/
            }
        }
        swc_file_ << "#end_neurite " << gid << "." << neurite.first << "\n";
    }
    swc_file_ << "#end_neuron " << gid << "\n";
}

Swc::~Swc() { swc_file_.close(); }
} // namespace growth
