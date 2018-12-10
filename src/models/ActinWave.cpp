#include "ActinWave.hpp"

#include "Neurite.hpp"
#include "Node.hpp"
#include "kernel_manager.hpp"

#include <cassert>
#include <cmath>
#include <memory>


namespace growth
{

ActinWave::ActinWave(TNodePtr target, double actin_content,
                     NeuritePtr ownNeurite)
    : ownNeurite_(ownNeurite)
    , targetNode_(target)
    , distanceToTarget_(targetNode_->get_distance_parent())
    , actin_content_(actin_content)
    , actin_content_tau_(ACTIN_CONTENT_TAU)
    , actin_wave_speed_(ACTIN_WAVE_SPEED)
    , lateral_branching_proba_(exp(-actin_content_ / actin_content_tau_))
{
    initialize_AW_distributions();
}

ActinWave::ActinWave(const ActinWave &copy)
    : ownNeurite_(copy.ownNeurite_)
    , targetNode_(copy.targetNode_)
    , distanceToTarget_(copy.distanceToTarget_)
    , actin_content_(copy.actin_content_)
    , actin_wave_speed_(copy.actin_wave_speed_)
    , lateral_branching_proba_(copy.lateral_branching_proba_)

{
    initialize_AW_distributions();
}

//~ ActinWave::~ActinWave() {}

void ActinWave::step(mtPtr rnd_engine, double substep)
{
    distanceToTarget_ -= actin_wave_speed_ * substep;
    if (distanceToTarget_ < 0)
    {
        if (targetNode_->has_child())
        {
            actin_on_node();
        }
        else
        {
            if (targetNode_->support_AW())
            {
                actin_on_growth_cone();
            }
        }
    }
    else
    {
        assert(0 < lateral_branching_proba_ < 1);
        if (branch_distribution_(*(rnd_engine).get()) >
            lateral_branching_proba_)
        {
            actin_make_branch(rnd_engine);
        }
    }
}
void ActinWave::actin_on_node()
{
    // the first action direct itself to the child 0
    NodePtr myNode    = std::dynamic_pointer_cast<Node>(targetNode_);
    targetNode_       = myNode->get_child(0);
    distanceToTarget_ = myNode->get_distance_parent();
    actin_content_ *= 0.5;
    lateral_branching_proba_ = exp(-actin_content_ / actin_content_tau_);
    ownNeurite_->add_actin(shared_from_this());

    // the second is created copied with half of actin content.
    ActinPtr newWave  = std::make_shared<ActinWave>(*this);
    targetNode_       = myNode->get_child(1);
    distanceToTarget_ = targetNode_->get_distance_parent();
    ownNeurite_->add_actin(newWave);
}

void ActinWave::actin_on_growth_cone() {}
void ActinWave::actin_make_branch(mtPtr rnd_engine)
{
    assert(actin_wave_speed_ > 0);
    assert(actin_content_ >= 0);
    assert(actin_content_tau_ >= 0);
    double new_length = 1.; // @todo: compute from actinContent
    //@TODO
    size_t TODO = 1;
    NodePtr todo_ptr;
    ownNeurite_->lateral_branching(targetNode_, TODO, todo_ptr, rnd_engine);
}

void ActinWave::initialize_AW_distributions()
{
    assert(0 < lateral_branching_proba_ < 1);
    branch_distribution_ = std::uniform_real_distribution<double>(0, 1);
}

void ActinWave::get_geometry(BPoint &xy, double &angle, mtPtr rnd_engine) const
{
    size_t id_x     = targetNode_->get_branch()->size();
    double distance = targetNode_->get_branch()->at(id_x)[2];
    while (distance <= distanceToTarget_)
    {
        id_x--;
        distance = targetNode_->get_branch()->at(id_x)[2];
    }
    BPoint xy_0 = BPoint(targetNode_->get_branch()->at(id_x)[0],
                       targetNode_->get_branch()->at(id_x)[1]);
    BPoint xy_1 = BPoint(targetNode_->get_branch()->at(id_x + 1)[0],
                       targetNode_->get_branch()->at(id_x + 1)[1]);
    double direction;
    if (xy_1.y() - xy_0.y())
    {
        direction = atan((xy_1.x() - xy_0.x()) / (xy_1.y() - xy_0.y()));
    }
    else
    {
        direction = M_PI / 2.;
    }
    xy    = BPoint(targetNode_->get_branch()->at(id_x)[0],
               targetNode_->get_branch()->at(id_x)[1]);
    angle = get_angle(rnd_engine, direction);
}

double get_angle(mtPtr rnd_engine, double direction)
{
    std::bernoulli_distribution d(0.5);
    std::normal_distribution<double> angle(direction, 1);
    if (d(*(rnd_engine).get()))
    {
        return angle(*(rnd_engine).get()) - M_PI / 2.;
    }
    else
    {
        return angle(*(rnd_engine).get()) + M_PI / 2.;
    }
}

void ActinWave::set_status(const statusMap &status)
{
    get_param(status, names::actin_content, actin_content_);
    get_param(status, names::actin_content_tau, actin_content_tau_);
    get_param(status, names::actin_wave_speed, actin_wave_speed_);
    lateral_branching_proba_ = exp(-actin_content_ / actin_content_tau_);
}

void ActinWave::get_status(statusMap &status) const {}
} // namespace growth
