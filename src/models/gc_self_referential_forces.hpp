#ifndef GC_SRF_H

#define GC_SRF_H
#include "GrowthCone.hpp"
#include "elements_types.hpp"
namespace growth
{

// Name of the class and inherited element
class GrowthCone_SelfReferentialForces : public virtual GrowthCone
{

    // typedef void (GrowthCone_SelfReferentialForces::*model_function)(void);
    typedef struct SelfAvoidance
    {
        double force;
        double decay;
    } SelfAvoidance;

    typedef struct SomaTropic
    {
        double force;
        double decay;
    } SomaTropic;

    typedef struct Inertial
    {
        double force;
        double decay;
    } Inertial;

  private:
    double gaussian_x, gaussian_y  ;
    SelfAvoidance selfavoidance_;
    SomaTropic somatropic_;
    Inertial inertial_;

  public:
    // public constructor required for models (consider moving to class?)
    GrowthCone_SelfReferentialForces();
    GrowthCone_SelfReferentialForces(const GrowthCone_SelfReferentialForces &);

    virtual GCPtr clone(BaseWeakNodePtr parent, NeuritePtr neurite,
                        double distanceToParent, std::string binaryID,
                        const Point &position, double angle) override;

    void prepare_for_split() override;
    void after_split() override;
    void initialize_SRF();

    virtual Point compute_target_position(
        const std::vector<double> &directions_weights, mtPtr rnd_engine,
        double &substep, double &new_angle) override;

    virtual void set_status(const statusMap &status) override;
    virtual void get_status(statusMap &status) const override;
};
    double generateGaussianNoise(const double& mean, const double &stdDev);
}
#endif
