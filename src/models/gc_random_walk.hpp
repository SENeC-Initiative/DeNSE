#ifndef GC_RW_H
#define GC_RW_H

#include "GrowthCone.hpp"
#include "elements_types.hpp"


namespace growth
{

// Name of the class and inherited element
class GrowthCone_RandomWalk : public virtual GrowthCone
{

    // typedef void (GrowthCone_RandomWalk::*model_function)(void);
    typedef struct Corr_rw
    {
        double tau;
        double previous_gaussian;
        double f_coeff;
        double sqrt_f_coeff;
        double det_delta;
    } Corr_rw;

    typedef struct Memory
    {
        double tau;
        double alpha_temp_sum;
        double alpha_coeff;
        double effective_angle;
    } Memory;

  private:
    bool _persistence_set_;
    double deterministic_angle_;
    double persistence_length_;
    Corr_rw corr_rw_;
    Memory memory_;

  public:
    // required for models
    GrowthCone_RandomWalk();

    GrowthCone_RandomWalk(const GrowthCone_RandomWalk &);

    virtual GCPtr clone(BaseWeakNodePtr parent, NeuritePtr neurite,
                        double distanceToParent, std::string binaryID,
                        const Point &position, double angle) override;

    void prepare_for_split() override;
    void after_split() override;
    void initialize_RW();

    virtual Point compute_new_position(
        const std::vector<double> &directions_weights, mtPtr rnd_engine,
        double substep, double frac, int n, int omp_id) override;

    virtual void set_status(const statusMap &status) override;
    virtual void get_status(statusMap &status) const override;
    virtual double get_state(const char *observable) const override;
};

} /* namespace */

#endif /* GC_RW_H */
