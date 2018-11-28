#ifndef MODELS_MANAGER_H
#define MODELS_MANAGER_H

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "GrowthCone.hpp"


namespace growth
{

class GrowthCone;
typedef std::shared_ptr<GrowthCone> GCPtr;


class ModelManager
{
  public:
    const std::vector<std::string> elongation_types;
    const std::vector<std::string> steering_methods;
    const std::vector<std::string> direction_selection;

    const std::unordered_map<std::string, std::string> specials;

    const std::unordered_map<std::string, std::string> abbrev_to_full;
    const std::unordered_map<std::string, std::string> full_to_abbrev;

    const std::string default_model;

  private:
    std::unordered_map<std::string, GCPtr> models_;

  public:
    ModelManager();

    void init_models();

    std::vector<std::string> get_elongation_types();
    std::vector<std::string> get_steering_methods();
    std::vector<std::string> get_direction_selection_methods();

    void get_models(std::unordered_map<std::string, std::string> &models, bool abbrev);
    GCPtr get_model(const std::string& name);
    GCPtr get_default_model();
};

}
#endif /* MODELS_MANAGER_H */
