#include "Area.hpp"

// kernel include
#include "kernel_manager.hpp"


namespace growth
{

Area::Area(GEOSGeom area, GEOSContextHandle_t handler, double height,
           const std::string &name,
           std::unordered_map<std::string, double> properties)
    : height_(height)
    , name_(name)
    , properties_(properties)
{
    assert(0 != area);
    assert(GEOSisValid_r(handler, area));

    for (int i = 0; i < kernel().parallelism_manager.get_num_local_threads();
         i++)
    {
        shape_.push_back(GeomPtr(GEOSGeom_clone_r(handler, area)));
        prepared_area_.push_back(GEOSPrepare_r(handler, area));
        const GEOSGeom border = GEOSBoundary_r(handler, area);
        prepared_border_.push_back(GEOSPrepare_r(handler, border));
        assert(prepared_area_[i] != 0);
        assert(prepared_border_[i] != 0);
    }
}


Area::~Area()
{
    for (const GEOSPreparedGeometry *shape : prepared_area_)
    {
        delete shape;
    }
    prepared_area_.clear();

    for (const GEOSPreparedGeometry *border : prepared_border_)
    {
        delete border;
    }
    prepared_border_.clear();
}


/**
 * @brief return the value of a property modifier for the substrate of the area
 *
 * @param name : name of a property
 *
 * @returns 1. if the name is not contained in the `Area.properties_` (will not
 * change the property) or a positive real that will modify the behavior of the
 * growth on this area.
 */
double Area::get_property(const std::string &name) const
{
    auto it = properties_.find(name);
    if (it == properties_.end())
    {
        return 1.;
    }
    return properties_.at(name);
}


double Area::get_height() const { return height_; }


std::string Area::get_name() const { return name_; }


void Area::get_properties(
    std::unordered_map<std::string, double> &properties) const
{
    properties.insert(properties_.cbegin(), properties_.cend());
}


const GEOSPreparedGeometry *Area::get_area(int omp_id) const
{
    return prepared_area_[omp_id];
}


const GEOSPreparedGeometry *Area::get_border(int omp_id) const
{
    return prepared_border_[omp_id];
}


GeomPtr Area::get_shape(int omp_id) const { return shape_.at(omp_id); }

} /* namespace */
