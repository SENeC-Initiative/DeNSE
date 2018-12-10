#include "Area.hpp"


namespace growth
{

Area::Area(BMultiPolygonPtr area, double height, const std::string &name,
           std::unordered_map<std::string, double> properties)
  : height_(height)
  , name_(name)
  , properties_(properties)
  , shape_(area)
{
    for (const auto& polygon : *(area.get()))
    {
        boundary_.push_back(BLineString());

        auto &ls = boundary_.back();
        ls.insert(ls.end(), polygon.outer().begin(), polygon.outer().end());

        for (const auto& inner : polygon.inners())
        {
            boundary_.push_back(BLineString());

            auto &l = boundary_.back();
            l.insert(l.end(), inner.begin(), inner.end());
        }
    }
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


const BMultiPolygonPtr Area::get_area() const
{
    return shape_;
}


const BMultiLineString& Area::get_boundary() const
{
    return boundary_;
}

} // namespace growth
