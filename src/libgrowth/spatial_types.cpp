#include "spatial_types.hpp"

#include <stdexcept>
#include <string>

namespace growth
{

Point::Point() {}

Point::Point(double x, double y)
    : x_(x)
    , y_(y)
{
}

Point::Point(const Point &pt)
    : x_(pt.at(0))
    , y_(pt.at(1))
{
}

double Point::operator[](const int idx)
{
    if (idx == 0 || idx == -2)
        return x_;
    else if (idx == 1 || idx == -1)
        return y_;
    else
    {
        throw std::runtime_error("Index " + std::to_string(idx) +
                                 "is too"
                                 "large for Point.");
    }
}

double Point::at(const int idx) const
{
    if (idx == 0 || idx == -2)
        return x_;
    else if (idx == 1 || idx == -1)
        return y_;
    else
    {
        throw std::runtime_error("Index " + std::to_string(idx) +
                                 "is too"
                                 "large for Point.");
    }
}
}
