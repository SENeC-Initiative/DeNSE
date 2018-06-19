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


void Point::shift(const double module, const double direction)
{
    x_ = x_ + cos(direction) * module;
    y_ = y_ + sin(direction) * module;
}


bool Point::operator==(const Point &other) const
{
    return (x_ == other.x_ and y_ == other.y_);
}


Point Point::operator-(const Point &point) const
{
    double dx = x_ - point.x_;
    double dy = y_ - point.y_;
    return Point(dx, dy);
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
} // namespace growth
