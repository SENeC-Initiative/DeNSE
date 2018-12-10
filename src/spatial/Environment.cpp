#include "Environment.hpp"


namespace growth
{

Environment::Environment(BMultiPolygonPtr environment)
  : environment_(environment)
{
    for (const auto& polygon : *(environment_.get()))
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


const BMultiPolygonPtr Environment::get_environment() const
{
    return environment_;
}


const BMultiLineString& Environment::get_boundary() const
{
    return boundary_;
}

} // namespace growth
