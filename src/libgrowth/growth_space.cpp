#include "growth_space.hpp"

namespace growth
{

const int Space::DEFAULT_DIM = 2;
int Space::dim_              = Space::DEFAULT_DIM;

Space::Space() {}

int Space::get_dimension() { return dim_; }

void Space::set_dimension_(int dim) { dim_ = dim; }

Space::Shape::Shape() {}
} // namespace growth
