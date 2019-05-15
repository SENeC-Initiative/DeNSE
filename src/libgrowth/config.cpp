#include "config.hpp"


namespace growth
{

Property::Property()
    : data_type(INT)
    , i(0)
{
}

Property::Property(int i_)
    : data_type(INT)
    , i(i_)
{
}

Property::Property(bool b_)
    : data_type(BOOL)
    , b(b_)
{
}

Property::Property(char const *arr)
    : Property(std::string(arr))
{
}

Property::Property(double d_)
    : data_type(DOUBLE)
    , d(d_)
{
}

Property::Property(const std::vector<long> &v)
    : data_type(VEC_LONG)
    , l(v)
{
}

Property::Property(std::string s_)
    : data_type(STRING)
    , s(s_)
{
}

Property::Property(const Property &prop)
{
    switch (data_type = prop.data_type)
    {
    case VEC_LONG:
        new (&l) std::vector<long>(prop.l);
        break;
    case INT:
        i = prop.i;
        break;
    case DOUBLE:
        d = prop.d;
        break;
    case BOOL:
        b = prop.b;
        break;
    case STRING:
        new (&s) std::string(prop.s);
        break;
    }
}

Property &Property::operator=(const Property &prop)
{
    if (data_type == VEC_LONG)
        l.~vector<long>();
    else if (data_type == STRING)
        s.~basic_string();
    switch (data_type = prop.data_type)
    {
    case VEC_LONG:
        new (&l) std::vector<long>(prop.l);
        break;
    case INT:
        i = prop.i;
        break;
    case DOUBLE:
        d = prop.d;
        break;
    case BOOL:
        b = prop.b;
        break;
    case STRING:
        new (&s) std::string(prop.s);
        break;
    }
    return *this;
}

Property::~Property()
{
    if (data_type == VEC_LONG)
    {
        l.~vector<long>();
    }
    else if (data_type == STRING)
        s.~basic_string();
}
}
