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


Property::Property(size_t ul_)
    : data_type(SIZE)
    , ul(ul_)
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
    , ll(v)
{
}


Property::Property(const std::vector<size_t> &v)
    : data_type(VEC_SIZE)
    , uu(v)
{
}


Property::Property(std::string s_)
    : data_type(STRING)
    , s(s_)
{
}


Property::Property(const std::vector<std::string> &v)
    : data_type(VEC_STRING)
    , ss(v)
{
}


Property::Property(const Property &prop)
{
    data_type = prop.data_type;

    switch (data_type)
    {
    case VEC_LONG:
        new (&ll) std::vector<long>(prop.ll);
        break;
    case VEC_SIZE:
        new (&uu) std::vector<size_t>(prop.uu);
        break;
    case VEC_STRING:
        new (&ss) std::vector<std::string>(prop.ss);
        break;
    case INT:
        i = prop.i;
        break;
    case SIZE:
        ul = prop.ul;
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
    // call destructor on old content if required
    switch (data_type)
    {
    case VEC_SIZE:
        uu.~vector<size_t>();
        break;
    case VEC_LONG:
        ll.~vector<long>();
        break;
    case VEC_STRING:
        ss.~vector<std::string>();
        break;
    case STRING:
        s.~basic_string();
        break;
    }

    // switch data_type
    data_type = prop.data_type;

    // init new content
    switch (prop.data_type)
    {
    case VEC_SIZE:
        new (&uu) std::vector<size_t>(prop.uu);
        break;
    case VEC_LONG:
        new (&ll) std::vector<long>(prop.ll);
        break;
    case VEC_STRING:
        new (&ss) std::vector<std::string>(prop.ss);
        break;
    case INT:
        i = prop.i;
        break;
    case SIZE:
        ul = prop.ul;
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
    switch (data_type)
    {
    case VEC_SIZE:
        uu.~vector<size_t>();
        break;
    case VEC_LONG:
        ll.~vector<long>();
        break;
    case VEC_STRING:
        ss.~vector<std::string>();
        break;
    case STRING:
        s.~basic_string();
        break;
    }
}
}
