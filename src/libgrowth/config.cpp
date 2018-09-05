#include "config.hpp"


namespace growth
{

Property::Property()
    : data_type(INT)
    , dimension("")
    , i(0)
{
}


Property::Property(int i_, const std::string& dim)
    : data_type(INT)
    , dimension(dim)
    , i(i_)
{
}


Property::Property(size_t ul_, const std::string& dim)
    : data_type(SIZE)
    , dimension(dim)
    , ul(ul_)
{
}


Property::Property(bool b_, const std::string& dim)
    : data_type(BOOL)
    , dimension(dim)
    , b(b_)
{
}


Property::Property(char const *arr, const std::string& dim)
    : Property(std::string(arr), dim)
{
}


Property::Property(double d_, const std::string& dim)
    : data_type(DOUBLE)
    , dimension(dim)
    , d(d_)
{
}


Property::Property(const std::vector<long> &v, const std::string& dim)
    : data_type(VEC_LONG)
    , dimension(dim)
    , ll(v)
{
}


Property::Property(const std::vector<size_t> &v, const std::string& dim)
    : data_type(VEC_SIZE)
    , dimension(dim)
    , uu(v)
{
}


Property::Property(std::string s_, const std::string& dim)
    : data_type(STRING)
    , dimension(dim)
    , s(s_)
{
}


Property::Property(const std::vector<std::string> &v, const std::string& dim)
    : data_type(VEC_STRING)
    , dimension(dim)
    , ss(v)
{
}


Property::Property(const std::unordered_map<std::string, double> &v,
                   const std::string& dim)
    : data_type(MAP_DOUBLE)
    , dimension(dim)
    , md(v)
{
}


Property::Property(const Property &prop)
{
    data_type = prop.data_type;
    dimension = prop.dimension;

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
    case MAP_DOUBLE:
        new (&md) std::unordered_map<std::string, double>(prop.md);
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
    case MAP_DOUBLE:
        md.~unordered_map<std::string, double>();
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
    dimension = prop.dimension;

    // init new content
    switch (prop.data_type)
    {
    case VEC_SIZE:
        new (&uu) std::vector<size_t>(prop.uu);
        break;
    case VEC_LONG:
        new (&ll) std::vector<long>(prop.ll);
        break;
    case MAP_DOUBLE:
        new (&md) std::unordered_map<std::string, double>(prop.md);
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
    case MAP_DOUBLE:
        md.~unordered_map<std::string,double>();
        break;
    case VEC_STRING:
        ss.~vector<std::string>();
        break;
    case STRING:
        s.~basic_string();
        break;
    }
}
} // namespace growth
