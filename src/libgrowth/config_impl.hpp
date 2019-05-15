#ifndef CONFIG_H_IMPL
#define CONFIG_H_IMPL

#include <vector>

#include "config.hpp"

namespace growth
{

// implementation get_param (get parameter from user statusMap)

bool get_param(const statusMap &map, const std::string &key, int &result)
{
    auto res = map.find(key);
    if (res != map.end())
    {
        result = res->second.i;
        return true;
    }
    else
    {
        return false;
    }
}

bool get_param(const statusMap &map, const std::string &key, bool &result)
{
    auto res = map.find(key);
    if (res != map.end())
    {
        result = res->second.b;
        return true;
    }
    else
    {
        return false;
    }
}

bool get_param(const statusMap &map, const std::string &key, double &result)
{
    auto res = map.find(key);
    if (res != map.end())
    {
        result = res->second.d;
        return true;
    }
    else
    {
        return false;
    }
}

bool get_param(const statusMap &map, const std::string &key,
               std::vector<long> &result)
{
    auto res = map.find(key);
    if (res != map.end())
    {
        result = std::vector<long>(res->second.l);
        return true;
    }
    else
    {
        return false;
    }
}

bool get_param(const statusMap &map, const std::string &key,
               std::string &result)
{
    auto res = map.find(key);
    if (res != map.end())
    {
        result = std::string(res->second.s);
        return true;
    }
    else
    {
        return false;
    }
}

// implementation set_param (write to statusMap to show it to the user)

void set_param(statusMap &map, const std::string &key, const int &value)
{
    map.insert({key, Property(value)});
}

void set_param(statusMap &map, const std::string &key, const bool &value)
{
    map.insert({key, Property(value)});
}

void set_param(statusMap &map, const std::string &key, const double &value)
{
    map.insert({key, Property(value)});
}

void set_param(statusMap &map, const std::string &key,
               const std::vector<long> &value)
{
    map.insert({key, Property(value)});
}

void set_param(statusMap &map, const std::string &key, const std::string &value)
{
    map.insert({key, Property(value)});
}
}

#endif /* CONFIG_H_IMPL */
