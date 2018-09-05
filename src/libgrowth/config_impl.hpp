#ifndef CONFIG_H_IMPL
#define CONFIG_H_IMPL

#include <vector>

#include "config.hpp"

namespace growth
{

// implementation get_param (get parameter from user statusMap)

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
               std::vector<size_t> &result)
{
    auto res = map.find(key);
    if (res != map.end())
    {
        result = std::vector<size_t>(res->second.uu);
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
        result = std::vector<long>(res->second.ll);
        return true;
    }
    else
    {
        return false;
    }
}


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


bool get_param(const statusMap &map, const std::string &key, size_t &result)
{
    auto res = map.find(key);
    if (res != map.end())
    {
        result = res->second.ul;
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


bool get_param(const statusMap &map, const std::string &key, char *&result)
{
    auto res = map.find(key);
    if (res != map.end())
    {
        std::string s(res->second.s);
        result = const_cast<char *>(s.c_str()); // problem here
        return true;
    }
    else
    {
        return false;
    }
}


bool get_param(const statusMap &map, const std::string &key,
               std::vector<std::string> &result)
{
    auto res = map.find(key);
    if (res != map.end())
    {
        result = std::vector<std::string>(res->second.ss);
        return true;
    }
    else
    {
        return false;
    }
}


bool get_param(const statusMap &map, const std::string &key,
               std::unordered_map<std::string, double> &result)
{
    auto res = map.find(key);
    if (res != map.end())
    {
        result = std::unordered_map<std::string, double>(res->second.md);
        return true;
    }
    else
    {
        return false;
    }
}


// implementation set_param (write to statusMap to show it to the user)

void set_param(statusMap &map, const std::string &key, const bool &value,
               const std::string& dim)
{
    map[key] = Property(value, dim);
}


void set_param(statusMap &map, const std::string &key, const double &value,
               const std::string& dim)
{
    map[key] = Property(value, dim);
}


void set_param(statusMap &map, const std::string &key, const int &value,
               const std::string& dim)
{
    map[key] = Property(value, dim);
}


void set_param(statusMap &map, const std::string &key, const size_t &value,
               const std::string& dim)
{
    map[key] = Property(value, dim);
}


void set_param(statusMap &map, const std::string &key,
               const std::vector<size_t> &value, const std::string& dim)
{
    map[key] = Property(value, dim);
}


void set_param(statusMap &map, const std::string &key,
               const std::vector<long> &value, const std::string& dim)
{
    map[key] = Property(value, dim);
}


void set_param(statusMap &map, const std::string &key, const std::string &value,
               const std::string& dim)
{
    map[key] = Property(value, dim);
}


void set_param(statusMap &map, const std::string &key,
               const std::vector<std::string> &value, const std::string& dim)
{
    map[key] = Property(value, dim);
}


void set_param(statusMap &map, const std::string &key,
               const std::unordered_map<std::string, double> &value,
               const std::string& dim)
{
    map[key] = Property(value, dim);
}

} // namespace growth

#endif /* CONFIG_H_IMPL */
