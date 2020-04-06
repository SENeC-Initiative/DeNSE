/*
 * config_impl.hpp
 *
 * This file is part of DeNSE.
 *
 * Copyright (C) 2019 SeNEC Initiative
 *
 * DeNSE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * DeNSE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with DeNSE. If not, see <http://www.gnu.org/licenses/>.
 */

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
               std::vector<stype> &result)
{
    auto res = map.find(key);
    if (res != map.end())
    {
        result = std::vector<stype>(res->second.uu);
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


bool get_param(const statusMap &map, const std::string &key, stype &result)
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
               const std::string &dim)
{
    map[key] = Property(value, dim);
}


void set_param(statusMap &map, const std::string &key, const double &value,
               const std::string &dim)
{
    map[key] = Property(value, dim);
}


void set_param(statusMap &map, const std::string &key, const int &value,
               const std::string &dim)
{
    map[key] = Property(value, dim);
}


void set_param(statusMap &map, const std::string &key, const stype &value,
               const std::string &dim)
{
    map[key] = Property(value, dim);
}


void set_param(statusMap &map, const std::string &key,
               const std::vector<stype> &value, const std::string &dim)
{
    map[key] = Property(value, dim);
}


void set_param(statusMap &map, const std::string &key,
               const std::vector<long> &value, const std::string &dim)
{
    map[key] = Property(value, dim);
}


void set_param(statusMap &map, const std::string &key, const std::string &value,
               const std::string &dim)
{
    map[key] = Property(value, dim);
}


void set_param(statusMap &map, const std::string &key,
               const std::vector<std::string> &value, const std::string &dim)
{
    map[key] = Property(value, dim);
}


void set_param(statusMap &map, const std::string &key,
               const std::unordered_map<std::string, double> &value,
               const std::string &dim)
{
    map[key] = Property(value, dim);
}

} // namespace growth

#endif /* CONFIG_H_IMPL */
