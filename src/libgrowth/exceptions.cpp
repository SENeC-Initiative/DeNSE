/*
 * exceptions.cpp
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

#include "exceptions.hpp"

namespace growth
{

InvalidArg::InvalidArg()
    : std::invalid_argument("")
    , msg_("")
{
}


InvalidArg::InvalidArg(const std::string &msg, const char *func,
                       const char *file, unsigned int line)
    : std::invalid_argument("@" + std::string(func) + " in " +
                            std::string(file) + ":" + std::to_string(line) +
                            "\n" + msg)
{
    name_ = "InvalidArgument";
}


const char *InvalidArg::what() const noexcept
{
    std::string msg = name_ + " " + std::string(std::invalid_argument::what());
    char *char_msg  = new char[msg.size() + 1];
    std::copy(msg.begin(), msg.end(), char_msg);
    char_msg[msg.size()] = '\0'; // don't forget the terminating 0
    return char_msg;
}


BadPropertyName::BadPropertyName(const std::string &key, const char *func,
                                 const char *file, unsigned int line)
    : InvalidArg("This configuration property does not exist: '" + key + "'.",
                 func, file, line)
{
    name_ = "BadPropertyName";
}


BadPropertyType::BadPropertyType(const std::string &key,
                                 const std::string &expected,
                                 const std::string &received, const char *func,
                                 const char *file, unsigned int line)
    : InvalidArg("Wrong type for configuration property '" + key +
                     "': expected '" + expected + "', received '" + received +
                     "'.",
                 func, file, line)
{
    name_ = "BadPropertyType";
}


InvalidParameter::InvalidParameter(const std::string &name,
                                   const std::string &value,
                                   const std::string &condition,
                                   const char *func, const char *file,
                                   unsigned int line)
    : InvalidArg("Invalid value " + value + "for parameter '" + name +
                     "': must be " + condition + ".",
                 func, file, line)
{
    name_ = "BadPropertyType";
}


InvalidParameter::InvalidParameter(const std::string &message, const char *func,
                                   const char *file, unsigned int line)
    : InvalidArg(message, func, file, line)
{
    name_ = "InvalidParameter";
}


InvalidTime::InvalidTime(const char *func, const char *file, unsigned int line)
    : InvalidArg("Negative time obtained.", func, file, line)
{
    name_ = "InvalidTime";
}

/*
 * Logic Error
 */

LogicError::LogicError()
    : std::logic_error("")
    , msg_("")
{
}


LogicError::LogicError(const std::string &msg, const char *func,
                       const char *file, unsigned int line)
    : std::logic_error("@" + std::string(func) + " in " + std::string(file) +
                       ":" + std::to_string(line) + "\n" + msg)
{
    name_ = "LogicError";
}


const char *LogicError::what() const noexcept
{
    std::string msg = name_ + " " + std::string(std::logic_error::what());
    char *char_msg  = new char[msg.size() + 1];
    std::copy(msg.begin(), msg.end(), char_msg);
    char_msg[msg.size()] = '\0'; // don't forget the terminating 0
    return char_msg;
}

} // namespace growth
