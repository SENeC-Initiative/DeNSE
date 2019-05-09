/*
 * exceptions.hpp
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

#ifndef EXCEPT_H
#define EXCEPT_H

#include <stdexcept>
#include <string>


namespace growth
{

class InvalidArg : public std::invalid_argument
{
  public:
    InvalidArg();

    InvalidArg(const std::string &msg, const char *func, const char *file,
               unsigned int line);

    const char *what() const noexcept override;

    virtual ~InvalidArg() throw() {}

    const char *name() const { return name_.c_str(); }

  protected:
    std::string msg_;
    std::string name_;
};


class BadPropertyName : public InvalidArg
{
  public:
    BadPropertyName(const std::string &key, const char *func, const char *file,
                    unsigned int line);
};


class BadPropertyType : public InvalidArg
{
  public:
    BadPropertyType(const std::string &key, const std::string &expected,
                    const std::string &received, const char *func,
                    const char *file, unsigned int line);
};


class InvalidParameter : public InvalidArg
{
  public:
    InvalidParameter(const std::string &name, const std::string &value,
                     const std::string &condition, const char *func,
                     const char *file, unsigned int line);
    InvalidParameter(const std::string &message, const char *func,
                     const char *file, unsigned int line);
};


class InvalidTime : public InvalidArg
{
  public:
    InvalidTime(const char *func, const char *file, unsigned int line);
};


class LogicError : public std::logic_error
{
  public:
    LogicError();

    LogicError(const std::string &msg, const char *func, const char *file,
               unsigned int line);

    const char *what() const noexcept override;

    virtual ~LogicError() throw() {}

    const char *name() const { return name_.c_str(); }

  protected:
    std::string msg_;
    std::string name_;
};

} // namespace growth

#endif // EXCEPT_H
