#ifndef EXCEPT_H
#define EXCEPT_H

#include <stdexcept>
#include <string>

namespace growth
{

class BadPropertyName : public std::runtime_error
{
  public:
    BadPropertyName(const std::string &key);
};

class BadPropertyType : public std::runtime_error
{
  public:
    BadPropertyType(const std::string &key, const std::string &expected,
                    const std::string &received);
};

class InvalidParameter : public std::runtime_error
{
  public:
    InvalidParameter(const std::string &name, const std::string &value,
                     const std::string &condition);
    InvalidParameter(const std::string &location, const std::string &message);
};

class InvalidTime : public std::runtime_error
{
  public:
    InvalidTime();
};
}

#endif // EXCEPT_H
