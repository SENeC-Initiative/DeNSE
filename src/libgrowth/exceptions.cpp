#include "exceptions.hpp"

namespace growth
{

BadPropertyName::BadPropertyName(const std::string &key)
    : std::runtime_error("This configuration property does not exist: '" + key +
                         "'.")
{
}

BadPropertyType::BadPropertyType(const std::string &key,
                                 const std::string &expected,
                                 const std::string &received)
    : std::runtime_error("Wrong type for configuration property '" + key +
                         "': expected '" + expected + "' but received '" +
                         received + "'.")
{
}

InvalidParameter::InvalidParameter(const std::string &name,
                                   const std::string &value,
                                   const std::string &condition)
    : std::runtime_error("Invalid value " + value + "for parameter '" + name +
                         "': must be " + condition)
{
}

InvalidParameter::InvalidParameter(const std::string &location,
                                   const std::string &message)
    : std::runtime_error("In " + location + ". " + message)
{
}

InvalidTime::InvalidTime()
    : std::runtime_error("Negative time obtained.")
{
}
}
