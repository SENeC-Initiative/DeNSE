#include "exceptions.hpp"

namespace growth
{

InvalidArg::InvalidArg()
    : std::runtime_error("")
    , msg_("")
{
}


InvalidArg::InvalidArg(const std::string &msg, const char* func,
                       const char* file, unsigned int line)
    : std::runtime_error(
        "@" + std::string(func) + " in " + std::string(file) + ":" +
        std::to_string(line) + "\n" + msg)
{
}

//~ const char* InvalidArg::what() const throw()
//~ {
    //~ std::string message = name_ + ": " + msg_;
    //~ return message.c_str();
//~ }


BadPropertyName::BadPropertyName(const std::string &key, const char* func,
                                 const char* file, unsigned int line)
    : InvalidArg(
        "This configuration property does not exist: '" + key + "'.",
        func, file, line)
{
    name_ = "BadPropertyName";
}


BadPropertyType::BadPropertyType(const std::string &key,
                                 const std::string &expected,
                                 const std::string &received,
                                 const char* func, const char* file,
                                 unsigned int line)
    : InvalidArg(
        "Wrong type for configuration property '" + key + "': expected '" +
        expected + "', received '" + received + "'.", func, file, line)
{
    name_ = "BadPropertyType";
}


InvalidParameter::InvalidParameter(const std::string &name,
                                   const std::string &value,
                                   const std::string &condition,
                                   const char* func, const char* file,
                                   unsigned int line)
    : InvalidArg(
        "Invalid value " + value + "for parameter '" + name + "': must be " +
        condition + ".", func, file, line)
{
    name_ = "BadPropertyType";
}

InvalidParameter::InvalidParameter(const std::string &message,
                                   const char* func, const char* file,
                                   unsigned int line)
    : InvalidArg(message, func, file, line)
{
    name_ = "InvalidParameter";
}

InvalidTime::InvalidTime(const char* func, const char* file, unsigned int line)
    : InvalidArg("Negative time obtained.", func, file, line)
{
    name_ = "InvalidParameter";
}

} /* namespace */
