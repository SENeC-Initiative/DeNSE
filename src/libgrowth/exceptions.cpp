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
    : std::invalid_argument(
        "@" + std::string(func) + " in " + std::string(file) + ":" +
        std::to_string(line) + "\n" + msg)
{
    name_ = "InvalidArgument";
}


const char* InvalidArg::what() const noexcept
{
    std::string msg = name_ + " " + std::string(std::invalid_argument::what());
    char * char_msg = new char[msg.size() + 1];
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
    : std::logic_error(
        "@" + std::string(func) + " in " + std::string(file) + ":" +
        std::to_string(line) + "\n" + msg)
{
    name_ = "LogicError";
}


const char* LogicError::what() const noexcept
{
    std::string msg = name_ + " " + std::string(std::logic_error::what());
    char * char_msg = new char[msg.size() + 1];
    std::copy(msg.begin(), msg.end(), char_msg);
    char_msg[msg.size()] = '\0'; // don't forget the terminating 0
    return char_msg;
}

} /* namespace */
