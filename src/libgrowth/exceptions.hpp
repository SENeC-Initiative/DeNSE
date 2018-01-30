#ifndef EXCEPT_H
#define EXCEPT_H

#include <stdexcept>
#include <string>

namespace growth
{

class InvalidArg : public std::runtime_error
{
  public:
    InvalidArg();

    InvalidArg(const std::string &msg, const char *func, const char *file,
               unsigned int line);

    //~ virtual const char* what() const throw()
    //~ {
    //~ std::string message = name_ + ": " + msg_;
    //~ return message.c_str();
    //~ }

    const char *what()
    {
        std::string message = name_ + ": " + msg_;
        return message.c_str();
    }

    virtual ~InvalidArg() throw() {}

    //~ virtual const char* what() const throw();

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


//~ struct BadPropertyName : std::exception
//~ {
//~ BadPropertyName(const std::string &key, const char* func,
//~ const char* file, unsigned int line)
//~ {
//~ msg_ = "@" + std::string(func) + " in " + std::string(file) + ":" +
//~ std::to_string(line) + ": This configuration property does not "
//~ "exist: '" + key + "'.";
//~ }

//~ const char* what() const noexcept
//~ {
//~ return msg_.c_str();
//~ }

//~ private:
//~ std::string msg_;
//~ };


class BadPropertyType : public InvalidArg
{
  public:
    BadPropertyType(const std::string &key, const std::string &expected,
                    const std::string &received, const char *func,
                    const char *file, unsigned int line);
};


//~ struct BadPropertyType : std::exception
//~ {
//~ BadPropertyType(const std::string &key, const std::string &expected,
//~ const std::string &received, const char* func,
//~ const char* file, unsigned int line)
//~ {
//~ msg_ = "@" + std::string(func) + " in " + std::string(file) + ":" +
//~ std::to_string(line) + ": Wrong type for configuration "
//~ "property '" + key + "': expected '" + expected + "', received "
//~ "'" + received + "'.";
//~ }

//~ const char* what() const noexcept
//~ {
//~ return msg_.c_str();
//~ }

//~ private:
//~ std::string msg_;
//~ };


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

} /* namespace */

#endif // EXCEPT_H
