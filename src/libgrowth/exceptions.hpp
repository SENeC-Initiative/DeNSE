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

    const char* what() const noexcept override;

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

} /* namespace */

#endif // EXCEPT_H
