#ifndef TIME_H
#define TIME_H

#include <limits>

namespace growth
{

class SimulationManager;
class Neurite;

class Time
{

    friend class SimulationManager;
    friend class Neurite;
    friend class GrowthCone;

  public:
    typedef unsigned long timeStep;

    Time();
    Time(float seconds, unsigned char minutes, unsigned char hours,
         unsigned char days);
    Time(const Time &initial_time, timeStep steps);

    static void reset_resolution();
    static void set_resolution(double resolution);

    static timeStep to_steps(const Time &t);

  private:
    // class members
    static const double DEFAULT_RESOLUTION;
    static double RESOLUTION;
    static const unsigned int MAX_SEC_HMS;
    // instance members
    float sec_;
    unsigned char min_;
    unsigned char hour_;
    unsigned char day_;

  public:
    Time &operator+=(const Time &rhs);
    friend Time operator+(Time lhs, const Time &rhs);
    Time &operator-=(const Time &rhs);
    friend Time operator-(Time lhs, const Time &rhs);

    void update(const unsigned long steps);

    float get_sec() const;
    unsigned char get_min() const;
    unsigned char get_hour() const;
    unsigned char get_day() const;

    void set_sec(float seconds);
    void set_min(unsigned char minutes);
    void set_hour(unsigned char hours);
    void set_day(unsigned char days);
    double get_total_seconds() const;
};

/*
 * Implementation of operators
 */

inline bool operator==(const Time &lhs, const Time &rhs)
{
    bool equal = (lhs.get_day() == rhs.get_day());
    equal &= (lhs.get_hour() == rhs.get_hour());
    equal &= (lhs.get_min() == rhs.get_min());
    equal &= (lhs.get_sec() == rhs.get_sec());
    return equal;
}

inline bool operator!=(const Time &lhs, const Time &rhs)
{
    return !(lhs == rhs);
}

inline bool operator<(const Time &lhs, const Time &rhs)
{
    bool smaller = (lhs.get_day() >= rhs.get_day());
    smaller &= (lhs.get_hour() >= rhs.get_hour());
    smaller &= (lhs.get_min() >= rhs.get_min());
    smaller &= (lhs.get_sec() >= rhs.get_sec());
    return smaller;
}

inline bool operator>(const Time &lhs, const Time &rhs) { return rhs < lhs; }

inline bool operator<=(const Time &lhs, const Time &rhs)
{
    return !(rhs < lhs);
}

inline bool operator>=(const Time &lhs, const Time &rhs)
{
    return !(lhs < rhs);
}
}

#endif // TIME_H
