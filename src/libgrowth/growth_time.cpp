#include "growth_time.hpp"

#include <cmath>
#include <cstdlib>
#include <stdio.h>

#include "exceptions.hpp"

namespace growth
{

const double Time::DEFAULT_RESOLUTION(1.);   // in seconds
double Time::RESOLUTION(DEFAULT_RESOLUTION); // in seconds
const unsigned int Time::MAX_SEC_HMS(86399); // 23 h 59 min 59 s (in seconds)

void Time::reset_resolution() { RESOLUTION = DEFAULT_RESOLUTION; }

void Time::set_resolution(double resolution) { RESOLUTION = resolution; }

// Time

Time::Time()
    : sec_(0.)
    , min_(0)
    , hour_(0)
    , day_(0)
{
}

Time::Time(float seconds, unsigned char minutes, unsigned char hours,
           unsigned char days)
    : sec_(0.)
    , min_(0)
    , hour_(0)
    , day_(0)
{
    set_day(days);
    set_hour(hours);
    set_min(minutes);
    set_sec(seconds);
}

Time::Time(const Time &initial_time, Time::timeStep steps = 0L)
    : sec_(initial_time.get_sec())
    , min_(initial_time.get_min())
    , hour_(initial_time.get_hour())
    , day_(initial_time.get_day())
{
    update(steps);
}

void Time::update(Time::timeStep steps)
{
    if (steps != 0L)
    {
        std::ldiv_t dv{};
        Time::timeStep int_resolution = (Time::timeStep)std::floor(RESOLUTION);
        float fraction_of_seconds     = (RESOLUTION - int_resolution) * steps;
        Time::timeStep total_steps =
            int_resolution * steps +
            (Time::timeStep)std::floor(fraction_of_seconds);
        fraction_of_seconds -= std::floor(fraction_of_seconds);
        // seconds
        dv = std::div(total_steps, 60L);
        sec_ += (float)dv.rem + fraction_of_seconds;
        // minutes
        dv = std::div(dv.quot, 60L);
        min_ += (char)dv.rem;
        // hours
        dv = std::div(dv.quot, 24L);
        hour_ += (char)dv.rem;
        // day
        day_ += (char)dv.quot;
    }
}

// getters

double Time::get_total_seconds() const
{
    return sec_ + min_ * 60 + hour_ * 3600 + day_ * 86400;
}


float Time::get_sec() const { return sec_; }

unsigned char Time::get_min() const { return min_; }

unsigned char Time::get_hour() const { return hour_; }

unsigned char Time::get_day() const { return day_; }

// setters

void Time::set_sec(float seconds)
{
    sec_ = seconds;
    if (seconds >= 60.)
    {
        int int_part = ((int)std::floor(seconds)) / 60;
        min_ += int_part;
        sec_ -= int_part * 60.;
    }
}

void Time::set_min(unsigned char minutes)
{
    min_ = minutes;
    if (minutes >= 60)
    {
        minutes /= 60;
        min_ -= minutes * 60;
        hour_ += minutes;
    }
}

void Time::set_hour(unsigned char hours)
{
    hour_ = hours;
    if (hours >= 24)
    {
        hours /= 24;
        hour_ -= hours * 24;
        day_ += hours;
    }
}

void Time::set_day(unsigned char days) { day_ = days; }

// convert time to steps

Time::timeStep Time::to_steps(const Time &t)
{
    timeStep steps = 0L;

    timeStep L_MAX = std::numeric_limits<unsigned long>::max();
    double DBL_MAX = std::numeric_limits<double>::max();

    // days are special (could get too large)
    timeStep sec_in_days = t.get_day() * 86400;
    steps += (timeStep)std::floor(sec_in_days / RESOLUTION);
    double remainder = (double)sec_in_days - steps * RESOLUTION;

    // steps in hours, min, sec
    timeStep seconds = t.get_hour() * 3600 + t.get_min() * 60 + t.get_sec();
    remainder += seconds / RESOLUTION;
    steps += (timeStep)std::floor(remainder);

    return steps;
}

// operators

Time &Time::operator+=(const Time &rhs)
{
    this->set_day(day_ + rhs.get_day());
    this->set_hour(hour_ + rhs.get_hour());
    this->set_min(min_ + rhs.get_min());
    this->set_sec(sec_ + rhs.get_sec());
    return *this;
}

Time operator+(Time lhs, const Time &rhs)
{
    lhs += rhs;
    return lhs;
}

Time &Time::operator-=(const Time &rhs)
{
    // define integer and decimal parts
    unsigned int int_part;
    float dec_part;
    char signed_tmp;
    // subtract and convert
    signed_tmp = this->day_ - rhs.get_day(); // day
    if (signed_tmp < 0)
    {
        this->day_ = 0;
        signed_tmp = this->hour_ + signed_tmp * 24; // add because negative
    }
    else
    {
        this->day_ = signed_tmp;
        signed_tmp = this->hour_;
    }
    signed_tmp -= rhs.get_hour(); // hour
    if (signed_tmp < 0)
    {
        this->hour_ = 0;
        signed_tmp  = this->min_ + signed_tmp * 60;
    }
    else
    {
        this->hour_ = signed_tmp;
        signed_tmp  = this->min_;
    }
    signed_tmp -= rhs.get_min(); // min
    if (signed_tmp < 0)
    {
        this->min_ = 0;
        this->sec_ += signed_tmp * 60;
    }
    else
    {
        this->min_ = signed_tmp;
    }
    this->sec_ -= rhs.get_sec();
    if (this->sec_ < 0.)
        throw InvalidTime();
    return *this;
}

Time operator-(Time lhs, const Time &rhs)
{
    lhs -= rhs;
    return lhs;
}
}
