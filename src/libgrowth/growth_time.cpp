#include "growth_time.hpp"

#include <cmath>
#include <cstdlib>
#include <stdio.h>

#include "exceptions.hpp"

namespace growth
{

const double Time::DEFAULT_RESOLUTION(1.);   // in minutes
double Time::RESOLUTION(DEFAULT_RESOLUTION); // in minutes

void Time::reset_resolution() { RESOLUTION = DEFAULT_RESOLUTION; }


void Time::set_resolution(double resolution) { RESOLUTION = resolution; }


Time Time::from_steps(size_t step, double substep)
{
    Time t = Time();
    t.update(step, substep);

    return t;
}


// Time

Time::Time()
    : sec_(0.)
    , min_(0)
    , hour_(0)
    , day_(0)
{
}


Time::Time(double seconds, unsigned char minutes, unsigned char hours,
           size_t days)
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
    update(steps, 0.);
}


void Time::update(Time::timeStep steps, double substeps)
{
    if (steps != 0L or substeps != 0.)
    {
        // prepare the quotient and remainder struct for integer division
        std::ldiv_t dv{};

        // total minutes in the step duration
        double total_min       = steps*RESOLUTION + substeps;
        Time::timeStep int_min = std::floor(total_min);
        double frac_min        = total_min - int_min;
        double frac_sec        = frac_min*60. - (int)(frac_min*60.);

        // seconds
        dv   = std::div(sec_ + frac_min*60., 60L);
        sec_ = (float)(dv.rem + frac_sec);
        // minutes
        dv   = std::div(min_ + int_min + dv.quot, 60L);
        min_ = (char)dv.rem;
        // hours
        dv    = std::div(hour_ + dv.quot, 24L);
        hour_ = (char)dv.rem;
        // day
        day_ += (size_t)dv.quot;
    }
}

// getters

double Time::get_total_seconds() const
{
    return sec_ + min_ * 60 + hour_ * 3600 + day_ * 86400;
}


double Time::get_total_minutes() const
{
    return sec_ / 60. + min_ + hour_ * 60 + day_ * 1440;
}


double Time::get_total_hours() const
{
    return sec_ / 3600. + min_ / 60. + hour_ + day_ * 24;
}


double Time::get_total_days() const
{
    return sec_ / 86400. + min_ / 1440. + hour_ / 24. + day_;
}


double Time::get_sec() const { return sec_; }


unsigned char Time::get_min() const { return min_; }


unsigned char Time::get_hour() const { return hour_; }


size_t Time::get_day() const { return day_; }


// setters

void Time::add_seconds(double seconds)
{
    sec_ += seconds;
    if (sec_ >= 60.)
    {
        int int_part = ((int)std::floor(sec_)) / 60;
        sec_ -= int_part * 60.;
        set_min(min_ + int_part);
    }
}


void Time::set_sec(double seconds)
{
    sec_ = seconds;
    if (seconds >= 60.)
    {
        int int_part = ((int)std::floor(seconds)) / 60;
        sec_ -= int_part * 60.;
        set_min(min_ + int_part);
    }
}


void Time::set_min(unsigned char minutes)
{
    min_ = minutes;
    if (minutes >= 60)
    {
        minutes /= 60;
        min_ -= minutes * 60;
        set_hour(hour_ + minutes);
    }
}


void Time::set_hour(unsigned char hours)
{
    hour_ = hours;
    if (hours >= 24)
    {
        hours /= 24;
        hour_ -= hours * 24;
        set_day(day_ + hours);
    }
}


void Time::set_day(size_t days) { day_ = days; }


// convert time to steps

void Time::to_steps(const Time &t, timeStep &steps, double &substep)
{
    steps = 0L;

    // days are special (could get too large)
    timeStep min_in_days = t.get_day() * 1440;
    steps += (timeStep)std::floor(min_in_days / RESOLUTION);
    double remainder = (double)min_in_days - steps * RESOLUTION;

    // steps in hours, min, sec
    timeStep minutes = t.get_hour() * 60 + t.get_min() + t.get_sec() / 60.;
    remainder += minutes;
    steps += (timeStep)std::floor(remainder / RESOLUTION);

    substep = remainder/RESOLUTION - std::floor(remainder/RESOLUTION);
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
    size_t int_part;
    double dec_part;
    long int signed_tmp;
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
        throw InvalidTime(__FUNCTION__, __FILE__, __LINE__);
    return *this;
}

Time operator-(Time lhs, const Time &rhs)
{
    lhs -= rhs;
    return lhs;
}
} // namespace growth
