/*
 * growth_time.cpp
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


Time Time::from_steps(stype step, double substep)
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
           stype days)
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
        double total_min       = steps * RESOLUTION + substeps;
        Time::timeStep int_min = std::floor(total_min);
        double frac_min        = total_min - int_min;
        double frac_sec        = frac_min * 60. - (int)(frac_min * 60.);

        // seconds
        dv   = std::div(sec_ + frac_min * 60., 60L);
        sec_ = (float)(dv.rem + frac_sec);
        // minutes
        dv   = std::div(min_ + int_min + dv.quot, 60L);
        min_ = (char)dv.rem;
        // hours
        dv    = std::div(hour_ + dv.quot, 24L);
        hour_ = (char)dv.rem;
        // day
        day_ += (stype)dv.quot;
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


stype Time::get_day() const { return day_; }


// setters

void Time::set_sec(double seconds)
{
    sec_ = seconds;
    if (seconds >= 60.)
    {
        int int_part = seconds / 60;
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


void Time::set_day(stype days) { day_ = days; }


// convert time to steps

void Time::to_steps(const Time &t, timeStep &steps, double &substep)
{
    steps   = (timeStep)std::floor(t.get_total_minutes() / RESOLUTION);
    substep = t.get_total_minutes() - steps * RESOLUTION;
}


// operators

Time &Time::operator+=(const Time &rhs)
{
    this->set_sec(sec_ + rhs.get_sec());
    this->set_min(min_ + rhs.get_min());
    this->set_hour(hour_ + rhs.get_hour());
    this->set_day(day_ + rhs.get_day());
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
    stype int_part;
    double dec_part;
    long int signed_tmp;
    int carry_over = 0;

    // subtract and convert
    sec_ -= rhs.get_sec();
    if (sec_ < 0)
    {
        sec_       = sec_ + 60. - rhs.get_sec();
        carry_over = 1;
    }

    signed_tmp = min_ - rhs.get_min() - carry_over; // min
    if (signed_tmp < 0)
    {
        min_       = min_ + 60 - rhs.get_min() - carry_over;
        carry_over = 1;
    }
    else
    {
        min_       = signed_tmp;
        carry_over = 0;
    }

    signed_tmp = hour_ - rhs.get_hour() - carry_over; // hour
    if (signed_tmp < 0)
    {
        hour_      = 24 + hour_ - rhs.get_hour() - carry_over;
        carry_over = 1;
    }
    else
    {
        hour_      = signed_tmp;
        carry_over = 0;
    }

    signed_tmp = day_ - rhs.get_day() - carry_over; // day
    if (signed_tmp < 0)
    {
        throw InvalidTime(__FUNCTION__, __FILE__, __LINE__);
    }

    day_ = signed_tmp;

    return *this;
}

Time operator-(Time lhs, const Time &rhs)
{
    lhs -= rhs;
    return lhs;
}
} // namespace growth
