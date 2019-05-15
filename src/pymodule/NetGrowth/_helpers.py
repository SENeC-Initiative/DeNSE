#!/usr/bin/env python
#-*- coding:utf-8 -*-

""" Helper functions """

from math import modf


def format_time(seconds=0., minutes=0, hours=0, days=0):
    '''
    Format a set of intervals values into a correct time.

    Parameters
    ----------
    seconds : float or int, optional (default: 0.)
        Number of seconds.
    minutes : float or int, optional (default: 0)
        Number of minutes.
    hours : float or int, optional (default: 0)
        Number of hours.
    days : float or int, optional (default: 0)
        Number of days.

    Returns
    -------
    (seconds, minutes, hours, days) formatted correctly, with maximal value
    for `days`, then remaining maximal value for `hours`, etc.
    '''
    if seconds < 0. or minutes < 0 or hours < 0 or days < 0:
        raise RuntimeError("All time units should be positive")
    if isinstance(days, float):
        dec, trunc = modf(days)
        days = trunc
        hours += 24 * dec
    if isinstance(hours, float):
        dec, trunc = modf(hours)
        hours = trunc
        minutes += 60 * dec
    if isinstance(minutes, float):
        dec, trunc = modf(minutes)
        minutes = trunc
        seconds += 60 * dec
    if seconds >= 60.:
        add_mins = seconds // 60
        minutes = int(minutes + add_mins)
        seconds -= 60 * add_mins
    if minutes >= 60:
        add_hours = minutes // 60
        hours = int(hours + add_hours)
        minutes -= 60 * add_hours
    if hours >= 24:
        add_days = hours // 24
        days = int(days + add_days)
        hours -= 24 * add_days
    return (seconds, minutes, hours, days)
