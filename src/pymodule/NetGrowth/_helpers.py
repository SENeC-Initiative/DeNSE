#!/usr/bin/env python
#-*- coding:utf-8 -*-

"""
Helper functions and data.

Contains functions to:
* format time,
* check for container,
* make the hash ids for simulations.

Also declares the dict specifying which levels are valid for each of the
possible observables, to check `levels` at recorder creation.
Same dict is declared to deduce the type of event.
"""

from math import modf
from collections import defaultdict
try:
    from collections.abc import Container as _container
except:
    from collections import Container as _container

import hashlib, json


# ---- #
# Time #
# ---- #

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


# --------- #
# Container #
# --------- #

def nonstring_container(obj):
    '''
    Returns true for any iterable which is not a string or byte sequence.
    '''
    if not isinstance(obj, _container):
        return False
    try:
        if isinstance(obj, unicode):
            return False
    except NameError:
        pass
    if isinstance(obj, bytes):
        return False
    if isinstance(obj, str):
        return False
    return True


# ---------- #
# Hash tools #
# ---------- #

def HashID(*args):
    '''
    Return the hash ID of an experiment.
    '''
    experiment_dict = {}
    for num, dict_ in enumerate(args):
        experiment_dict[num] = dict_
    return _hash_dict(experiment_dict)


def _hash_dict(_dict):
    sha=hashlib.sha1()
    sha.update(str(json.dumps(_dict, sort_keys =True)).encode('utf-8'))
    hash_name = sha.hexdigest()
    return hash_name[:16]


# --------------- #
# Valid arguments #
# --------------- #

valid_levels = {
    "neuron": ["length", "speed", "num_growth_cones", "stopped"],
    "neurite": ["length", "speed", "num_growth_cones", "A", "stopped"],
    "growth_cone": [
        "length", "speed", "resource", "angle", "persistence_angle", "stopped"
    ],
}

# default dict for event types: set discrete events manually
def default_continuous():
    return "continuous"

ev_type = defaultdict(default_continuous)
ev_type["num_growth_cones"] = "discrete"

# unsettable variables used for GetDefaults

unsettables = {
    "growth_cone": ["observables"],
    "neuron": ["observables"],
    "neurite": ["observables"],
}
