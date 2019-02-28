#!/usr/bin/env python
#-*- coding:utf-8 -*-

"""
Helper functions and data.

Contains functions to:
* format time,
* deal with spatial environments,
* check for container,
* deal with units,
* check and parse parameters,
* make the hash ids for simulations.

Also declares the dict specifying which levels are valid for each of the
possible observables, to check `levels` at recorder creation.
Same dict is declared to deduce the type of event.
"""

import sys
from math import modf
from collections import defaultdict, Iterable, namedtuple
try:
    from collections.abc import Container as _container
except:
    from collections import Container as _container

try:
    import IPython
    ip_support = True
except:
    ip_support = False

import numpy as np
from shapely.geometry import Point, MultiLineString
from shapely.ops import cascaded_union, linemerge

from .units import ureg, _cpp_units
from .environment import Shape


DICT_IS_ORDERED = sys.version_info >= (3, 6)


# ------------------------ #
# Classes (Time and Model) #
# ------------------------ #

time_units = ("day", "hour", "minute", "second")

class Time(namedtuple("Time", time_units)):

    __slots__ = ()

    def _repr_pretty_(self, p, cycle):
        str_pretty = ""
        for u in time_units:
            val = getattr(self, u)
            if isinstance(val, ureg.Quantity):
                val = val.m
            if val != 0:
                if str_pretty:
                    str_pretty += " "
                str_pretty += "{} {}{}".format(val, u, "s" if val > 1 else "")
        if not str_pretty:
            str_pretty = "0 minute"
        p.text(str_pretty)


# Models

model_blocks = ("extension", "steering", "direction_selection")

class Model(namedtuple("Model", model_blocks)):

    def __str__(self):
        return "_".join(self)

    def __repr__(self):
        return "Model(elongation={el}, steering={st}, direction={dir})".format(
            el=self.elongation_type, st=self.steering_method,
            dir=self.direction_selection
        )

    def _repr_pretty_(self, p, cycle):
        p.begin_group(4, "Model(")
        p.breakable("")
        p.text("extension={},".format(self.extension))
        p.breakable()
        p.text("steering={},".format(self.steering))
        p.breakable()
        p.text("direction_selection={}".format(self.direction_selection))
        p.end_group(4, '')
        p.breakable("")
        p.text(")")
        


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


# ----- #
# Space #
# ----- #

def get_neurite_angles(pos, soma_size, area, max_neurites):
    '''
    Return the mean angles where neurites will start extending.

    This fuction is used to model the growth of neurites on patterned
    substrates.
    '''
    angles = {}

    try:
        p = Point(*pos.m)
    except AttributeError:
        p = Point(pos[0].m, pos[1].m)

    circle = p.buffer(soma_size.m).exterior

    from .environment import plot_shape
    import matplotlib.pyplot as plt

    intersect = area.intersection(circle)

    if isinstance(intersect, MultiLineString):
        intersect = linemerge(intersect.geoms)

    if isinstance(intersect, MultiLineString):
        lines = intersect.geoms
        for i, l in enumerate(lines):
            p_mid = l.interpolate(0.5, normalized=True)
            theta = np.arctan2(p_mid.y - p.y, p_mid.x - p.x)
            if not angles:
                angles["axon"] = theta*ureg.rad
            elif len(angles) < max_neurites:
                angles["dendrite_{}".format(i)] = theta*ureg.rad
            else:
                break
    elif not intersect.is_empty:
        p_mid = intersect.interpolate(0.5, normalized=True)
        theta = np.arctan2(p_mid.y - p.y, p_mid.x - p.x)
        if not angles:
            angles["axon"] = theta*ureg.rad
        elif len(angles) < max_neurites:
            angles["dendrite_{}".format(len(angles))] = theta*ureg.rad

    return angles


def get_area(area_info, culture):
    area = None
    if isinstance(area_info, str):
        area = culture.areas[area_info]
    elif isinstance(area_info, Shape):
        area = area_info
    elif is_iterable(area_info):
        from shapely.ops import cascaded_union
        areas = []
        for a in area_info:
            if isinstance(a, str):
                areas.append(culture.areas[a])
            else:
                areas.append(a)
        area = cascaded_union(areas)
    return area


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


# ----- #
# Units #
# ----- #

old_hook = sys.displayhook


def format_dict(obj, start=None):
    strout  = "{\n"
    max_len = 0
    start   = "  " if start is None else start + "  "

    for k, v in obj.items():
        max_len = max(max_len, len(str(k)))

    for k, v in obj.items():
        wlen = max_len-len(str(k))
        if isinstance(v, str):
            # data    = ["{}'{}'{}: {},\n".format(start, k, ' '*wlen, repr(v))]
            # idx_lim = data[-1].find(":")
            # # cut at long lines
            # while len(data[-1]) > 82:
            #    tmp = data[-1]
            #    # locate space before 80
            #    idx_sp    = tmp[:80].rfind(" ")
            #    idx_paren = tmp[:80].rfind("(")
            strout += "{}'{}'{}: {},\n".format(start, k, ' '*wlen, repr(v))
        elif isinstance(v, dict):
            strout += "{}'{}'{}: {},\n".format(start, k, ' '*wlen,
                                               format_dict(v), start=start)
        else:
            strout += "{}'{}'{}: {},\n".format(start, k, ' '*wlen, str(v))

    strout += "}"

    return strout


def displayhook(obj):
    """Custom displayhook for the exec in default(), which prevents
    assignment of the _ variable in the builtins.
    """
    # reproduce the behavior of the standard displayhook, not printing None
    if isinstance(obj, ureg.Quantity):
        old_hook(str(obj))
    elif isinstance(obj, dict):
        print(format_dict(obj))
    elif obj is not None:
        old_hook(obj)
    return None


def _sorted_for_pprint(items):
    """
    Sort the given items for pretty printing. Since some predictable
    sorting is better than no sorting at all, we sort on the string
    representation if normal sorting fails.
    """
    items = list(items)
    try:
        return sorted(items)
    except Exception:
        try:
            return sorted(items, key=str)
        except Exception:
            return items


def dict_formatter(obj, p, cycle):
    if cycle:
        return p.text('{...}')
    start = '{'
    end   = '}'
    step  = 2
    p.begin_group(step, start)
    keys = obj.keys()
    max_len = 0
    for k, v in obj.items():
        max_len = max(max_len, len(str(k)))
    # if dict isn't large enough to be truncated, sort keys before displaying
    # From Python 3.7, dicts preserve order by definition, so we don't sort.
    if not DICT_IS_ORDERED \
            and not (p.max_seq_length and len(obj) >= p.max_seq_length):
        keys = _sorted_for_pprint(keys)

    if obj:
        p.breakable()
    for idx, key in enumerate(keys):
        if idx:
            p.text(',')
            p.breakable()
        p.begin_group(max_len + step + 2, "")
        p.pretty(key)
        wlen = max_len-len(str(key))
        p.text(' '*wlen + ': ')
        p.pretty(obj[key])
        p.end_group(max_len + step + 2, "")
    if obj:
        p.end_group(step, '')
        p.breakable()
        p.text(end)
    else:
        p.end_group(step, end)


def _set_pprinter_factory(start, end):
    """
    Factory that returns a pprint function useful for sets and frozensets.
    """
    def inner(obj, p, cycle):
        if cycle:
            return p.text(start + '...' + end)
        if len(obj) == 0:
            p.text(start + end)
        else:
            step = 2
            p.begin_group(step, start)
            p.breakable()
            # Like dictionary keys, we will try to sort the items
            # if there aren't too many
            if not (p.max_seq_length and len(obj) >= p.max_seq_length):
                items = _sorted_for_pprint(obj)
            else:
                items = obj
            for idx, x in p._enumerate(items):
                if idx:
                    p.text(',')
                    p.breakable()
                p.pretty(x)
            p.end_group(step, '')
            p.breakable()
            p.text(end)
    return inner


# update the system's display
try:
    if __IPYTHON__ and ip_support:
        ip = IPython.get_ipython()
        formatter = ip.display_formatter.formatters['text/plain']
        formatter.for_type(dict, dict_formatter)
        formatter.for_type(list, _set_pprinter_factory("[", "]"))
        formatter.for_type(tuple, _set_pprinter_factory("( ", ")"))
        formatter.for_type(set, _set_pprinter_factory("{ ", "}"))
    else:
        sys.displayhook = displayhook
except NameError:
        sys.displayhook = displayhook


def to_cppunit(val, valname):
    ''' Convert a value to the correct cppunit '''
    if isinstance(val, ureg.Quantity):
        correct_unit  = 1
        dimdict       = val.dimensionality
        unitcontainer = val.units.dimensionality

        if unitcontainer in _cpp_units:
            correct_unit = _cpp_units[unitcontainer]
        else:
            for name, power in dimdict.items():
                if "name" == "[length]" and power == -3:
                    correct_unit /= L
                else:
                    correct_unit *= _cpp_units[name]**power

        # check for dimensionless
        if correct_unit == 1:
            correct_unit = "dimensionless"

        try:
            return val.to(correct_unit).magnitude
        except Exception as e:
            str_err = "For " + valname + ": " + str(e)
            e.message = str_err
            print(str_err)
            print(correct_unit)
            print(val.dimensionality)
            print(val, val.units, val.units.dimensionality)
            print(_cpp_units.keys())
            raise e
    return val


# ----------------- #
# Parameter parsing #
# ----------------- #

def is_string(value):
    try:
        return isinstance(value, (str, bytes, unicode))
    except:
        return isinstance(value, (str, bytes))


def is_scalar(value):
    value = value.magnitude if isinstance(value, ureg.Quantity) else value

    return (is_string(value)
            or isinstance(value, dict)
            or not isinstance(value, Iterable))


def is_quantity(value):
    return isinstance(value, ureg.Quantity)


def is_iterable(obj):
    if is_quantity(obj):
        return is_iterable(obj.m)
    return isinstance(obj, Iterable) or nonstring_container(obj)


def neuron_param_parser(param, culture, n, on_area=None, rnd_pos=True):
    if culture is None:
        if rnd_pos:
            raise RuntimeError("Cannot seed neurons randomly in space when "
                               "no spatial environment exists. Please set the "
                               "'position' entry in `params`.")
        if "position" not in param:
                raise RuntimeError("`position` entry required in `params` if "
                                   "no `culture` is provided.")
    elif rnd_pos:
        container = culture
        if on_area is not None:
            if is_scalar(on_area):
                container = culture.areas[on_area]
            else:
                areas     = [culture.areas[name] for name in on_area]
                container = cascaded_union(areas)

        sradius = 0.
        if "soma_radius" in param:
            if nonstring_container(param["soma_radius"]):
                sradius = np.max(param["soma_radius"])
            else:
                sradius = param["soma_radius"]

        xy = culture.seed_neurons(
            container=container, neurons=n, soma_radius=sradius,
            return_quantity=True)

        param["position"] = xy

    if "description" not in param:
        param["description"] = "generic_neuron"

    return param


# --------------- #
# Valid arguments #
# --------------- #

valid_levels = {
    "neuron": ["length", "speed", "num_growth_cones"],
    "neurite": ["length", "speed", "num_growth_cones", "A"],
    "growth_cone": [
        "length", "speed", "resource", "angle", "persistence_angle",
        "retraction_time", "status", "stepping_probability",
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
