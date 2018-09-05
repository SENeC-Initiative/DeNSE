#!/usr/bin/env python
#-*- coding:utf-8 -*-

""" Helper functions for geometry"""

from .geometry import Shape


def _get_wall_area(area, name, culture, env_buffer, width):
    ''' Return wall area '''
    wall_buffer = Shape([])
    for other_name, other_area in culture.areas.items():
        if other_name != name:
            is_higher = other_area.height > area.height
            if is_higher:
                other_buffer = area.intersection(
                    other_area.exterior.buffer(width))
                wall_buffer  = wall_buffer.union(other_buffer)

    wall_buffer = wall_buffer.union(area.intersection(env_buffer))

    return wall_buffer
