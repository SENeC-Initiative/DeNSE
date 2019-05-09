# -*- coding: utf-8 -*-
#
# _helpers_geom.py
#
# This file is part of DeNSE.
#
# Copyright (C) 2019 SeNEC Initiative
#
# DeNSE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# DeNSE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with DeNSE. If not, see <http://www.gnu.org/licenses/>.


""" Helper functions for environment"""

from .environment import Shape


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
