#!/usr/bin/env python
#-*- coding:utf-8 -*-
#
# This file is part of the PyNCulture project, which aims at providing tools to
# easily generate complex neuronal cultures.
# Copyright (C) 2017 SENeC Initiative
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

'''
Geometry utility functions.
'''


_di_mag = {
    'um': 1e-6,
    'mm': 1e-3,
    'cm': 1e-2,
    'dm': 0.1,
    'm': 1.
}


def conversion_magnitude(source_unit, target_unit):
    '''
    Returns the magnitude necessary to convert from values in `source_unit` to
    values in `target_unit`.
    '''
    return _di_mag[source_unit] / _di_mag[target_unit]
