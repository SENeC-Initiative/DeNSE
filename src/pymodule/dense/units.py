# -*- coding: utf-8 -*-
#
# units.py
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

import sys
import warnings


import pint
from pint import Quantity, UnitRegistry, set_application_registry

import dense


# hide warnings

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    Quantity([])


# check for the registry

ureg = pint._APP_REGISTRY

if ureg == pint._DEFAULT_REGISTRY:
    ureg = UnitRegistry()
    set_application_registry(ureg)


class FormattedQuantity(ureg.Quantity):

    def _repr_pretty_(self, p, circle):
        p.text(str(self))

    def _repr_latex_(self):
        return '{:L}'.format(self)

    def _repr_html_(self):
        return '{:H}'.format(self)

    def __repr__(self):
        return str(self)

ureg.Quantity = FormattedQuantity
Q_   = ureg.Quantity


# length

m   = ureg.meter
cm  = ureg.cm
mm  = ureg.mm
um  = ureg.micrometer


# time

day    = ureg.day
hour   = ureg.hour
minute = ureg.min
second = ureg.second


# Volume

L  = ureg.L
mL = ureg.mL
uL = ureg.uL
nL = ureg.nL


# frequency

cps = ureg.count / ureg.second
cpm = ureg.count / ureg.minute
cph = ureg.count / ureg.hour
cpd = ureg.count / ureg.day


# concentration

M  = ureg.mol / ureg.L
mM = ureg.millimol / ureg.L
uM = ureg.micromol / ureg.L
nM = ureg.nanomol / ureg.L


# angles

rad = ureg.rad
deg = ureg.deg


# default units

_cpp_units = {
    # string definitions (differenciate dimensionless from angular)
    "[length]": um,
    "[time]": minute,
    "[substance]": ureg.micromol,
    "[degree]": rad,
    "[frequency]": cpm,
    # unit container definition
    um.dimensionality: um,
    minute.dimensionality: minute,
    uM.dimensionality: uM,
    cpm.dimensionality: cpm,
}

dense._units.update({
    # string definitions (differenciate dimensionless from angular)
    "[length]": um,
    "[time]": minute,
    "[substance]": ureg.micromol,
    "[angle]": deg,
    "[frequency]": cpm,
    # unit container definition
    um.dimensionality: um,
    minute.dimensionality: minute,
    uM.dimensionality: uM,
    cpm.dimensionality: cpm,
})
