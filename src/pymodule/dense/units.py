import sys
import pint
from pint import UnitRegistry, set_application_registry

import dense


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
