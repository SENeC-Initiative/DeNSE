import pint
from pint import UnitRegistry, set_application_registry

import dense


ureg = UnitRegistry()
set_application_registry(ureg)

Q_   = ureg.Quantity


# length

m   = ureg.meter
cm  = ureg.cm
mm  = ureg.mm
mum = ureg.micrometer


# time

minute = ureg.min
hour   = ureg.hour
second = ureg.second


# frequency

cps = ureg.count / ureg.second
cpm = ureg.count / ureg.minute
cph = ureg.count / ureg.hour


# concentration

M   = ureg.mol / ureg.L
mM  = ureg.millimol / ureg.L
muM = ureg.micromol / ureg.L


# angles

rad = ureg.rad
deg = ureg.deg


# default units

dense._cpp_units.update({
    # string definitions (differenciate dimensionless from angular)
    "[length]": mum,
    "[time]": minute,
    "[substance]": muM,
    "[degree]": deg,
    "[frequency]": cpm,
    # unit container definition
    mum.dimensionality: mum,
    minute.dimensionality: minute,
    muM.dimensionality: muM,
    cpm.dimensionality: cpm,
})

dense._units.update({
    # string definitions (differenciate dimensionless from angular)
    "[length]": mum,
    "[time]": minute,
    "[substance]": muM,
    "[angle]": deg,
    "[frequency]": cpm,
    # unit container definition
    mum.dimensionality: mum,
    minute.dimensionality: minute,
    muM.dimensionality: muM,
    cpm.dimensionality: cpm,
})
