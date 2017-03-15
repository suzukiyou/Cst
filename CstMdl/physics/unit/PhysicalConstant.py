#coding:utf-8

from sympy import Rational, pi
import sympy.physics.units as u
from sympy.core import mul as symul
import numpy as np
from CstMdl.physics.unit import numericalunits as nu
from CstMdl.physics.unit import PhysicalUnit as PU
from CstMdl.physics.unit import PhysicalQuantity as PQ
from CstMdl.physics.unit import PhysicalMean as PM

class PhysicalConstant(PQ.PhysicalQuantity):
    def __init__(self,unit):
        self.unit=unit
        self.mean=PM.PhysicalMean(unit.getName())
        self.quantity=np.float32(1)
    def getStr(self):
        return "<"+self.mean.getMean()+">"+str(self.quantity)

c = c_constant = speed_of_light = PhysicalConstant(PU.c_unit)
G0 = G0_constant = gravitational_constant =  PhysicalConstant(PU.G0_unit)
u0 = magnet_constant = PhysicalConstant(PU.u0_unit)
e0 = electric_constant = PhysicalConstant(PU.e0_unit)
Z0 = vacuum_impedance = PhysicalConstant(PU.Z0_unit)
planck = planck_constant = PhysicalConstant(PU.planck_unit)
hbar = PhysicalConstant(PU.hbar_unit)
Patm = atomospher_pressure = PhysicalConstant(PU.atm)

Na = avogadro = avogadro_constant = PhysicalConstant(PU.avogadro_unit)
kB = boltzmann = boltzmann_constant = PhysicalConstant(PU.boltzmann_unit)
gee = gee_constant = PhysicalConstant(PU.gee_unit)
sigmaSB = sigmaSB_constant = PhysicalConstant(PU.sigmaSB_unit)
alphaFS = alphaFS_constant = PhysicalConstant(PU.alphaFS_unit)
Rgas = Rgas_constant = PhysicalConstant(PU.Rgas_unit)
elementary_charge = elementary_charge_constant = PhysicalConstant(PU.elementary_charge_unit)
uBohr = uBor_constant = PhysicalConstant(PU.uBor_unit)
uNuc = uNuc_constant = PhysicalConstant(PU.uNuc_unit)
aBohr = aBohr_constant = PhysicalConstant(PU.aBohr_unit)

me = me_constant = PhysicalConstant(PU.me_unit)
ma = ma_constant = PhysicalConstant(PU.ma_unit)
mn = mn_constant = PhysicalConstant(PU.mn_unit)

Rinf = Rinf_constant = PhysicalConstant(PU.Rinf_unit)
Ry = Ry_constant = PhysicalConstant(PU.Ry_unit)
Eh = Eh_constant = PhysicalConstant(PU.Eh_unit)
ARichardson = PhysicalConstant(PU.ARichardson_unit)
Phi0 = Phi0_constant =PhysicalConstant(PU.Phi0_unit)
KJos = KJos_constant = PhysicalConstant(PU.KJos_unit)
RKlitz = RKlitz_constant =PhysicalConstant(PU.RKlitz_unit)

REarth = PhysicalConstant(PU.REarth_unit)
MEarth = PhysicalConstant(PU.MEarth_unit)
Msolar = PhysicalConstant(PU.MEarth_unit)

def find_constant(constant=None):

    import PhysicalConstant as x
    rv = []
    if constant==None:
        for i in dir(x):
            if isinstance(eval("x."+i),x.PhysicalConstant):
                rv.append(i)
    elif isinstance(constant, str):
        rv = [i for i in dir(x) if constant in i]
    else:
        pass
    return sorted(rv, key=len)

if __name__=="__main__":
    

    m1=PQ.PhysicalQuantity(5.9722e24,PU.kg)
    m2=PQ.PhysicalQuantity(1,PU.kg)
    r1=PQ.PhysicalQuantity(6371,PU.km)
    F1=G0*m1*m2/r1/r1
    print F1.setUnit(PU.N).getStr()
    print find_constant()
    
