#coding:utf-8

from sympy import Rational, pi
import sympy.physics.units as u
from sympy.core import mul as symul
import numpy as np
from CstMdl.physics.unit import numericalunits as nu
from CstMdl.physics.unit import PhysicalUnit as PU
from CstMdl.physics.unit import PhysicalMean as PM

class PhysicalQuantity(object):
    mean=PM.PhysicalMean("")
    unit=PU.U
    def __init__(self,num,unit=None):
        print unit,isinstance(unit,PU.PhysicalUnit)
        if isinstance(unit,PU.PhysicalUnit):
            self.unit=unit
        self.quantity=np.float32(num)
    @classmethod
    def cast(cls,q,unit=None):
        if isinstance(unit,PU.PhysicalUnit):
            rq=cls(1,unit)
        else:
            rq=cls(1)
        q.setUnit(rq.unit)
        rq.quantity=q.quantity
        return rq
        
    def getUnit(self):
        return self.unit
    def getMean(self):
        return self.mean
    def getQuantity(self):
        return self.quantity
    def getStr(self):
        return "<"+self.mean.getMean()+":"+self.unit.getName()+">"+str(self.quantity)
    def __repr__(self):
        return self.getStr()
    def setUnit(self,postunit):    
        preunit=self.unit
        postunit=postunit
        t=self.quantity*PU._convert(preunit,postunit)
        self.quantity=t
        self.unit=postunit
        return self
    def absoluteSetUnit(self,postunit): 
        preunit=self.unit
        postunit=postunit
        t=(self.quantity-preunit.getAbsolutePoint())*PU._convert(preunit,postunit)+postunit.getAbsolutePoint()
        self.quantity=t
        self.unit=postunit
        return self
    def __add__(self,pq):
        assert isinstance(pq,PhysicalQuantity)
        assert PU._eqDimention(self.unit,pq.unit)
        t=self.quantity+pq.setUnit(self.unit).quantity
        return self.__class__(t,self.unit)
    def __sub__(self,pq):
        assert isinstance(pq,PhysicalQuantity)
        assert PU._eqDimention(self.unit,pq.unit)
        t=self.quantity-pq.setUnit(self.unit).quantity
        return self.__class__(t,self.unit)
    def __div__(self,pq):
        assert isinstance(pq,PhysicalQuantity)
        u=self.unit.getMul()/pq.unit.getMul()
        t=self.quantity/pq.quantity
        return PhysicalQuantity(t,PU.PhysicalUnit(u,""))
    def __mul__(self,pq):
        assert isinstance(pq,PhysicalQuantity)
        u=self.unit.getMul()*pq.unit.getMul()
        t=self.quantity*pq.quantity
        return PhysicalQuantity(t,PU.PhysicalUnit(u,""))
    def __cmp__(self,pq):
        pqunit=pq.unit
        assert isinstance(pq,PhysicalQuantity)
        assert PU._eqDimention(self.unit,pq.unit)
        t=self.quantity-pq.setUnit(self.unit).quantity
        pq.setUnit(pqunit)
        return t
    def __len__(self):
        if isinstance(self.quantity,np.ndarray):
            return len(self.quantity)
        else:
            raise TypeError("quantity is not array")
    def __getitem__(self,key):
        if isinstance(self.quantity,np.ndarray):
            return self.quantity.__getitem__(key)
        else:
            raise TypeError("quantity is not array")
        
class Dimentionless(PhysicalQuantity):
    mean=PM.PhysicalMean("Dimentionless")
    unit=PU.U
class Length(PhysicalQuantity):
    mean=PM.PhysicalMean("Length")
    unit=PU.m
class Weight(PhysicalQuantity):
    mean=PM.PhysicalMean("Weight")
    unit=PU.kg
class Time(PhysicalQuantity):
    mean=PM.PhysicalMean("Time")
    unit=PU.s
class Temperature(PhysicalQuantity):
    mean=PM.PhysicalMean("Temperature")
    unit=PU.K
        
class Current(PhysicalQuantity):
    mean=PM.PhysicalMean("ElectricCurrent")
    unit=PU.A
class Luminance(PhysicalQuantity):
    mean=PM.PhysicalMean("Luminance")
    unit=PU.cd
class Radian(PhysicalQuantity):
    mean=PM.PhysicalMean("Radian")
    unit=PU.rad
class Stradian(PhysicalQuantity):
    mean=PM.PhysicalMean("Stradian")
    unit=PU.sr
class Amount(PhysicalQuantity):
    mean=PM.PhysicalMean("Amount")
    unit=PU.mol
class Area(PhysicalQuantity):
    mean=PM.PhysicalMean("Area")
    unit=PU.m2
class Volume(PhysicalQuantity):
    mean=PM.PhysicalMean("Volume")
    unit=PU.m3
class SecondMoment(PhysicalQuantity):
    mean=PM.PhysicalMean("SecondMoment")
    unit=PU.m4
class WaveNumber(PhysicalQuantity):
    mean=PM.PhysicalMean("WaveNumber")
    unit=PU.wave_number
class Frequency(PhysicalQuantity):
    mean=PM.PhysicalMean("Frequency")
    unit=PU.Hz
class Acceleration(PhysicalQuantity):
    mean=PM.PhysicalMean("Acceleration")
    unit=PU.m_per_s2
class Force(PhysicalQuantity):
    mean=PM.PhysicalMean("Force")
    unit=PU.N
class Pressure(PhysicalQuantity):
    mean=PM.PhysicalMean("Pressure")
    unit=PU.Pa
class Energy(PhysicalQuantity):
    mean=PM.PhysicalMean("Energy")
    unit=PU.J
class Power(PhysicalQuantity):
    mean=PM.PhysicalMean("Power")
    unit=PU.W
class Momento(PhysicalQuantity):
    mean=PM.PhysicalMean("Moment")
    unit=PU.kgm_per_s
class Velocity(PhysicalQuantity):
    mean=PM.PhysicalMean("Velocity")
    unit=PU.m_per_s
class VolumeFlux(PhysicalQuantity):
    mean=PM.PhysicalMean("VolumeFlux")
    unit=PU.m3_per_s
class MassFlux(PhysicalQuantity):
    mean=PM.PhysicalMean("MassFlux")
    unit=PU.kg_per_s
class LineMass(PhysicalQuantity):
    mean=PM.PhysicalMean("LineMass")
    unit=PU.kg_per_m
class AreaMass(PhysicalQuantity):
    mean=PM.PhysicalMean("AreaMass")
    unit=PU.kg_per_m2
class Dense(PhysicalQuantity):
    mean=PM.PhysicalMean("Dense")
    unit=PU.kg_per_m3
class RotateVelocity(PhysicalQuantity):
    mean=PM.PhysicalMean("RotateVelocity")
    unit=PU.rad_per_s
class RotateMoment(PhysicalQuantity):
    mean=PM.PhysicalMean("RotateAcceleration")
    unit=PU.kgm2
class ThermalExpantion(PhysicalQuantity):
    mean=PM.PhysicalMean("ThermalExpantion")
    unit=PU.per_K
class HeatCapacity(PhysicalQuantity):
    mean=PM.PhysicalMean("HeatCapacity")
    unit=PU.J_per_kgK
class SpecificHeat(PhysicalQuantity):
    mean=PM.PhysicalMean("SpecificHeat")
    unit=PU.J_per_kg
class ThermalConductivity(PhysicalQuantity):
    mean=PM.PhysicalMean("ThermalConductivity")
    unit=PU.W_per_mK
class HeatTransferCoefficient(PhysicalQuantity):
    mean=PM.PhysicalMean("HeatTransferCoefficient")
    unit=PU.W_per_m2K
class Voltage(PhysicalQuantity):
    mean=PM.PhysicalMean("Voltage")
    unit=PU.V
class Capacitance(PhysicalQuantity):
    mean=PM.PhysicalMean("Capacitance")
    unit=PU.C
class Resistance(PhysicalQuantity):
    mean=PM.PhysicalMean("ElectricResistance")
    unit=PU.ohm
class Conductance(PhysicalQuantity):
    mean=PM.PhysicalMean("ElectricConductance")
    unit=PU.S
class Capatitance(PhysicalQuantity):
    mean=PM.PhysicalMean("ElectricCapacitance")
    unit=PU.F
class MagneticFlux(PhysicalQuantity):
    mean=PM.PhysicalMean("MagneticFlux")
    unit=PU.Wb
class MagneticFluxDensity(PhysicalQuantity):
    mean=PM.PhysicalMean("MagneticFluxDensity")
    unit=PU.T
class Inductance(PhysicalQuantity):
    mean=PM.PhysicalMean("MagneticInductance")
    unit=PU.H


if __name__=="__main__":
    v1=Velocity(5)
    print v1.getStr()
    print v1.setUnit(PU.km_per_h).getStr()
    v2=Velocity(3,PU.m_per_s)
    v3=v2+v1
    print v3.getStr()
    v4=v2-v1
    print v4.getStr()

    n1=v2/v1
    print n1.getStr()
    n2=v2*v1
    print n2.getStr()
    print v1<v2
    print v1>v2

    v5=Voltage(220,PU.V)
    p1=Power(0.75,PU.kW)
    phi1=Dimentionless(0.85)
    i1=Current.cast(p1/phi1/v5)
    print i1
