#coding:utf-8

from CstMdl.physics.unit import PhysicalUnit as PU
from CstMdl.physics.unit import PhysicalMean as PM
from CstMdl.physics.unit import PhysicalQuantity as PQ
from CstMdl.physics.unit import PhysicalConstant as PC
import numpy as np

Patm = PC.Patm.setUnit(PU.kPa)
kg_per_kgDA = PU.PhysicalUnit(PU.U.getMul(),"kg/kgDA")
g_per_kgDA = PU.PhysicalUnit(kg_per_kgDA.getMul()*PU.u.milli,"g/kgDA")
cpa = PQ.HeatCapacity(1.006,PU.kJ_per_kgK)
cpa.mean=PM.PhysicalMean("DryAir HeatCapasity")
cpq = PQ.HeatCapacity(1.805,PU.kJ_per_kgK)
cpq.mean=PM.PhysicalMean("Vapor HeatCapasity")
cpw = PQ.HeatCapacity(4.187,PU.kJ_per_kgK)
cpw.mean=PM.PhysicalMean("Water HeatCapasity")
hfg = PQ.SpecificHeat(2501,PU.kJ_per_kg)
hfg.mean=PM.PhysicalMean("Water evaporation Enthalpy")

class Elevation(PQ.Length):
    mean=PM.PhysicalMean("Elevation")
    unit=PU.m    

class AtompsphericPressure(PQ.Pressure):
    mean=PM.PhysicalMean("AtomsphericPressure")
    unit=PU.kPa
    @classmethod
    def fpatm(cls,h,tdb):
        patm=Patm.getQuantity()*np.power(1-0.0065*h.getQuantity()/tdb.absoluteSetUnit(PU.K).getQuantity(),5.258)
        tdb.absoluteSetUnit(PU.degC)
        return cls(patm)
    @classmethod
    def converth(cls,h,tdb):
        self.quantity=self.getQuantity()*np.power(1-0.0065*h.getQuantity()/tdb.absoluteSetUnit(PU.K).getQuantity(),5.258)
        return self
    
class DryAirDensity(PQ.Dense):
    mean=PM.PhysicalMean("DryAirDensity")
    unit=PU.kg_per_m3
    @classmethod
    def frho(cls,tdb,patm=Patm):
        x=patm.setUnit(PU.atm).getQuantity()
        patm.setUnit(PU.kPa)
        rho=353.25/tdb.absoluteSetUnit(PU.K).getQuantity()
        tdb.absoluteSetUnit(PU.degC)
        return cls(rho)
    
class WetAirDensity(PQ.Dense):
    mean=PM.PhysicalMean("WetAirDensity")
    unit=PU.kg_per_m3
    @classmethod
    def frho(cls,tdb,pw,patm=Patm):
        x=patm.setUnit(PU.atm).getQuantity()
        patm.setUnit(PU.kPa)
        rho=353.25*x/tdb.absoluteSetUnit(PU.K).getQuantity()*(1-0.378*pw.getQuantity()/patm.getQuantity())
        tdb.absoluteSetUnit(PU.degC)
        return cls(rho)
    
class SpecificEnthalpy(PQ.SpecificHeat):
    mean=PM.PhysicalMean("SpecificEnthalpy")
    unit=PU.kJ_per_kg
    @classmethod
    def fhair(cls,sh,lh):
        fhair = sh.getQuantity()+lh.getQuantity()
        return cls(fhair)
class SensibleSpecificEnthalpy(PQ.SpecificHeat):
    mean=PM.PhysicalMean("SensibleSpecificEnthalpy")
    unit=PU.kJ_per_kg
    @classmethod
    def fhash(cls,tdb):
        fhash = cpa.getQuantity()*tdb.getQuantity()
        return cls(fhash)
class LatentSpecificEnthalpy(PQ.SpecificHeat):
    mean=PM.PhysicalMean("LatentSpecificEnthalpy")
    unit=PU.kJ_per_kg
    @classmethod
    def fhalh(cls,tdb,w):
        fhalh = w.getQuantity()*(cpq.getQuantity()*tdb.getQuantity()+hfg.getQuantity())
        return cls(fhalh)   
class SaturationEnthalpy(PQ.SpecificHeat):
    mean=PM.PhysicalMean("SaturationEnthalpy")
    unit=PU.kJ_per_kg
class RelativeHumidity(PQ.Dimentionless):
    mean=PM.PhysicalMean("RelativeHumidity")
    unit=PU.percent
    @classmethod
    def fphi(cls,tdb,w,patm=Patm):
        pw = VaporPartialPressure.fpww(w,patm)
        ps = SaturationVaporPartialPressure.fpws(tdb)
        t=ps.getQuantity()
        if isinstance(t,np.ndarray):
            t[t>0]=100*t[t>0]*pw.getQuantity()[t>0]/t[t>0]
        else:
            if t>0:
                t=0
            else:
                t=100*pw.getQuantity()/ps.getQuantity()
        return cls(t)
class SaturationVaporPartialPressure(PQ.Pressure):
    mean=PM.PhysicalMean("SaturationVaporPartialPressure")
    unit=PU.kPa
    @classmethod
    def fpws(cls,tdb):
        c1 = -5800.2206
        c2 = 1.3914993
        c3 = -0.048640239
        c4 = 0.000041764768
        c5 = -0.000000014452093
        c6 = 6.5459673
        td = tdb.absoluteSetUnit(PU.K).getQuantity()
        tdb.absoluteSetUnit(PU.degC)
        fpws = np.exp(c1 / td + c2 + c3 * td + c4 * td * td + c5 * td * td * td + c6 * np.log(td)) / 1000
        return cls(fpws)
class VaporPartialPressure(PQ.Pressure):
    mean=PM.PhysicalMean("VaporPartialPressure")
    unit=PU.kPa
    @classmethod
    def fpww(cls,w,patm=Patm):
        return cls(patm.getQuantity()*w.getQuantity()/(0.62198+w.getQuantity()))
class DryBulbTemperature(PQ.Temperature):
    mean=PM.PhysicalMean("DryBulbTemperature")
    unit=PU.degC
    @classmethod
    def ftdb(cls,w,ha):
        return cls(ha.getQuantity() - hfg.getQuantity()*w.getQuantity()) / (cpa.getQuantity() + cpq.getQuantity() * w.getQuantity())
class DewPointTemperature(PQ.Temperature):
    mean=PM.PhysicalMean("DewPointTemperature")
    unit=PU.degC
    """
    @classmethod
    def ftdew(cls,w,patm=Patm):
        c0 = 6.54
        c1 = 14.526
        c2 = 0.7389
        c3 = 0.09486
        c4 = 0.4569

        ps = VaporPartialPressure.fpww(w,patm).getQuantity()
        if isinstance(ps,np.ndarray):
            dew=np.zeros(len(ps))
            print dew,[ps<0.00001]
            dew[ps<0.00001]=0
            print dew,[(0.611213>=ps)*(ps>=0.00001)]
            alpha=np.log(ps[(0.611213>=ps)*(ps>=0.00001)])
            dew[(0.611213>=ps)*(ps>=0.00001)]=6.09+alpha*(12.608+alpha*0.4959)
            print dew,[ps>0.611213]
            alpha=np.log(ps[ps>0.611213])
            dew[ps>0.611213]= c0+alpha*(c1+alpha*(c2+alpha*c3))+c4*np.power(alpha,0.1984)
            print dew
            ps=dew
        else:
            if ps<0.00001:
                ps=0
            else:
                alpha=np.log10(ps)
                if ps>0.611213:
                    ps=c0+alpha*(c1+alpha*(c2+alpha*c3))+c4*np.power(alpha,0.1984)
                else:
                    ps=6.09 + alpha * (12.608 + alpha * 0.4959)
        return cls(ps)
    """
    @classmethod
    def ftdew(cls,w,patm=Patm):
        c0 = 6.54
        c1 = 14.526
        c2 = 0.7389
        c3 = 0.09486
        c4 = 0.4569

        ps = VaporPartialPressure.fpww(w,patm).getQuantity()
        alpha=np.log(ps/0.611213)
        if isinstance(ps,np.ndarray):
            tdew=np.zeros(len(ps))
            if1=alpha>=0
            tdew[if1]=13.715*alpha[if1]\
                       +8.4262e-1*np.power(alpha[if1],2)\
                       +1.9048e-2*np.power(alpha[if1],3)\
                       +7.8158e-3*np.power(alpha[if1],4)
            if2=alpha<0
            tdew[if2]=13.7204*alpha[if2]\
                       +7.36631e-1*np.power(alpha[if2],2)\
                       +3.32136e-2*np.power(alpha[if2],3)\
                       +7.78591e-4*np.power(alpha[if2],4)
        else:
             if alpha>=0:
                 tdew=13.715*alpha\
                       +8.4262e-1*np.power(alpha,2)\
                       +1.9048e-2*np.power(alpha,3)\
                       +7.8158e-3*np.power(alpha,4)
             else:
                 tdew=13.7204*alpha\
                       +7.36631e-1*np.power(alpha,2)\
                       +3.32136e-2*np.power(alpha,3)\
                       +7.78591e-4*np.power(alpha,4)
                 
        return cls(tdew)

class SaturationTemperature(PQ.Temperature):
    mean=PM.PhysicalMean("SaturationTemperature")
    unit=PU.degC
class WetBulbTemperature(PQ.Temperature):
    mean=PM.PhysicalMean("WetBulbTemperature")
    unit=PU.degC
class AbsoluteHumidity(PQ.Dimentionless):
    mean=PM.PhysicalMean("AbsoluteHumidity")
    unit=kg_per_kgDA
    @classmethod
    def fwpw(cls,pw,patm=Patm):
        return cls(0.62198*pw.getQuantity()/(patm.getQuantity() - pw.getQuantity()))
    @classmethod
    def fwphi(cls,tdb,phi,patm=Patm):
        ps=SaturationVaporPartialPressure.fpws(tdb).getQuantity()
        pw=VaporPartialPressure(phi.setUnit(PU.U).getQuantity()*ps)
        phi.setUnit(PU.percent)
        return cls.fwpw(pw,patm)
    @classmethod
    def fwtwb(cls,tdb,twb,patm=Patm):
        pstwb = AbsoluteHumidity.fpws(twb)
        ws = AbsoluteHumidity.fwpw(pstwb, patm)
        fwtwb = (ws.getQuantity()*(hfg.getQuantity()+(cpq.getQuantity()-cpw.getQuantity())*twb.getQuantity())\
                 -cpa.getQuantity()*(tdb.getQuantity()-twb.getQuantity())) \
                /(hfg.getQuantity()+cpq.getQuantity()*tdb.getQuantity()-cpw.getQuantity()*twb.getQuantity())
        return cls(fwtwb)
    @classmethod
    def fwha(cls,ha):
        fwha = (ha.getQuantity()-cpa.getQuantity()*tdb.getQuantity()) / (cpq.getQuantity()*tdb.getQuantity()+hfg.getQuantity())
        return cls(fwha)

t=DryBulbTemperature(20)
ps=SaturationVaporPartialPressure.fpws(t)
phi=RelativeHumidity(65)
w=AbsoluteHumidity.fwphi(t,phi)
pw=VaporPartialPressure.fpww(w)
tdew=DewPointTemperature.ftdew(w)
rho=WetAirDensity.frho(t,w)
#del t,ps,phi,w,pw,tdew
    
class AirVolume(PQ.Volume):
    mean=PM.PhysicalMean("AirVolume")
    unit=PU.m3
class AirFlux(PQ.VolumeFlux):
    mean=PM.PhysicalMean("AirFlux")
    unit=PU.m3_per_h
class SensibleAirHeat(PQ.Power):
    mean=PM.PhysicalMean("SensibleAirHeat")
    unit=PU.kW
class LatentAirHeat(PQ.Power):
    mean=PM.PhysicalMean("LatentAirHeat")
    unit=PU.kW
class TotalAirHeat(PQ.Power):
    mean=PM.PhysicalMean("TotalAirHeat")
    unit=PU.kW

    
if __name__=="__main__":
    t1=np.array([31,28])
    t1=DryBulbTemperature(t1)
    print t1.getStr()
    ps1=SaturationVaporPartialPressure.fpws(t1)
    print ps1.getStr()
    phi1=np.array([60,50])
    phi1=RelativeHumidity(phi1)
    print phi1.getStr()
    w1=AbsoluteHumidity.fwphi(t1,phi1)
    print w1.getStr()
    pw1=VaporPartialPressure.fpww(w1)
    print pw1.getStr()
    tdew1=DewPointTemperature.ftdew(w1)
    print tdew1.getStr()
    hash1=SensibleSpecificEnthalpy.fhash(t1)
    print hash1.getStr()
    halh1=LatentSpecificEnthalpy.fhalh(t1,w1)
    print halh1.getStr()
    ha1=SpecificEnthalpy.fhair(hash1,halh1)
    print ha1.getStr()
    rho1=WetAirDensity.frho(t1,pw1)
    print rho1
    Q1=AirFlux(1.25)
    Sh1=SensibleAirHeat.cast(hash1*rho1*Q1)
    print Sh1
    Lh1=LatentAirHeat.cast(halh1*rho1*Q1)
    print Lh1
    Th1=TotalAirHeat.cast(Sh1+Lh1)
    print Th1
    
