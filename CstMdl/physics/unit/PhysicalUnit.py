#coding:utf-8

from sympy import Rational, pi
import sympy.physics.units as u
from sympy.core import mul as symul
import numpy as np
from CstMdl.physics.unit import numericalunits as nu

class PhysicalUnit(object):
    _Mul=None
    _name=None
    def __init__(self,unit,name):
        assert isinstance(unit,u.Unit) or (hasattr(unit,"is_real") and unit.is_real) or isinstance(unit,symul.Mul)
        assert isinstance(name,basestring)
        self._name=name
        self._Mul=unit
    def getName(self):
        return self._name
    def getMul(self):
        return self._Mul
    def __repr__(self):
        return "<Physical Unit>:"+self._name
    def setAbsolutePoint(self,num):
        self._absolutePoint=num
    def getAbsolutePoint(self):
        return self._absolutePoint

def _eqDimention(preunit,postunit):
    n=(preunit.getMul()/postunit.getMul())
    try:
        float(n)
        return True
    except:
        return False
    
def _convert(preunit,postunit):
    n=(preunit.getMul()/postunit.getMul())
    return float(n)    


#dimentionless unit
U=PhysicalUnit(Rational(1),"U")
ppc = percent = percents = PhysicalUnit(U.getMul()*u.centi,"%")
ppk = permilli = parts_per_kilo = PhysicalUnit(U.getMul()*u.milli,"ppk")
ppm = parts_per_million = PhysicalUnit(U.getMul()*u.micro,"ppm")

# Base units (sympy.physics.unit)
length = m = meter = meters =PhysicalUnit(u.m,"m")
mass = kg = kilogram = kilograms =PhysicalUnit(u.kg,"kg")
time = s = second = seconds =PhysicalUnit(u.s,"s")
current = A = ampere = amperes = PhysicalUnit(u.A,"A")
temperature = K = kelvin = kelvins = PhysicalUnit(u.K,"K")
K.setAbsolutePoint(0)
amount = mol = mole = moles = PhysicalUnit(u.mol,"mol")
luminosity = cd = candela = candelas = PhysicalUnit(u.cd,"cd")

rad = radian = radians = PhysicalUnit(U.getMul(),"rad")
rot = rotation = PhysicalUnit(rad.getMul()*2*pi,"rot")
deg = degree = degrees = PhysicalUnit(pi/180*U.getMul(),"deg")
sr = steradian = steradians = PhysicalUnit(U.getMul(),"sr")

# Derived units
frequency = Hz = hz = hertz =  PhysicalUnit(1/u.s,"Hz")
force = N = newton = newtons =PhysicalUnit(u.kg*u.m/u.s**2,"N")
energy = J = joule = joules = PhysicalUnit(N.getMul()*u.m,"J")
power = W = watt = watts = PhysicalUnit(J.getMul()/u.s,"W")
pressure = Pa = pa = pascal = pascals = PhysicalUnit(N.getMul()/u.m**2,"Pa")
charge = C = coulomb = coulombs = PhysicalUnit(u.A*u.s,"C")
voltage = V = volt = volts = PhysicalUnit(W.getMul()/u.A,"V")
resistance = ohm = ohms = PhysicalUnit(V.getMul()/u.A,"ohm")
conductance = S = siemens = mho = mhos = PhysicalUnit(u.A/V.getMul(),"S")
capacitance = F = farad = farads = PhysicalUnit(C.getMul()/V.getMul(),"F")
magnetic_flux = Wb = wb = weber = webers = PhysicalUnit(J.getMul()/u.A,"Wb")
magnetic_flux_density = T = tesla = teslas = PhysicalUnit(V.getMul()*u.s/u.m**2,"T")
inductance = H = henry = henrys = PhysicalUnit(V.getMul()*u.s/u.A,"H")
velocity = speed = m_per_s =  PhysicalUnit(u.m/u.s,"m/s")
acceleration = m_per_s2 = PhysicalUnit(u.m/u.s**2,"m/s2")
momentum = kgm_per_s = PhysicalUnit(u.kg*u.m/u.s,"kgm/s")
angular_momentum = kgm2_per_s = PhysicalUnit(u.kg*u.m**2/u.s,"kgm2/s")
inertia = kgm2 = PhysicalUnit(u.kg*u.m**2,"kgm2")
dense = density =kg_per_m3= PhysicalUnit(u.kg/u.m**3,"kg/m3")
wave_number = optical_power = dioptre = D =  PhysicalUnit(1/u.m,"m-1")
lm = lumen = PhysicalUnit(sr.getMul()*u.cd,"lm")
illuminance = lux = lx = PhysicalUnit(sr.getMul()*u.cd/u.m**2,"lx")
rad_per_s = PhysicalUnit(rad.getMul()/u.s,"rad/s")

# Physical constants
c_constant = Rational(299792458)
c_unit = speed_of_light_unit = PhysicalUnit(c_constant*u.m/u.s,"speed of light:m/s")
G0_constant = Rational('6.67428')*u.ten**-11
G0_unit = gravitational_constant_unit = PhysicalUnit(G0_constant*u.m**3/u.kg /u.s**2,"gravitational constant:m3/kgs2")
u0_unit = magnetic_constant_unit = PhysicalUnit(4*pi*u.ten**-7*N.getMul()/A.getMul()**2,"magnetic constant:N/A2")
e0_unit = electric_constant_unit = PhysicalUnit(1/(u0_unit.getMul()*c_unit.getMul()**2),"electric constant:1/u0c2")
Z0_unit = vacuum_impedance_unit = PhysicalUnit(u0_unit.getMul()*c_unit.getMul(),"vacuum impedance:W/A2")

planck_constant = Rational('6.62606896')*u.ten**-34
planck_unit = PhysicalUnit(planck_constant*J.getMul()*u.s,"planck:Js")
hbar_unit = PhysicalUnit(planck_unit.getMul()/(2*pi),"hbar:Js")

avogadro_constant = Rational('6.02214179')*u.ten**23
avogadro_unit = PhysicalUnit(avogadro_constant/mol.getMul(),"avogadro:1/mol")
boltzmann_constant = Rational('1.3806505')*u.ten**-23
boltzmann_unit = PhysicalUnit(boltzmann_constant*J.getMul()/u.K,"boltzman:J/K")

gee_constant = Rational('9.80665')
gee_unit = PhysicalUnit(gee_constant*u.m/u.s**2,"g:m/s2")

# Common time units
ms = millisecond = milliseconds = PhysicalUnit(u.milli*u.s,"ms")
us = microsecond = microseconds = PhysicalUnit(u.micro*u.s,"us")
ns = nanosecond = nanoseconds = PhysicalUnit(u.nano*u.s,"ms")
ps = picosecond = picoseconds = PhysicalUnit(u.pico*u.s,"ms")

minute = minutes = PhysicalUnit(60*u.s,"min")
h = hour = hours = PhysicalUnit(60*minute.getMul(),"h")
day = days = PhysicalUnit(24*hour.getMul(),"day")

anomalistic_year = anomalistic_years = PhysicalUnit(Rational('365.259636')*day.getMul(),"anomalistic_year")
sidereal_year = sidereal_years = PhysicalUnit(Rational('31558149.540')*u.s,"sidereal_year")
tropical_year = tropical_years = PhysicalUnit(Rational('365.24219')*day.getMul(),"tropical_year")
common_year = common_years = PhysicalUnit(Rational('365')*day.getMul(),"common_year")
julian_year = julian_years = PhysicalUnit(Rational('365.25')*day.getMul(),"julian_year")
draconic_year = draconic_years = PhysicalUnit(Rational('346.62')*day.getMul(),"draconic_year")
gaussian_year = gaussian_years = PhysicalUnit(Rational('365.2568983')*day.getMul(),"gaussian_year")
full_moon_cycle = full_moon_cycles = PhysicalUnit(Rational('411.78443029')*day.getMul(),"full_moon_cycle")

year = years = common_year

# rotation

rot_per_s = PhysicalUnit(rot.getMul()/u.s,"rot/s")
rot_per_min = PhysicalUnit(rot.getMul()/minute.getMul(),"rot/min")

# Common length units
km = kilometer = kilometers = PhysicalUnit(u.kilo*u.m,"km")
dm = decimeter = decimeters = PhysicalUnit(u.deci*u.m,"dm")
cm = centimeter = centimeters = PhysicalUnit(u.centi*u.m,"cm")
mm = millimeter = millimeters = PhysicalUnit(u.milli*u.m,"mm")
um = micrometer = micrometers = micron = microns = PhysicalUnit(u.micro*u.m,"um")
nm = nanometer = nanometers = PhysicalUnit(u.nano*u.m,"nm")
pm = picometer = picometers = PhysicalUnit(u.pico*u.m,"pm")

angstrom = PhysicalUnit(u.m*nu.angstrom,"angstrome")

ft = foot = feet = PhysicalUnit(Rational('0.3048')*u.m,"ft")
inch = inches = PhysicalUnit(Rational('25.4')*u.mm,"inch")
yd = yard = yards = PhysicalUnit(3*ft.getMul(),"yd")
mi = mile = miles = PhysicalUnit(5280*ft.getMul(),"mi")

lightyear = PhysicalUnit(c_unit.getMul()*julian_year.getMul(),"lightyear")
astro_unit = PhysicalUnit( 149597870691*u.m,"astro_unit")
pc = parsec = PhysicalUnit(648000/pi*astro_unit.getMul(),"parsec")
kpc = PhysicalUnit(pc.getMul()*u.kilo,"kpc")
Mpc = PhysicalUnit(pc.getMul()*u.mega,"Mpc")
Gpc = PhysicalUnit(pc.getMul()*u.giga,"Gpc")

#common area units
m2 = area = PhysicalUnit(u.m**2,"m2")
cm2 = PhysicalUnit(cm.getMul()**2,"cm2")
mm2 = PhysicalUnit(mm.getMul()**2,"mm2")
km2 = PhysicalUnit(km.getMul()**2,"km2")

inch2= PhysicalUnit(inch.getMul()**2,"inch2")

#common volumn units       
m3 = PhysicalUnit(u.m**3,"m3")
cm3 = PhysicalUnit(cm.getMul()**3,"cm3")
mm3 = PhysicalUnit(mm.getMul()**3,"mm3")

l = liter = liters = PhysicalUnit(m3.getMul()/1000,"l")
dl = PhysicalUnit(l.getMul()*u.deci,"dl")
cl = PhysicalUnit(l.getMul()*u.centi,"cl")
ml = PhysicalUnit(l.getMul()*u.milli,"ml")

inch3 = PhysicalUnit(inch.getMul()**3,"inch3")

m4 = PhysicalUnit(u.m**4,"m4")
cm4 = PhysicalUnit(cm.getMul()**4,"cm4")
mm4 = PhysicalUnit(mm.getMul()**4,"mm4")

#common frequency units
Hz=PhysicalUnit(1/u.s,"Hz")
kHz=PhysicalUnit(Hz.getMul()*u.kilo,"kHz")
MHz=PhysicalUnit(Hz.getMul()*u.mega,"MHz")
GHz=PhysicalUnit(Hz.getMul()*u.giga,"GHz")

#common velocity units
m_per_min = PhysicalUnit(u.m/minute.getMul(),"m/min")
km_per_h = PhysicalUnit(km.getMul()/hour.getMul(),"km/h")

#common mass units
g = gram = grams = PhysicalUnit(u.kg/u.kilo,"g")
mg = milligram = milligrams = PhysicalUnit(g.getMul()*u.milli,"mg")
ug = microgram = micrograms = PhysicalUnit(g.getMul()*u.micro,"ug")

ton = tonne = PhysicalUnit(u.kg*1000,"ton")
amu = amus = atomic_mass_unit = PhysicalUnit(gram.getMul()/avogadro_unit.getMul()/mol.getMul(),"atomic mass unit")
Da = Dalton =  PhysicalUnit(amu.getMul(),"Dalton") #Dalton
kDa = PhysicalUnit(Da.getMul()*u.kilo,"kDa")

lbm = PhysicalUnit(0.45359237*u.kg,"pound mass") # pound mass (international avoirdupois pound)


#common dense unit

kg_per_m3= PhysicalUnit(u.kg/u.m**3,"kg/m3")
ton_per_m3=PhysicalUnit(tonne.getMul()/u.m**3,"t/m3")
g_per_cm3=PhysicalUnit(g.getMul()/cm.getMul()**3,"g/cm3")
kg_per_l=PhysicalUnit(u.kg/l.getMul(),"kg/l")

kg_per_m= PhysicalUnit(u.kg/u.m,"kg/m")
kg_per_m2= PhysicalUnit(u.kg/u.m**2,"kg/m2")

#volume flux unit
m3_per_s = PhysicalUnit(m3.getMul()/u.s,"m3/s")
m3_per_min = PhysicalUnit(m3.getMul()/minute.getMul(),"m3/min")
m3_per_h = PhysicalUnit(m3.getMul()/hour.getMul(),"m3/h")
l_per_s = PhysicalUnit(l.getMul()/u.s,"l/s")
l_per_min = PhysicalUnit(l.getMul()/minute.getMul(),"l/min")
l_per_h = PhysicalUnit(l.getMul()/hour.getMul(),"l/h")

#mass flux unit
kg_per_s = PhysicalUnit(kg.getMul()/u.s,"kg/s")
kg_per_min = PhysicalUnit(kg.getMul()/minute.getMul(),"kg/min")
kg_per_h = PhysicalUnit(kg.getMul()/hour.getMul(),"kg/h")

#mol units

M = molar = PhysicalUnit(mol.getMul()/l.getMul(),"molar")
mmu = mmus = PhysicalUnit(gram.getMul()/mol.getMul(),"mol mass unit")

#common Force units
kN = PhysicalUnit(N.getMul()*u.kilo,"kN")
kgf = PhysicalUnit(u.kg*gee_unit.getMul(),"dyn")
dyn = PhysicalUnit(N.getMul()*u.ten**-5,"dyn")
lbf = PhysicalUnit(lbm.getMul()*gee_unit.getMul(),"lbf")

#common Energy units
mJ = PhysicalUnit(J.getMul()*u.milli,"mJ")
uJ = PhysicalUnit(J.getMul()*u.micro,"uJ")
nJ = PhysicalUnit(J.getMul()*u.nano,"nJ")
pJ = PhysicalUnit(J.getMul()*u.pico,"pJ")
fJ = PhysicalUnit(J.getMul()*u.femto,"fJ")
kJ = PhysicalUnit(J.getMul()*u.kilo,"kJ")
MJ = PhysicalUnit(J.getMul()*u.mega,"MJ")
GJ = PhysicalUnit(J.getMul()*u.giga,"GJ")

erg = PhysicalUnit(J.getMul()*u.ten**-7,"erg")

eV = PhysicalUnit(J.getMul()*1.602176487*u.ten**-19,"eV")
meV = PhysicalUnit(eV.getMul()*u.milli,"meV")
keV = PhysicalUnit(eV.getMul()*u.kilo,"keV")
MeV = PhysicalUnit(eV.getMul()*u.mega,"MeV")
GeV = PhysicalUnit(eV.getMul()*u.giga,"GeV")
TeV = PhysicalUnit(eV.getMul()*u.tera,"TeV")

kgfm = PhysicalUnit(u.kg*u.m*gee_unit.getMul(),"kgfm")

btu = british_thermal_unit = PhysicalUnit(J.getMul()*1055.056,"btu") #British thermal unit
cal = smallcal = PhysicalUnit(J.getMul()*4.184,"cal") #small calorie ("gram calorie")
kcal = PhysicalUnit(4184*J.getMul(),"kcal") #kilocalorie ("large Calorie", "dietary Calorie")

Wh = PhysicalUnit(3600*J.getMul(),"Wh") #watt-hour
kWh = PhysicalUnit(Wh.getMul()*u.kilo,"kWh") # kilowatt-hour

#common Work units
kW = PhysicalUnit(W.getMul()*u.kilo,"kW")
MW = PhysicalUnit(W.getMul()*u.mega,"MW")
GW = PhysicalUnit(W.getMul()*u.giga,"GW")
TW = PhysicalUnit(W.getMul()*u.tera,"TW")

VA = PhysicalUnit(W.getMul(),"VA")
kVA = PhysicalUnit(W.getMul()*u.kilo,"kVA")
MVA = PhysicalUnit(W.getMul()*u.mega,"MVA")
GVA = PhysicalUnit(W.getMul()*u.giga,"GVA")
TVA = PhysicalUnit(W.getMul()*u.tera,"TVA")

HP = PhysicalUnit(kW.getMul()*0.746,"HP")
cal_per_s = PhysicalUnit(cal.getMul()/u.s,"cal/s")
cal_per_h = PhysicalUnit(cal.getMul()/hour.getMul(),"cal/h")
kcal_per_s = PhysicalUnit(kcal.getMul()/u.s,"kcal/s")
kcal_per_h = PhysicalUnit(kcal.getMul()/hour.getMul(),"kcal/h")
kgfm_per_s = PhysicalUnit(kgfm.getMul()/u.s,"kgfm/s")

#common Pressure unit
hPa = PhysicalUnit(Pa.getMul()*u.ten**2,"hPa")
kPa = PhysicalUnit(Pa.getMul()*u.kilo,"kPa")
MPa = PhysicalUnit(Pa.getMul()*u.mega,"MPa")
GPa = PhysicalUnit(Pa.getMul()*u.giga,"GPa")

bar = PhysicalUnit(Pa.getMul()*u.ten**5,"bar")
mbar = PhysicalUnit(bar.getMul()*u.milli,"mbar")
kbar = PhysicalUnit(bar.getMul()*u.kilo,"kbar")
Mbar = PhysicalUnit(bar.getMul()*u.mega,"mbar")

atm = PhysicalUnit(101325*Pa.getMul(),"atm")
torr = PhysicalUnit(atm.getMul()/760,"torr")
mmHg = PhysicalUnit(torr.getMul(),"mmHg")
mtorr = PhysicalUnit(torr.getMul()*u.milli,"mtorr")
psi = PhysicalUnit(lbf.getMul()/inch.getMul()**2,"psi")

mmAq = mmH2O = PhysicalUnit(gee_unit.getMul()*1000*g.getMul()/cm.getMul()**3*mm.getMul(),"mmAq")
mAq = mH2O = PhysicalUnit(mmAq.getMul()/u.milli,"mAq")

kgf_per_cm2 = PhysicalUnit(kgf.getMul()/cm.getMul()**2,"kgf/cm2")
N_per_cm2 = PhysicalUnit(N.getMul()/cm.getMul()**2,"N/cm2")
kN_per_cm2 = PhysicalUnit(N.getMul()/cm.getMul()**2,"N/cm2")
N_per_mm2 = PhysicalUnit(N.getMul()/cm.getMul()**2,"N/cm2")
kN_per_mm2 = PhysicalUnit(N.getMul()/cm.getMul()**2,"N/cm2")

#common degree unit
degF = PhysicalUnit(1.8*u.K,"degF")
degF.setAbsolutePoint(-459.67)
degC = PhysicalUnit(u.K,"degC")
degC.setAbsolutePoint(-273.15)

#common termal unit
per_K = thermal_expantion = PhysicalUnit(1/u.K,"1/K")

J_per_kgK = heat_capacity = PhysicalUnit(J.getMul()/u.kg/u.K,"J/kgK")
kJ_per_kgK = PhysicalUnit(J_per_kgK.getMul()*u.kilo,"kJ/kgK")
cal_per_gK = PhysicalUnit(cal.getMul()/g.getMul()/u.K,"cal/gK")
kcal_per_gK = PhysicalUnit(kcal.getMul()/g.getMul()/u.K,"kcal/gK")

J_per_kg = heat_capacity = PhysicalUnit(J.getMul()/u.kg,"J/kg")
kJ_per_kg = PhysicalUnit(J_per_kg.getMul()*u.kilo,"kJ/kg")
cal_per_g = PhysicalUnit(cal.getMul()/g.getMul(),"cal/g")
kcal_per_g = PhysicalUnit(kcal.getMul()/g.getMul(),"kcal/g")

W_per_mK = thermal_conductivity = PhysicalUnit(W.getMul()/u.m*u.K,"W/mK")
kW_per_mK  = PhysicalUnit(kW.getMul()/u.m*u.K,"kW/mK")
cal_per_hmK = PhysicalUnit(cal_per_h.getMul()/u.m*u.K,"cal/mK")
kcal_per_hmK = PhysicalUnit(kcal_per_h.getMul()/u.m*u.K,"kcal/mK")

W_per_m2K = heat_transefer_coefficient = PhysicalUnit(W.getMul()/u.m**2*u.K,"W/m2K")
kW_per_m2K  = PhysicalUnit(kW.getMul()/u.m**2*u.K,"kW/m2K")
cal_per_hm2K = PhysicalUnit(cal_per_h.getMul()/u.m**2*u.K,"cal/m2K")
kcal_per_hm2K = PhysicalUnit(kcal_per_h.getMul()/u.m**2*u.K,"kcal/m2K")

#common electric unit
mC = PhysicalUnit(C.getMul()*u.milli,"mC")
uC = PhysicalUnit(C.getMul()*u.micro,"uC")
nC = PhysicalUnit(C.getMul()*u.nano,"nC")
Ah = PhysicalUnit(A.getMul()*3600*u.s,"Ah")
mAh = PhysicalUnit(Ah.getMul()*u.milli,"mAh")

mA = PhysicalUnit(A.getMul()*u.milli,"mA")
uA = PhysicalUnit(A.getMul()*u.micro,"uA")
nA = PhysicalUnit(A.getMul()*u.nano,"nA")
pA = PhysicalUnit(A.getMul()*u.pico,"pA")
fA = PhysicalUnit(A.getMul()*u.femto,"fA")

mV = PhysicalUnit(V.getMul()*u.milli,"mV")
uV = PhysicalUnit(V.getMul()*u.micro,"uV")
nV = PhysicalUnit(V.getMul()*u.nano,"nV")
kV = PhysicalUnit(V.getMul()*u.kilo,"kV")
MV = PhysicalUnit(V.getMul()*u.mega,"MV")
GV = PhysicalUnit(V.getMul()*u.giga,"GV")
TV = PhysicalUnit(V.getMul()*u.tera,"TV")

mohm = PhysicalUnit(ohm.getMul()*u.milli,"mohm")
kohm = PhysicalUnit(ohm.getMul()*u.kilo,"kohm")
Mohm = PhysicalUnit(ohm.getMul()*u.mega,"Mohm")
Gohm = PhysicalUnit(ohm.getMul()*u.giga,"Gohm")

mS = PhysicalUnit(S.getMul()*u.milli,"mS")
uS = PhysicalUnit(S.getMul()*u.micro,"uS")
nS = PhysicalUnit(S.getMul()*u.nano,"nS")

mT = PhysicalUnit(T.getMul()*u.milli,"mT")
uT = PhysicalUnit(T.getMul()*u.micro,"uT")
nT = PhysicalUnit(T.getMul()*u.nano,"nT")
G = PhysicalUnit(T.getMul()*u.ten**-4,"G")
mG = PhysicalUnit(G.getMul()*u.milli,"mG")
uG = PhysicalUnit(G.getMul()*u.micro,"uG")
kG = PhysicalUnit(G.getMul()*u.kilo,"kG")
Oe = PhysicalUnit(1000/4/pi*A.getMul()/u.m,"Oe")

mF = PhysicalUnit(F.getMul()*u.milli,"mF")
uF = PhysicalUnit(F.getMul()*u.micro,"uF")
nF = PhysicalUnit(F.getMul()*u.nano,"nF")
pF = PhysicalUnit(F.getMul()*u.pico,"pF")
fF = PhysicalUnit(F.getMul()*u.femto,"fF")
aF = PhysicalUnit(F.getMul()*u.atto,"aF")

mH = PhysicalUnit(H.getMul()*u.milli,"mH")
uH = PhysicalUnit(H.getMul()*u.micro,"uH")
nH = PhysicalUnit(H.getMul()*u.nano,"nH")

#some constant
sigmaSB_constant = Rational(5.670373*u.ten**-8)
sigmaSB_unit = PhysicalUnit(sigmaSB_constant*W.getMul()/(u.m**2*u.K**4),"Stefan-Boltzmann constant:W/m2K4")
alphaFS_constant = Rational(7.2973525698*u.ten**-3)  
alphaFS_unit = PhysicalUnit(alphaFS_constant*U.getMul(),"fine-structure constant:U")

#Constants--chemistry, atomic physics, electrons
Rgas_unit = PhysicalUnit(boltzmann_unit.getMul()*avogadro_unit.getMul(),"ideal gas constant:J/Kmol")
elementary_charge_unit = PhysicalUnit(1.602176565*u.ten**-19*C.getMul(),"charge of proton:C")
uBor_unit = PhysicalUnit(9.27400968*u.ten**-24*J.getMul()/T.getMul(),"Bohr magnetron:J/T")
uNuc_unit = PhysicalUnit(5.05078353*u.ten**-27*J.getMul()/T.getMul(),"nuclear magneton:J/T")
aBohr_unit = PhysicalUnit(0.52917721092*u.ten**-10*u.m,"Bohr radius:m")

me_unit = PhysicalUnit(9.10938291*u.ten**-31*u.kg,"electron mass:kg")
ma_unit = PhysicalUnit(1.672621777*u.ten**-27*u.kg,"proton mass:kg")
mn_unit = PhysicalUnit(1.674927351*u.ten**-27*u.kg,"neutron mass:kg")
Rinf_unit = PhysicalUnit(10973731.568539/u.m,"Rydberg constant:1/m")
Ry_unit = PhysicalUnit(me_unit.getMul()*(elementary_charge_unit.getMul())**4/8/e0_unit.getMul()**2/planck_unit.getMul()**2,"Rydberg energy constant:J") 
Eh_unit = PhysicalUnit(2*Ry_unit.getMul(),"hartree energy constant:J")
ARichardson_unit = PhysicalUnit(4*pi*elementary_charge_unit.getMul()*me_unit.getMul()*boltzmann_unit.getMul()**2/planck_unit.getMul()**3,"Richardson constant:A/m2K2")
Phi0_unit = PhysicalUnit(planck_unit.getMul()/2/elementary_charge_unit.getMul(),"magnetic flux quantum:Wb")
KJos_unit = PhysicalUnit(Phi0_unit.getMul()**-1,"Josephson constant:Hz/V")
RKlitz_unit = PhysicalUnit(planck_unit.getMul()/elementary_charge_unit.getMul()**2,"von Klitzing constant:ohm")

#Constants--astronomical and properties of earth
REarth_unit = PhysicalUnit(6371*km.getMul(),"radius of earth")
MEarth_unit = PhysicalUnit(5.9722*u.ten**24*u.kg,"mass of earth")
Msolar_unit = PhysicalUnit(1.98892*u.ten**30*u.kg,"mass of earth")

def find_unit(unit=None,getinstance=False):
    """
    Return a list of matching units names.
    if quantity is a string -- units containing the string `quantity`
    if quantity is a unit -- units having matching base units

    Examples
    ========

    >>> from sympy.physics import units as u
    >>> u.find_unit('charge')
    ['charge']
    >>> u.find_unit(u.charge)
    ['C', 'charge', 'coulomb', 'coulombs']
    >>> u.find_unit('volt')
    ['volt', 'volts', 'voltage']
    >>> u.find_unit(u.inch**3)[:5]
    ['l', 'cl', 'dl', 'ml', 'liter']
    """
    import PhysicalUnit as PU
    rv = []
    if unit == None:
        for i in dir(PU):
            tu=eval("PU."+i)
            if isinstance(tu,PU.PhysicalUnit):
                if getinstance:
                    rv.append(tu)
                else:
                    rv.append(i)
    elif isinstance(unit, str):
        if getinstance:
            for i in dir(PU):
                if unit in i:
                    rv.append(eval("PU."+i))
        else:
            rv = [i for i in dir(u) if unit in i]
    else:
        for i in dir(PU):
            try:
                tu=eval("PU."+i)
                if _eqDimention(unit,tu):
                    if getinstance:
                        rv.append(tu)
                    else:
                        rv.append(str(i))
            except Exception:
                pass
    return rv

def test_isolative_unit_list():
    
    import PhysicalUnit as PU
    rv = fancy_isolative_unit_list(True)
    for i in dir(PU):
        tu=eval("PU."+i)
        if not isinstance(tu,PU.PhysicalUnit):
            continue
        else:
            flag=True
            for j in rv:
                if _eqDimention(tu,j):
                    flag=False
            if flag==True:
                rv.append(tu)
    
    return [i.getName() for i in rv]

def base_unit_list():
    t=[U,m,kg,s,A,K,cd,mol]
    return [i.getName() for i in t]

def fancy_isolative_unit_list(flag=False):
    t=[U,m,kg,s,A,K,cd,mol,\
       m2,m3,m4,wave_number,\
       Hz,m_per_s2,N,Pa,J,W,kgm_per_s,\
       m_per_s,m3_per_s,kg_per_s,\
       kg_per_m,kg_per_m2,kg_per_m3,\
       rad_per_s,kgm2,kgm2_per_s,\
       per_K,J_per_kg,J_per_kgK,W_per_mK,W_per_m2K,\
       V,C,ohm,S,F,Wb,T,H,lm,lx]
    if flag:
        return t
    else:
        return [i.getName() for i in t]
if __name__=="__main__":     
    print find_unit(planck_unit)
    print base_unit_list()
    print funcy_isolative_unit_list()
