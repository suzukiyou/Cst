#coding:utf-8

import numpy as np
from tabulate import tabulate

class TabulateConfigure(object):
    @classmethod
    def repr(cls,obj):
        return u"<%s:%s>"%(unicode(obj.__class__.__name__),unicode(getattr(obj,"name")))
    @classmethod
    def description(cls,obj,index,operate=None):
        t=[[i,getattr(obj,i)] for i in index if getattr(obj,i)!=None]
        return tabulate(t)

class TabulateConfigureMixin(object):
    def __repr__(self):
        return self.tabulateConfigure.repr(self)
    def description(self,table=None,operate=None):
        if table==None:
            return self.tabulateConfigure.description(self,self.description_attribute,operate)
        else:
            return self.tabulateConfigure.description(self,table,operate)
    
class BuildHSection(TabulateConfigureMixin):
    description_attribute=("name,H,B_upper,tw,tf_upper,B_lower,tf_lower,"
                           +"A,Aw,Af_upper,Af_lower,"
                           +"Iy,Iz,Zy_upper,Zy_lower,Zz_upper,Zz_lower,iy,iz,"
                           +"Zyp,Zzp,J,Iw,xn_elas,xn_plas").split(",")
    tabulateConfigure=TabulateConfigure
    def __init__(self,name,H,B,tw,tf,B2=None,tf2=None):
        self.name=name
        self.H=np.float32(H)
        self.B_upper=np.float32(B)
        self.tw=np.float32(tw)
        self.tf_upper=np.float32(tf)
        if B2!=None:
            self.B_lower=np.float32(B2)
            self.tf_lower=np.float32(tf2)
        else:
            self.B_lower=np.float32(B)
            self.tf_lower=np.float32(tf)
        self.getArea()
        self.getElasticCentroid()
        self.getInertia()
        self.getPlasticCentroid()
        self.getPlasticInertia()
        self.getStVenant()
        self.getInertiaOfLateralBuckling()
    def getArea(self):
        self.Af_upper=self.B_upper*self.tf_upper
        self.Aw=(self.H-self.tf_upper-self.tf_lower)*self.tw
        self.Af_lower=self.B_lower*self.tf_lower
        self.A=self.Af_upper+self.Aw+self.Af_lower
    def getElasticCentroid(self):
        A=self.A
        uf=self.B_upper*(-self.tf_upper**2/2.+0)
        w=self.tw*(-(self.H-self.tf_lower)**2/2.+self.tf_upper**2/2.)
        lf=self.B_lower*(-self.H**2/2.+(self.H-self.tf_lower)**2/2.0)
        self.xn_elas=-(uf+w+lf)/A
        return self.xn_elas
    def getPlasticCentroid(self):
        r=self.Af_lower+(self.H-self.tf_lower)*self.tw
        l=self.Af_upper-self.tf_upper*self.tw
        self.xn_plas=(r-l)/2.0/self.tw
        return self.xn_plas
    def getIy(self):
        uf0=self.B_upper*np.power(self.tf_upper,3)/12.0
        uf_trans=np.power(self.tf_upper/2.0-self.xn_elas,2)\
                  *self.B_upper*self.tf_upper
        w0=self.tw*np.power(self.H-self.tf_upper-self.tf_lower,3)/12.0
        w_trans=np.power(self.H/2.-self.xn_elas,2)*self.tw\
                 *(self.H-self.tf_upper-self.tf_lower)
        lf0=self.B_lower*np.power(self.tf_lower,3)/12.0
        lf_trans=np.power(self.H-self.tf_lower/2.-self.xn_elas,2)\
                  *self.B_lower*self.tf_lower
        self.Iy=uf0+uf_trans+w0+w_trans+lf0+lf_trans
        return self.Iy
    def getIz(self):
        self.Iz=np.power(self.tw,3)*(self.H-self.tf_upper-self.tf_lower)/12.0\
                 +np.power(self.B_upper,3)*self.tf_upper/12.0\
                 +np.power(self.B_lower,3)*self.tf_lower/12.0
        return self.Iz
    def getInertia(self):
        self.getIy()
        self.getIz()
        self.Zy_upper=self.Iy/self.xn_elas
        self.Zy_lower=self.Iy/(self.H-self.xn_elas)
        self.Zz_upper=self.Iz/self.B_upper*2
        self.Zz_lower=self.Iz/self.B_lower*2
        self.iy=np.sqrt(self.Iy/self.A)
        self.iz=np.sqrt(self.Iz/self.A)
    def getZyp(self):
        x=np.array([0-self.xn_plas,\
                    self.tf_upper-self.xn_plas,\
                    self.H-self.tf_lower-self.xn_plas,\
                    self.H-self.xn_plas])
        x2=np.power(x,2)
        self.Zyp=(self.B_upper*(x2[0]-x2[1])\
                  +self.tw*(x2[1]+x2[2])\
                  +self.B_lower*(x2[3]-x2[2]))/2.0
        return self.Zyp
    def getZzp(self):
        self.Zzp=self.tf_upper*np.power(self.B_upper,2)/4.0\
                  +(self.H-self.tf_upper-self.tf_lower)*np.power(self.tw,2)/4.0\
                  +self.tf_lower*np.power(self.B_lower,2)/4.0
        return self.Zzp
    def getPlasticInertia(self):
        self.getZyp()
        self.getZzp()
    def getStVenant(self):
        Juf=self.B_upper*np.power(self.tf_upper,3)
        Jw=(self.H-self.tf_upper-self.tf_lower)*np.power(self.tw,3)
        Jlf=self.B_lower*np.power(self.tf_lower,3)
        self.J=(Juf+Jw+Jlf)/3.0
        return self.J
    def getInertiaOfLateralBuckling(self):
        Iu=self.tf_upper*np.power(self.B_upper,3)/12.0
        Il=self.tf_lower*np.power(self.B_lower,3)/12.0
        eu=self.xn_elas
        el=self.H-self.xn_elas
        e=(el*Il-eu*Iu)/(Il+Iu)
        self.Iw=np.power(self.H-self.tf_upper/2.0-self.tf_lower/2.0,2)\
                 *Iu*Il/(Iu+Il)
        return self.Iw
class Material(TabulateConfigureMixin):
    description_attribute=("name,F,E,G,T").split(",")
    tabulateConfigure=TabulateConfigure
    def __init__(self,name,F,E,G=None,T=None):
        self.name=name
        self.F=np.float32(F)
        self.E=np.float32(E)
        self.G=np.float32(G)
        self.T=np.float32(T)
        self.culLambda()
    def culLambda(self):
        self.Lambda=np.sqrt(np.pi*np.pi*self.E/0.6/self.F)
        return self.Lambda
class LimitState(TabulateConfigureMixin):
    description_attribute=["name"]
    tabulateConfigure=TabulateConfigure
    def __init__(self,name):
        self.name=name

class LimitState_S_AIJ2002(LimitState):
    name="AIJ2002"
    @classmethod
    def culft(cls,obj):
        elem=obj.element
        obj.ft_long=elem.material.F/1.5
        obj.ft_short=obj.ft_long*1.5
    @classmethod
    def culfb(cls,obj):
        elem=obj.element
        obj.C=1
        obj.Af=elem.section.Af_upper+elem.section.Aw/6.0
        obj.I_f=elem.section.Iz/2.0-np.power(elem.section.tw,3)\
                 *(elem.section.H-elem.section.tf_upper-elem.section.tf_lower)\
                 /12.0/3.0
        obj.i_f=np.sqrt(obj.I_f/obj.Af)
        obj.fb_long_1=(1-0.4*np.power(elem.lb/obj.i_f/elem.material.Lambda,2)\
                       /obj.C)*obj.ft_long
        obj.fb_long_2=0.433*elem.material.E/(elem.lb*elem.section.H/obj.Af/2.0)
        obj.fb_long=np.min([obj.ft_long,np.max([obj.fb_long_1,obj.fb_long_2])])
        obj.fb_short=obj.fb_long*1.5
"""
class LimitState_S_JSME2005(LimitState_S_AIJ2002):
    name="JSME2005"
""" 

class LimitState_S_AIJ2005(LimitState):
    name="AIJ2005"
    @classmethod
    def culft(cls,obj):
        elem=obj.element
        obj.ft_long=elem.material.F/1.5
        obj.ft_short=obj.ft_long*1.5
    @classmethod
    def culfb(cls,obj):
        elem=obj.element
        obj.culft()
        obj.C=1
        F=elem.material.F
        E=elem.material.E
        Iz=elem.section.Iz
        Iw=elem.section.Iw
        G=elem.material.G
        J=elem.section.J
        lb=elem.lb
        obj.Me=obj.C*np.sqrt(np.power(np.pi,4)*E*E*Iz*Iw/np.power(lb,4)\
                             +np.power(np.pi,2)*E*Iz*G*J/np.power(lb,2))
        obj.Mz=F*elem.section.Zz_upper
        obj.plambdab=0.3
        obj.elambdab=1.0/np.sqrt(0.6)
        obj.lambdab=np.sqrt(obj.Mz/obj.Me)
        obj.nu_fb=3./2.+2./3.*np.power(obj.lambdab/obj.elambdab,2)
        obj.fb_long_1=F/obj.nu_fb
        obj.fb_long_2=(1-0.4*(obj.lambdab-obj.plambdab)\
                       /(obj.elambdab-obj.plambdab))*F/obj.nu_fb
        obj.fb_long_3=1/np.power(obj.lambdab,2)*F/2.17
        if obj.lambdab<obj.plambdab:
            obj.fb_long=obj.fb_long_1
        elif obj.lambdab>obj.elambdab:
            obj.fb_long=obj.fb_long_3
        else:
            obj.fb_long=obj.fb_long_2
        obj.fb_short=obj.fb_long*1.5

class LimitState_S_EN1993(LimitState):
    name="EN1993"
    @classmethod
    def culCrossSectionClassification(cls,obj):
        elem=obj.element
        sec=elem.section
        cw=sec.H-sec.tf_upper-sec.tf_lower-2*sec.tw
        cf=(np.min([sec.B_upper,sec.B_lower])-sec.tw*3.)/2.0
        alpha=1.0/cw*(cw/2.0+1.0/2.0*elem.N/sec.tw/elem.material.F)
        obj.alpha=alpha
        obj.width_to_thickness_class=1
    @classmethod
    def culNplRd(cls,obj):
        elem=obj.element
        sec=elem.section
        obj.NplRd=sec.A*elem.material.F/1.0
    @classmethod
    def culMplRd(cls,obj):
        elem=obj.element
        sec=elem.section
        obj.MplyRd=sec.Zyp*elem.material.F/1.0
        obj.MplzRd=sec.Zzp*elem.material.F/1.0
        n=elem.M/obj.NplRd
        a=sec.Aw/sec.A
        obj.MNyRd=obj.MplyRd*(1-n)/(1-0.5*a)
        obj.MNzRd=obj.MplzRd*(1-np.power((n-a)/(1-a),2))
    @classmethod
    def culMbRd(cls,obj):
        elem=obj.element
        sec=elem.section
        G=elem.material.G
        E=elem.material.E
        Iz=sec.Iz
        It=sec.J
        Iw=sec.Iw
        Lb=elem.lb
        Mcr=np.pi/Lb*np.sqrt(np.pi*np.pi*E*Iw/Lb/Lb+G*It*E*Iz)
        lambdaLT=np.sqrt(sec.Zyp*elem.material.F/Mcr)
        alphaLT=0.21
        phiLT=0.5*(1+alphaLT*(lambdaLT-0.2)+np.power(lambdaLT,2))
        xiLT=np.min([1.0/(phiLT+np.sqrt(np.power(phiLT,2)\
                                        -np.power(lambdaLT,2)))\
                     ,1.0])
        MbRd=xiLT*sec.Zyp*elem.material.F/1.0

        obj.Mcr=Mcr
        obj.lambdaLT=lambdaLT
        obj.phiLT=phiLT
        obj.xiLT=xiLT
        obj.MbRd=MbRd

class LimitStateForElement(TabulateConfigureMixin):
    description_attribute=["name","element","LimitState"]
    tabulateConfigure=TabulateConfigure

    def __init__(self,element,ls):
        self.element=element
        self.LimitState=ls
        self.name=self.LimitState.name+"->"+self.element.name
    def cul(self):
        raise NotImplementedError
class LimitStateForElementOnAIJ(LimitStateForElement):
    description_attribute=["name","element","LimitState",\
                           "ft_long","fb_long","ft_short","fb_short"]

    def culft(self):
        self.LimitState.culft(self)
    def culfb(self):
        self.LimitState.culfb(self)
    def cul(self):
        self.culft()
        self.culfb()

class LimitStateForElementOnEN(LimitStateForElement):
    description_attribute=["name","element","LimitState",\
                           "NplRd","MplyRd","MplzRd","MbRd"]
    def culNplRd(self):
        self.LimitState.culNplRd(self)
    def culMplRd(self):
        self.LimitState.culMplRd(self)
    def culMbRd(self):
        self.LimitState.culMbRd(self)
    def cul(self):
        self.culNplRd()
        self.culMplRd()
        self.culMbRd()

class Element(TabulateConfigureMixin):
    description_attribute=["name","material","section","length","lb"]
    tabulateConfigure=TabulateConfigure

    def __init__(self,name,material=None,section=None,length=None,lb=None):
        self.name=name
        self.material=material
        self.section=section
        self.length=np.float32(length)
        self.lb=np.float32(lb)
       
        
if __name__=="__main__":
    s400=BuildHSection("H400x400x13x21",400,400,13,21)
    print s400.description()
    s200=BuildHSection("H200x100x5.5x8",200,100,5.5,8)
    print s200.description()
    s1=BuildHSection("ISP550x300",550,200,12,25,300,25)
    print s1.description()
    
    q1=1.35*30.+1.5*6.0
    M1=q1*8.0**2/8.0*1000.**2

    m1=Material("S355",355,210000,79000)
    e1=Element("B1",material=m1,section=s400,length=8000,lb=8000)
    e1.N=0.
    e1.M=M1

    l2002_e1=LimitStateForElementOnAIJ(e1,LimitState_S_AIJ2002)
    l2002_e1.cul()
    print l2002_e1.description()
    l2005_e1=LimitStateForElementOnAIJ(e1,LimitState_S_AIJ2005)
    l2005_e1.cul()
    print l2005_e1.description()
    EN3_e1=LimitStateForElementOnEN(e1,LimitState_S_EN1993)
    EN3_e1.cul()
    print EN3_e1.description()


