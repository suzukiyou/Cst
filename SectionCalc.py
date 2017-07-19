#coding:utf-8

import numpy as np
from tabulate import tabulate

class tabulateConfigure(object):
    @classmethod
    def repr(cls,obj):
        return u"<%s:%s>"%(unicode(obj.__class__.__name__),unicode(getattr(obj,"name")))
    @classmethod
    def description(cls,obj,index,operate=None):
        t=[[i,getattr(obj,i)] for i in index if getattr(obj,i)!=None]
        return tabulate(t)

class BuildHSection(object):
    description_attribute=("name,H,B_upper,tw,tf_upper,B_lower,tf_lower,"
                           +"A,Aw,Af_upper,Af_lower,"
                           +"Iy,Iz,Zy_upper,Zy_lower,Zz_upper,Zz_lower,iy,iz,"
                           +"Zyp,Zzp,J,Iw,xn_elas,xn_plas").split(",")
    tabulateConfigure=tabulateConfigure
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
    def __repr__(self):
        return self.tabulateConfigure.repr(self)
    def description(self,table=None,operate=None):
        if table==None:
            return self.tabulateConfigure.description(self,self.description_attribute,operate)
        else:
            return self.tabulateConfigure.description(self,table,operate)

if __name__=="__main__":
    s400=BuildHSection("H400x400x13x21",400,400,13,21)
    print s400.description()
    s200=BuildHSection("H200x100x5.5x8",200,100,5.5,8)
    print s200.description()
