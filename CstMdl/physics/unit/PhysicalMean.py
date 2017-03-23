#coding:utf-8
import numpy as np

class PhysicalMean(object):
    _mean=None
    ROUNDUP="ROUNDUP"
    ROUNDDOWN="ROUNDDOWN"
    ROUND="ROUND"
    def __init__(self,mean,symbol=None,prec=None,roundness=None):
        assert isinstance(mean,basestring)
        self._mean=mean
        self._symbol=symbol
        self._prec=prec
        self._round=roundness
    def getMean(self):
        return self._mean
    def setMean(self,mean):
        self._mean=mean
    def getSymbol(self):
        return self._symbol
    def setSymbol(self,symbol):
        self._symbol=symbol
    def __repr__(self):
        return "<Physical Mean>:"+self._mean
    def precDefined(self):
        return isinstance(self._prec,np.int)
    def roundDefined(self):
        return self._round in [self.ROUNDUP,self.ROUNDDOWN,self.ROUND]
    def figureDefined(self):
        return self.precDefined() and self.roundDefined()
    def getPrec(self):
        return self.prec
    def setPrec(self,prec):
        self._prec=prec
    def getRound(self):
        return self._round
    def getRoundFunc(self):
        if self._round==self.ROUNDUP:
            return np.ceil
        if self._round==self.ROUNDDOWN:
            return np.floor
        if self._round==self.ROUND:
            return np.round        
    def setRound(self,roundness):
        self._round=roundness
        
    
