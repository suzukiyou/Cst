#coding:utf-8

class PhysicalMean(object):
    _mean=None
    def __init__(self,mean):
        assert isinstance(mean,basestring)
        self._mean=mean
    def getMean(self):
        return self._mean
    def __repr__(self):
        return "<Physical Mean>:"+self._mean
