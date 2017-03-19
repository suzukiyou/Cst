#coding:utf-8

import codecs
import sys,os
import isbnlib
import requests
import xml.etree.ElementTree as ET
from StringIO import StringIO

dec=codecs.getdecoder("utf-8")
enc=codecs.getencoder("utf-8")
enc932=codecs.getencoder("cp932")

class BibIdentifier(object):
    ISBN10="ISBN10"
    ISBN13="ISBN13"
    ISBN="ISBN"
    JPNO="JPNO"
    identifier_types=(ISBN10,ISBN13,ISBN,JPNO)
    ISBNtypes=(ISBN10,ISBN13,ISBN)
    
    def __init__(self,identifier_type,id):
        if identifier_type in self.identifier_types:
            self.identifier_type=identifier_type
            self.identifier={self.identifier_type:id}
            for i in self.identifier_types:
                if not i == self.identifier_type:
                    self.identifier[i]=None
        else:
            raise KeyError("identifier not exist in list : (ISBN10,ISBN13,ISBN,JPNO,MSMARCNO)")
    @classmethod
    def create_by_ISBN(cls,ISBN):
        t=cls.ISBN_type_check(ISBN)
        return cls(t,ISBN)
    @classmethod
    def create_by_JPNO(cls,JPNO):
        return cls(cls.JPNO,JPNO)
    @classmethod
    def ISBN_type_check(cls,isbn): 
        if isbnlib.is_isbn10(isbn):
            return cls.ISBN10
        elif isbnlib.is_isbn13(isbn):
            return cls.ISBN10
        else:
            return cls.ISBN
    def get_identifier(self):
        return self.identifier[self.identifier_type]

    
class DCNDL(ET.Element):
    rdfns={"rdfs":"http://www.w3.org/2000/01/rdf-schema#",
           "dcndl":"http://ndl.go.jp/dcndl/terms/",
           "owl":"http://www.w3.org/2002/07/owl#",
           "dc":"http://purl.org/dc/elements/1.1/",
           "foaf":"http://xmlns.com/foaf/0.1/",
           "rdf":"http://www.w3.org/1999/02/22-rdf-syntax-ns#",
           "dcterms":"http://purl.org/dc/terms/"}
    def __init__(self,tree):
        self.tree=tree
        self.element=self.tree.getroot()
        self.JPNO=self.tree.find(".//dcterms:identifier[@rdf:datatype='http://ndl.go.jp/dcndl/terms/JPNO']",self.rdfns)
        self.ISBN=self.tree.find(".//dcterms:identifier[@rdf:datatype='http://ndl.go.jp/dcndl/terms/ISBN']",self.rdfns)
        self.title=self.tree.find(".//dcterms:title",self.rdfns)
        self.creator=self.tree.find(".//dcterms:creator/foaf:Agent/foaf:name",self.rdfns)
        self.publisher=self.tree.find(".//dcterms:publisher/foaf:Agent/foaf:name",self.rdfns)
        self.page=self.tree.find(".//dcterms:extent",self.rdfns)
        self.subtitle=self.tree.find(".//dcndl:partInformation/rdf:Description/dcterms:title",self.rdfns)

    def __getattr__(self,key):
        if not hasattr(self.__class__,key):
            return getattr(self.element,key)
        else:
            raise AttributeError
        
    def __repr__(self,sep=",",keys=["JPNO","ISBN","title","subtitle","creator","publisher","page"]):
        l=[]
        for i in [getattr(self,j) for j in keys]:
            if isinstance(i,ET.Element):
                l.append(i.text)
            else:
                l.append("")
        return sep.join(l)
    def local_path(self):
        return os.path.join(os.path.abspath(os.path.dirname(__file__)),"saved","jpno"+self.JPNO.text+".xml")            
    @classmethod
    def parse(cls,rdf):
        f=StringIO(enc(rdf)[0])
        return cls(ET.parse(f))

