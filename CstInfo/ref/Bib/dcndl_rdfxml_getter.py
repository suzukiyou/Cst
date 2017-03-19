
#coding:utf-8

import urllib2
import codecs
import sys,os
import traceback
import requests
from CstInfo.ref.Bib.BibIdentifier import BibIdentifier,DCNDL

dec=codecs.getdecoder("utf-8")
enc=codecs.getencoder("utf-8")
enc932=codecs.getencoder("cp932")

class DCNDL_getter(object):
    def __init__(self):
        self.identifier=None
        self.request=None
        self.rdf=None
        self.dcndl=None
    def get(self,identifier):
        self.request_by_identifier(identifier)
        source=self.request.text
        self.rdf=source.replace(u"&lt;",u"<").replace(u"&gt;",u">").replace(u"&quot;",u"\"")
        self.dcndl=DCNDL.parse(self.rdf)
        
    def request_by_identifier(self,identifier):

        self.identifier=identifier
        ndlurl=r"http://iss.ndl.go.jp/api/sru"
        apipram={"operation":"searchRetrieve",\
                 "version":"1.2",\
                 "recordSchema":"dcndl",\
                 "onlyBib":"True",\
                 "recordPacking":"string",\
                 "query":""}
        apipram["query"]=self.create_query()
        self.request=requests.get(ndlurl,params=apipram)
        
    def create_query(self):
        if self.identifier.identifier_type in BibIdentifier.ISBNtypes:
            return "".join(["isbn=\"",self.identifier.get_identifier(),\
                            "\" AND dpid=iss-ndl-opac"])
        elif self.identifier.identifier_type in BibIdentifier.JPNO:
            return "".join(["jpno=\"",self.identifier.get_identifier(),\
                            "\" AND dpid=iss-ndl-opac"])

if __name__=="__main__":
  
    message=u"""
dcndl_rdfxml_getter.py
国会図書館のAPIを叩いて文書情報を保存したりします。
get  XXXXXXXXXXXX:    apiを叩いてdcndl-rdfxmlを取得する
display :    取ってきたrdfxml要約(タイトル・著者)を表示する
dispall :    取ってきたrdfxml全文を表示する
save    :    取ってきたrdfxmlファイルを保存する
             ファイル名は./saved/${JPNO}.xmlです
upload  :    (未実装)アップロードする
quit    :    終了する
使い方
>>>get XXXXXXXXXX
->dcndlを取ってきてサマリ表示する
->XXX~にはISBN(10/13)を指定してください。
->JPNOを指定もできます。指定する際はJPNOXXXXXXXXとすること。
>>>save
->表示されたやつを/save/以下に保存する。
"""
    print message

    getter=None
    while True:
        t=raw_input("(dcndl_rdfxml_getter)>>>")
        if t=="quit":
            break
        cmd_param=t.split(" ")
        if cmd_param[0]=="get":
            try:
                tidf=cmd_param[1]
                if len(tidf)==10 or len(tidf)==13:
                    tidf=BibIdentifier.create_by_ISBN(tidf)
                elif tidf[0:4]=="JPNO":
                    tidf=BibIdentifier.create_by_JPNO(tidf[4:])
                else:
                    print "ignore identifier : %s"%tidf
                    continue
                getter=DCNDL_getter()
                getter.get(tidf)
                print enc932(getter.dcndl.__repr__())[0]
            except:
                info=sys.exc_info()
                tbinfo=traceback.format_tb(info[2])
                print 'Python Error.'.ljust( 80, '=' )
                for tbi in tbinfo:
                    print tbi
                print '  %s' % str( info[1] )
                print '\n'.rjust( 80, '=' )
                
        if cmd_param[0]=="display":
            if len(cmd_param)==1:
                print enc932(getter.dcndl.__repr__())[0]
            else:
                print enc932(getter.dcndl.__repr__(keys=cmd_param[1:]))[0]
        if cmd_param[0]=="dispall":
            print enc932(getter.rdf)[0]
        if cmd_param[0]=="save":
            filename=getter.dcndl.local_path()
            try:
                orfp=codecs.open(filename,"w","utf-8")
                orfp.write(getter.rdf)
                orfp.close()
                print "save as : %s"%filename
            except:
                info=sys.exc_info()
                tbinfo=traceback.format_tb(info[2])
                print 'Python Error.'.ljust( 80, '=' )
                for tbi in tbinfo:
                    print tbi
                print '  %s' % str( info[1] )
                print '\n'.rjust( 80, '=' )
        if cmd_param[0]=="help":
            print message
        if cmd_param[0]=="upload":
            pass


