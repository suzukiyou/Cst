import numpy as np
from scipy.optimize import root
from functools import partial
from collections import OrderedDict


class Graph(object):
    def __init__(self,name=""):
        self.name=name
        self.vertexList = OrderedDict()
        self.edgeList = OrderedDict()

    def add_vertex(self,vartex):
        self.vertexList[vartex.name] = vartex
        vartex.set_graph(self)
        return vartex

    def add_vertex2(self, name):
        v = Vertex(name)
        return self.add_vertex(v)

    def add_vertexes(self,vartexes):
        return [self.add_vertex(v) for v in vartexes]

    def add_vertexes2(self, names, delimiter=","):
        return [self.add_vertex2(i) for i in names.split(delimiter)]

    def add_edge(self,edge):
        self.edgeList[edge.toString()] = edge
        edge.set_graph(self)
        return edge

    def add_edge2(self, edgeStr, name=""):
        preStr, postStr = edgeStr.split("->")
        pre = self.vertexList[preStr]
        post = self.vertexList[postStr]
        e = Edge(pre, post)
        return self.add_edge(e)

    def add_edges(self,edges):
        return [self.add_edge(e) for e in edges]

    def add_edges2(self, edgesStr, delimiter=","):
        return [self.add_edge2(i) for i in edgesStr.split(",")]

    def set_index(self):
        self.vertexIndexes=[vStr for vStr in self.vertexList]
        self.edgeIndexes=[eStr for eStr in self.edgeList]

    def calc_linkedList(self):
        self.linkedListDownStream=[]
        self.linkedListUpStream=[]
        for vStr in self.vertexIndexes:
            tempLinkedListUpStream=[]
            tempLinkedListDownStream=[]
            v=self.vertexList[vStr]
            for eStr in self.edgeIndexes:
                e=self.edgeList[eStr]
                if e.pre is v:
                    tempLinkedListDownStream.append(e)
                elif e.post is v:
                    tempLinkedListUpStream.append(e)
            self.linkedListDownStream.append(tempLinkedListDownStream)
            self.linkedListUpStream.append(tempLinkedListUpStream)
            v.linkedListDownStream=tempLinkedListDownStream
            v.linkedListUpStream=tempLinkedListUpStream
        return [self.linkedListDownStream,self.linkedListUpStream]

    def calc_vertexAdjunction(self):
        self.vertexAdjunction=[]
        for vStr in self.vertexIndexes:
            base=np.zeros(len(self.vertexIndexes))
            vIndex=self.vertexIndexes.index(vStr)
            for upEdge in self.linkedListUpStream[vIndex]:
                base[self.vertexIndexes.index(upEdge.pre.name)]=-1
            for downEdge in self.linkedListDownStream[vIndex]:
                base[self.vertexIndexes.index(downEdge.post.name)]=1
            self.vertexAdjunction.append(base)
        self.vertexAdjunction=np.array(self.vertexAdjunction)
        return self.vertexAdjunction

    def calc_edgeAdjunction(self):
        self.edgeAdjunction=[]
        for vStr in self.vertexIndexes:
            base=np.zeros(len(self.edgeIndexes))
            vIndex=self.vertexIndexes.index(vStr)
            for upEdge in self.linkedListUpStream[vIndex]:
                base[self.edgeIndexes.index(upEdge.toString())]=-1
            for downEdge in self.linkedListDownStream[vIndex]:
                base[self.edgeIndexes.index(downEdge.toString())]=1
            self.edgeAdjunction.append(base)
        self.edgeAdjunction=np.array(self.edgeAdjunction)
        return self.edgeAdjunction


    def __getitem__(self, key):
        try:
            return self.vertexList[key]
        except KeyError:
            return self.edgeList[key]


class Vertex(object):
    def __init__(self, name):
        self.name = name
        self.linkedListUpStream=[]
        self.linkedListDownStream=[]
        

    def set_graph(self,graph):
        self.graph=graph
        return self

    def get_index(self):
        return list(self.graph.vertexList).index(self.name)

    def isUpTerminal(self):
        return len(self.linkedListUpStream)==0
        
    def isDownTerminal(self):
        return len(self.linkedListDownStream)==0
        
    def __repr__(self):
        return "<Vertex:{}>".format(self.name)


class Edge(object):
    def __init__(self, pre, post, name=""):
        self.name = name
        self.pre = pre
        self.post = post

    def toString(self):
        return self.pre.name+"->"+self.post.name

    def set_graph(self,graph):
        self.graph=graph
        return self
    
    def __repr__(self):
        if self.name != "":
            return "<Edge:{}>".format(self.name)
        else:
            return "<Edge:({}->{})>".format(self.pre.name, self.post.name)

class Fluid(object):
    def __init__(self,name,temperature=20):
        self.name=name
        self.temperature=temperature # セルシウス度
        self._rho=None #kg/m3
        self._mu=None #Pa・s=kg/ms
        self._nu=None #m2/s
    @property
    def rho(self):
        if self._rho is not None:
            return self._rho
        else:
            raise
    @property
    def mu(self):
        if self._mu is not None:
            return self._mu
        else:
            raise
    @property
    def nu(self):
        if self._nu is not None:
            return self._nu
        else:
            self._nu=self.mu/self.rho
            return self._nu
    def __repr__(self):
        return "<Fluid:{} temp={:.0f}>".format(self.name,self.temperature)

class WaterFluid(Fluid):
    @property
    def rho(self):
        if self._rho is not None:
            return self._rho
        else:
            t=self.temperature
            self._rho=(999.83952+16.945176*t-7.9870401*1e-3*t**2\
                      -46.170461*1e-6*t**3+105.56302*1e-9*t**4-280.54253*1e-12*t**5)\
                      /(1+16.879850*1e-3*t)
            return self._rho
    @property
    def mu(self):
        if self._mu is not None:
            return self._mu
        else:
            T=273.15+self.temperature
            self._mu=2.414*1e-5*np.power(10,247.8/(T-140))
            return self._mu
    def __repr__(self):
        return "<WaterFluid:{} temp={:.0f}>".format(self.name,self.temperature)

class Pipe(object):
    def __init__(self,name,epsilon):
        self.name=name
        self.epsilon=epsilon

class Circuit(object):
    def __init__(self,graph):
        self.graph=graph
        self.sectionList=OrderedDict()
        self.terminalList=OrderedDict()

    def add_section(self,section):
        self.sectionList[section.edge.toString()]=section
        section.circuit=self
        return section
        
    def add_section2(self,edgeStr,D,epsilon):
        edge=self.graph.add_edge2(edgeStr)
        section=Section(edge,D,epsilon)
        return self.add_section(section)
    def add_sections(self,sections):
        return [self.add_section(s) for s in sections]
    def add_sections2(self,sectionsStringList):
        svl=[]
        for i in sectionsStringList:
            l=i.split(",")
            edgeStr=l[0]
            D=float(l[1])
            epsilon=float(l[2])
            svl.append([edgeStr,D,epsilon])
        return [self.add_section2(edgeStr,D,epsilon) for edgeStr,D,epsilon in svl]

    def add_terminal(self,terminal):
        self.terminalList[terminal.vertex.name]=terminal
        terminal.circuit=self
        return terminal
    def add_terminal2(self,vertexStr,flux):
        v=self.graph.add_vertex2(vertexStr)
        t=Terminal(v,flux)
        return self.add_terminal(t)

    def set_fluid(self,fluid,indexes=None):
        if indexes is None:
            [self.sectionList[sStr].set_fluid(fluid) for sStr in self.sectionList]
        else:
            [self.sectionList[sStr].set_fluid(fluid) for sStr in indexes]

    def set_pipe(self,pipe,indexes=None):
        if indexes is None:
            [self.sectionList[sStr].set_pipe(pipe) for sStr in self.sectionList]
        else:
            [self.sectionList[sStr].set_pipe(pipe) for sStr in indexes]
    def calc_flux(self,fluxDirection=-1):
        #fluxDirection=-1のときTerminalで吐出方向を正としてfluxを読む
        adj=np.array(self.graph.edgeAdjunction)
        indefiniteList=[]
        definiteList=np.zeros(len(self.graph.vertexIndexes))
        for tStr in self.terminalList:
            t=self.terminalList[tStr]
            if t.isDefinite():
                definiteList[self.graph.vertexIndexes.index(t.vertex.name)]=t.flux*fluxDirection
            elif not t.isDefinite():
                indefiniteList.append(self.graph.vertexIndexes.index(t.vertex.name))
        self.fluxA=np.delete(adj,indefiniteList,axis=0)
        self.fluxb=np.delete(definiteList,indefiniteList,axis=0)
        self.flux=np.linalg.solve(self.fluxA,self.fluxb)
        for sStr,flux in zip(self.graph.edgeIndexes,self.flux):
            self.sectionList[sStr].flux=flux

class Terminal(object):
    def __init__(self,vertex,flux):
        self.vertex=vertex
        vertex.terminal=self
        self.flux=flux
    def isDefinite(self):
        return True
    def __repr__(self):
        return "<Treminal: {} : Q={:.3f}".format(self.vertex.name,self.flux)

class IndefiniteTerminal(Terminal):
    def __init__(self,vertex):
        self.vertex=vertex
        vertex.terminal=self
    def isDefinite(self):
        return False
    def __repr__(self):
        return "<Treminal: {} : Q={:.3f}".format(self.vertex.name,self.flux)

        
def calc_Colebrook(f, epsilon, D, Re):
    return 1.0/np.sqrt(f)+2.*np.log10(epsilon/D/3.7+2.51/Re/np.sqrt(f))


def calc_Reynolds_pipe(v, D, nu):
    return v*D/nu


class Section(object):
    def __init__(self,edge,D,epsilon):
        self.edge=edge
        edge.section=self
        self.D=D
        self.epsilon=epsilon
        self.flux=None
    def set_fluid(self,fluid):
        self.fluid=fluid
        return self
    def set_pipe(self,pipe):
        self.pipe=pipe
        return self
    def set_partsEnumerator(self,enum):
        self.partsEnumerator=enum
        enum.section=self
        return self
    def calc_Q(self,fluxUnit=1):
        self.Q=self.flux*fluxUnit
        return self
    def calc_Reynolds(self):
        self.A=self.D**2*np.pi/4.
        self.v=self.Q/self.A
        self.Re=calc_Reynolds_pipe(self.v,self.D,self.fluid.nu)
        return self
    def calc_f(self):
        func=partial(calc_Colebrook,epsilon=self.pipe.epsilon,D=self.D,Re=self.Re)
        self.resultf=root(func,0.01)
        self.f=self.resultf.x[0]
        return self
    def calc_unitPr(self):
        self.Pv=self.v**2*self.fluid.rho/2.
        unitPr=self.f/self.D*self.Pv
        self.unitPr=unitPr
        return self
    def calc_Pr(self):
        self.totalPr=self.partsEnumerator.calc_Pr()
        return self

    def __repr__(self):
        if self.flux is None:
            return "<Section: {} :D={:.3f}>".format(self.edge.toString(),self.D)
        else:
            return "<Section: {} :D={:.3f} :Q={:.3f}>".format(self.edge.toString(),self.D,self.flux)

class PartsEnumerator(object):
    pass

class PartsEnumeratorList(PartsEnumerator):
    def __init__(self):
        self.partsList=[]
    def add_part(self,part):
        self.partsList.append(part)
        part.section=self.section
    def calc_Pr(self):
        self.PrList=[part.Pr for part in self.partsList]
        return np.sum(self.PrList)

class PartsEnumeratorQuantity(PartsEnumerator):
    def __init__(self):
        self.partsDict=OrderedDict()
    def add_part(self,part,quantity):
        self.partsDict[part.name]=[part,quantity]
        part.section=self.section
    def calc_Pr(self):
        self.PrList=[]
        for name in self.partsDict:
            valueList=self.partsDict[name]
            part=valueList[0]
            quantity=valueList[1]
            self.PrList.append(part.Pr*quantity)
        return np.sum(self.PrList)

class Part(object):
    def __init__(self,name=""):
        self.name=name
    def set_section(self,section):
        self.section=section
    @property
    def length(self):
        raise
    @property
    def Pr(self):
        raise

class PartPipe(Part):
    def __init__(self,name,l):
        self.name=name
        self.l=l
    @property
    def length(self):
        return self.l
    @property
    def Pr(self):
        return self.length*self.section.unitPr

class PartJointZeta(Part):
    def __init__(self,name,zeta):
        self.name=name
        self.zeta=zeta
    @property
    def length(self):
        return self.Pr/self.section.unitPr
    @property
    def Pr(self):
        return self.section.Pv*self.zeta

class PartJointEffectiveLength(Part):
    def __init__(self,name,l):
        self.name=name
        self.l=l
    @property
    def length(self):
        return self.l
    @property
    def Pr(self):
        return self.section.unitPr*self.l


if __name__=="__main__":
    g=Graph("g1")
    g.add_vertexes2("A,B,C,D,E,F,G")
    #g.add_edges2("A->B,B->C,C->D,D->E,E->F,F->G")
    c=Circuit(g)
    c.add_sections2(["A->B,0.1,0.00001",
                     "B->C,0.1,0.00001",
                     "C->D,0.1,0.00001",
                     "D->E,0.1,0.00001",
                     "C->F,0.1,0.00001",
                     "F->G,0.1,0.00001"])
    c.add_terminal(Terminal(c.graph.vertexList["E"],300))
    c.add_terminal(Terminal(c.graph.vertexList["G"],400))
    c.add_terminal(IndefiniteTerminal(c.graph.vertexList["A"]))
    c.graph.set_index()
    c.graph.calc_linkedList()
    c.graph.calc_vertexAdjunction()
    c.graph.calc_edgeAdjunction()
    c.calc_flux(fluxDirection=-1)
    c.set_fluid(WaterFluid("Water20",temperature=20))
    c.set_pipe(Pipe("Su",epsilon=3e-6))
    for sStr in c.sectionList:
        s=c.sectionList[sStr]
        s.calc_Q(0.001/60.)
        s.calc_Reynolds()
        s.calc_f()
        s.calc_unitPr()
