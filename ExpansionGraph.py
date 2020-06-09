import numpy as np
from utils import plotGraph
from algorithm import canonicalForm,string2matrix


class ExpansionGraph():
    def __init__(self,matrixAdj_ : np.ndarray,topoGraphs_,graphs_,freqEdges_,theta_):
        self.graphs = graphs_
        self.theta = theta_
        self.matrixAdj = matrixAdj_
        self.subGraphs = topoGraphs_
        self.spaceGraphs = {np.array2string(matrixAdj_):topoGraphs_}
        self.canEdges = []
        self.associativeEdges = []
        self.setCandidateEdges(freqEdges_)
        self.setAssociativeEdge()


    def setCandidateEdges(self,freqEdges):
        mapEdges = {(e[0],e[1]):e[2] for e in freqEdges.keys()}
        indices = np.where(self.matrixAdj == 0)
        canEdges = []
        for i in range(len(indices[0])):
            iR = indices[0][i]
            iC = indices[1][i]
            k = tuple(sorted((self.matrixAdj[iR,iR],self.matrixAdj[iC,iC])))
            if k in mapEdges and iR <= iC:
                canEdges.append((iR,iC,mapEdges[k]))   
        self.canEdges = canEdges

    def setAssociativeEdge(self):
        for edge in self.canEdges:
            isAssociative = True
            for graph in self.subGraphs.keys():
                for sub in self.subGraphs[graph]:
                    if self.graphs[graph][sub[edge[0]],sub[edge[1]]] != edge[2]:
                        isAssociative = False
                        break
                if not isAssociative:
                    break
            if isAssociative:
                self.associativeEdges.append(edge)


    def joinEdge(self,graph: np.ndarray,edge):
        graph[edge[0],edge[1]] = edge[2]
        graph[edge[1],edge[0]] = edge[2]
        return graph

    def searchGraph(self,graph,canEdges):
        newTempGrapsearchGraphhs = {}
        encodeGraph = np.array2string(graph)

        #bottom-up pruning
        codeFullGraph = self.mergeToGraph(graph,canEdges)
        if codeFullGraph in self.spaceGraphs:
            if len(self.spaceGraphs[codeFullGraph].items()) >= self.theta:
                print("bottom-up aval",codeFullGraph)
                return {
                    codeFullGraph : self.spaceGraphs[codeFullGraph]
                }


        #end bottom-up  
        for i,edge in enumerate(canEdges):
            canGraph = self.joinEdge(graph.copy(),edge)
            embedCanGraph = np.array2string(canGraph)
            for j in self.spaceGraphs[encodeGraph].keys():
                topo = []
                for subGraph in self.spaceGraphs[encodeGraph][j]:
                    sNode = subGraph[edge[0]] # id source node
                    dNode = subGraph[edge[1]] # id destination node
                    if self.graphs[j][sNode,dNode] == edge[2]:
                        topo.append(subGraph)
                if len(topo) > 0:
                    if embedCanGraph not in self.spaceGraphs:
                        self.spaceGraphs[embedCanGraph] = {}
                    self.spaceGraphs[embedCanGraph][j] = topo
            if embedCanGraph in self.spaceGraphs:
                self.searchGraph(canGraph,canEdges[i+1:]) 
            else:
                self.searchGraph(graph,canEdges[i+1:])

        return

    def mergeToGraph(self,graph,canEdges):
        encodeGraph = np.array2string(graph)
        fullGraph = graph.copy()
        for i,edge in enumerate(canEdges):
            fullGraph = self.joinEdge(fullGraph,edge)

        codeFullGraph = np.array2string(fullGraph)
        for idGraph in self.spaceGraphs[encodeGraph].keys():
            topo = []
            for sub in self.spaceGraphs[encodeGraph][idGraph]:
                subGraph = sub.copy()
                flag = True
                for i,edge in enumerate(canEdges):
                    if  self.graphs[idGraph][subGraph[edge[0]],subGraph[edge[1]]] != edge[2]:
                        flag = False
                        break
                if flag:
                    topo.append(subGraph)
            if len(topo) > 0:
                if codeFullGraph not in self.spaceGraphs:
                    self.spaceGraphs[codeFullGraph] = {}
                self.spaceGraphs[codeFullGraph][idGraph] = topo
        return codeFullGraph

    def checkLethal(self):
        initialTree = self.matrixAdj.copy()
        for asEdge in self.associativeEdges:
            self.matrixAdj = self.joinEdge(self.matrixAdj,asEdge)
        
        if canonicalForm(initialTree)['code'] != canonicalForm(self.matrixAdj)['code']:
            return True
        
        self.mergeToGraph(initialTree,self.associativeEdges)
        return False 

    def eliminateAssEdges(self):
        newCans = []
        for edge in self.canEdges:
            if edge not in self.associativeEdges:
                newCans.append(edge)
        self.canEdges = newCans

    def expand(self):
        if self.checkLethal():
            return {}
        
        self.eliminateAssEdges()

        self.searchGraph(self.matrixAdj,self.canEdges)

        frequents = {}
        for k,v in self.spaceGraphs.items():
            if len(v.items()) >= self.theta:
                frequents[k] = v
        
        eqGraphClasses = {}
        canTree = canonicalForm(self.matrixAdj)['code']    
        if len(frequents.items()) > 0:
            for k,v in frequents.items():
                subGraph = string2matrix(k)
                cam = canonicalForm(subGraph)
                if cam['code'] == canTree:
                    eqGraphClasses[k] = v
        # print("eqGraphClass",eqGraphClasses)

        return eqGraphClasses
        



