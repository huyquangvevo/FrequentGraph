import numpy as np
from utils import plotGraph
from algorithm import canonicalForm


class ExpansionGraph():
    def __init__(self,matrixAdj_ : np.ndarray,topoGraphs_,graphs_,freqEdges_,theta_):
        self.graphs = graphs_
        self.theta = theta_
        self.matrixAdj = matrixAdj_
        self.spaceGraphs = {np.array2string(matrixAdj_):topoGraphs_}
        self.canEdges = self.getCandidateEdges(freqEdges_)

    def getCandidateEdges(self,freqEdges):
        mapEdges = {(e[0],e[1]):e[2] for e in freqEdges.keys()}
        # mapEdges = {(15,30) : 10}
        # print("freqEdges",freqEdges)
        indices = np.where(self.matrixAdj == 0)
        # print("matrix",self.matrixAdj)
        # print("indices",indices)
        canEdges = []
        for i in range(len(indices[0])):
            iR = indices[0][i]
            iC = indices[1][i]
            k = tuple(sorted((self.matrixAdj[iR,iR],self.matrixAdj[iC,iC])))
            if k in mapEdges and iR <= iC:
                canEdges.append((iR,iC,mapEdges[k]))   
        # print(canEdges)
        return canEdges

    def joinEdge(self,graph: np.ndarray,edge):
        graph[edge[0],edge[1]] = edge[2]
        graph[edge[1],edge[0]] = edge[2]
        return graph

    def searchGraph(self,graph,canEdges):
        newTempGraphs = {}
        encodeGraph = np.array2string(graph)
        print("len canEdges",len(canEdges))
        for i,edge in enumerate(canEdges):
            canGraph = self.joinEdge(graph.copy(),edge)
            embedCanGraph = np.array2string(canGraph)
            for j in self.spaceGraphs[encodeGraph].keys():
                topo = []
                for subGraph in self.spaceGraphs[encodeGraph][j]:
                    sNode = subGraph[edge[0],edge[0]] # id source node
                    dNode = subGraph[edge[1],edge[1]] # id destination node
                    if self.graphs[j][sNode,dNode] == edge[2]:
                        topo.append(self.joinEdge(subGraph,edge))
                if len(topo) > 0:
                    if embedCanGraph not in self.spaceGraphs:
                        self.spaceGraphs[embedCanGraph] = {}
                    self.spaceGraphs[embedCanGraph][j] = topo
            # self.searchGraph(canGraph,canEdges[i+1:]) if (embedCanGraph in self.spaceGraphs) else self.searchGraph(graph,canEdges[i+1:])
            if embedCanGraph in self.spaceGraphs:
                self.searchGraph(canGraph,canEdges[i+1:]) 
            else:
                self.searchGraph(graph,canEdges[i+1:])

            # print("lenCanEdges in loop",len(canEdges))
            # print("loop in space graphs",i)   
        print("return spaces",len(self.spaceGraphs.items()))
        return

    
    
    def expand(self):
        print("canEdges",len(self.canEdges))
        # for k,v in self.spaceGraphs.items():
            # for i,g in v.items():
                # plotGraph(g[0])

        # for k,v in self.spaceGraphs.items():
        #     for idGraph,g in v.items():
        #         # print("g",g[0])
        #         for i in range(g[0].shape[0]):
        #             g[0][i,i] = self.graphs[idGraph][g[0][i,i],g[0][i,i]]
        #         # print("g after",g[0])
        #         plotGraph(g[0],isShowedID=False)
        
        # self.searchGraph(self.matrixAdj,self.canEdges)
        self.searchGraph(self.matrixAdj,0)

        print("space graphs",self.spaceGraphs)
        frequents = {}
        for k,v in self.spaceGraphs.items():
            if len(v.items()) >= self.theta * len(self.graphs):
                frequents[k] = v
                # break
        
        eqGraphClasses = {}
        canTree = canonicalForm(self.matrixAdj)    
        if len(frequents.items()) > 0:
            # print("result",frequents)
            for k,v in frequents.items():
                idFirstGraph = list(v.keys())[0]
                subGraph = np.copy(v[idFirstGraph][0])
                # subGraph = np.copy(firstGraph[0])
                for i in range(subGraph.shape[0]):
                    subGraph[i,i] = self.graphs[idFirstGraph][subGraph[i,i],subGraph[i,i]]
                # print("subgraph",subGraph)
                # print("can subgraph",canonicalForm(subGraph))
                # print("tree",self.matrixAdj)
                # print("canTree",canTree)
                if canonicalForm(subGraph) == canTree:
                    eqGraphClasses[np.array2string(subGraph)] = v
            # exit(0)
        print("eqGraphClass",eqGraphClasses)
        # print("current form tree",canonicalForm(self.matrixAdj))

        return eqGraphClasses
        



