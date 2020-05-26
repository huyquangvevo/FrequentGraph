import numpy as np
from utils import plotGraph

class ExpansionGraph():
    def __init__(self,matrixAdj_ : np.ndarray,topoGraphs_,graphs_,freqEdges_,theta_):
        self.graphs = graphs_
        self.theta = theta_
        self.matrixAdj = matrixAdj_
        self.spaceGraphs = {np.array2string(matrixAdj_):topoGraphs_}
        self.canEdges = self.getCandidateEdges(freqEdges_)
        # print("canEdges",self.canEdges)

    def getCandidateEdges(self,freqEdges):
        mapEdges = {(e[0],e[1]):e[2] for e in freqEdges.keys()}
        # print("mapEdges",mapEdges)
        # mapEdges = {(15,30) : 10}
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
        # groupFreqGraphs = []
        encodeGraph = np.array2string(graph)
        # print("in search graph",graph)
        # print("can edges",canEdges)
        # print("len spaces graph",len(self.spaceGraphs.items()))
        # print("spaces graph",self.spaceGraphs.keys())
        # for k,v in self.spaceGraphs.items():
                # print(v[0][0])
            # plotGraph(v[0][0],isShowedID=False)
            # break
        # print("search graphs",canEdges)
        # print("space init graph",self.spaceGraphs)
        for i,edge in enumerate(canEdges):
            canGraph = self.joinEdge(np.ndarray.copy(graph),edge)
            embedCanGraph = np.array2string(canGraph)
            for j in self.spaceGraphs[encodeGraph].keys():
                topo = []
                for subGraph in self.spaceGraphs[encodeGraph][j]:
                    sNode = subGraph[edge[0],edge[0]] # id source node
                    dNode = subGraph[edge[1],edge[1]] # id destination node
                    # print("sNode",sNode,"dNode",dNode,"graph shape",self.graphs[j].shape)
                    if self.graphs[j][sNode,dNode] == edge[2]:
                        topo.append(self.joinEdge(subGraph,edge))
                if len(topo) > 0:
                    if embedCanGraph not in self.spaceGraphs:
                        self.spaceGraphs[embedCanGraph] = {}
                    self.spaceGraphs[embedCanGraph][j] = topo
            # print("tree",self.matrixAdj)
            # print("spaceGraphs\n",self.spaceGraphs)
            # for k,v in self.spaceGraphs.items():
                # print(v[0][0])
                # plotGraph(v[0][0],isShowedID=False)
            # exit(0)
            # if embedCanGraph in self.spaceGraphs:
            # groupFreqGraphs.extend(
                # self.searchGraph(canGraph,canEdges[i+1:]) if (embedCanGraph in self.spaceGraphs) else self.searchGraph(graph,canEdges[i+1:])
            # )
            # else:
                # print("not in")
                # groupFreqGraphs.extend(
                    # print(graph)
                    # self.searchGraph(graph,canEdges[i+1:])
                # )
            # if len(groupFreqGraphs) != 0:
                # print("groupFreqGraphs in loop",groupFreqGraphs)
                # exit(0)
        # if len(groupFreqGraphs) != 0:
            self.searchGraph(canGraph,canEdges[i+1:]) if (embedCanGraph in self.spaceGraphs) else self.searchGraph(graph,canEdges[i+1:])
            
        # print("groupFreqGraphs",len(groupFreqGraphs))
        # return groupFreqGraphs
        return
    
    def expand(self):
        # print("canEdges",self.canEdges)
        self.searchGraph(self.matrixAdj,self.canEdges)
        frequents = {}
        for k,v in self.spaceGraphs.items():
            if len(v.items()) > 9 and len(self.spaceGraphs.items()) > 20:
                frequents[k] = v
                break
        if len(frequents.items()) > 0:
            print("result",frequents)
            for k,v in frequents.items():
                for idGraph,g in v.items():
                    print("g",g[0])
                    for i in range(g[0].shape[0]):
                        g[0][i,i] = self.graphs[idGraph][g[0][i,i],g[0][i,i]]
                    print("g after",g[0])
                    plotGraph(g[0],isShowedID=False)
            exit(0)
        # print("S",S)
        return
        



