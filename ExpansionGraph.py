import numpy as np

class ExpansionGraph():
    def __init__(self,matrixAdj_ : np.ndarray,topoGraphs_,graphs_,freqEdges_,theta_):
        self.graphs = graphs_
        self.theta = theta_
        self.matrixAdj = matrixAdj_
        self.spaceGraphs = {np.array2string(matrixAdj_):topoGraphs_}
        self.canEdges = self.getCandidateEdges(freqEdges_)

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
        groupFreqGraphs = []
        encodeGraph = np.array2string(graph)
        print("search graphs",canEdges)
        for i,edge in enumerate(canEdges):
            canGraph = self.joinEdge(graph,edge)
            embedCanGraph = np.array2string(canGraph)
            for j in self.spaceGraphs[encodeGraph].keys():
                topo = []
                for subGraph in self.spaceGraphs[encodeGraph][j]:
                    sNode = subGraph[edge[0],edge[0]] # id source node
                    dNode = subGraph[edge[1],edge[1]] # id destination node
                    if self.graphs[i][sNode,dNode] == edge[2]:
                        topo.append(self.joinEdge(subGraph,edge))
                if len(topo) > 0:
                    if embedCanGraph not in self.spaceGraphs:
                        self.spaceGraphs[embedCanGraph] = {}
                    self.spaceGraphs[embedCanGraph][j] = topo
            print("spaceGraphs\n",self.spaceGraphs)
            exit(0)
            groupFreqGraphs.append(
                self.searchGraph(canGraph,canEdges[i+1:])
            )
        return groupFreqGraphs
    
    def expand(self):
        # print("canEdges",self.canEdges)
        S = self.searchGraph(self.matrixAdj,self.canEdges)
        return S
        



