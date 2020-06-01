import numpy as np
from utils import plotGraph
from algorithm import canonicalForm


class ExpansionGraph():
    def __init__(self,matrixAdj_ : np.ndarray,topoGraphs_,graphs_,freqEdges_,theta_):
        self.graphs = graphs_
        self.theta = theta_
        self.matrixAdj = matrixAdj_
        self.subGraphs = topoGraphs_
        self.spaceGraphs = {np.array2string(matrixAdj_):topoGraphs_}
        self.canEdges = [] #self.setCandidateEdges(freqEdges_)
        self.associativeEdges = [] #self.setAssociativeEdge()
        self.setCandidateEdges(freqEdges_)
        self.setAssociativeEdge()
        print("associate edges",self.associativeEdges)


    def setCandidateEdges(self,freqEdges):
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
        self.canEdges = canEdges
        # return canEdges

    def setAssociativeEdge(self):
        for edge in self.canEdges:
            isAssociative = True
            for graph in self.subGraphs.keys():
                for sub in self.subGraphs[graph]:
                    if self.graphs[graph][sub[edge[0],edge[0]],sub[edge[1],edge[1]]] != edge[2]:
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
            if len(self.spaceGraphs[codeFullGraph].items()) > self.theta:
                print("bottom-up aval",self.spaceGraphs[codeFullGraph])
                # exit(0)
                return {
                    codeFullGraph : self.spaceGraphs[codeFullGraph]
                }


        #end bottom-up  
            

        # print("len canEdges",len(canEdges))
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

        return

    def mergeToGraph(self,graph,canEdges):
        encodeGraph = np.array2string(graph)
        fullGraph = graph.copy()
        for i,edge in enumerate(canEdges):
            fullGraph = self.joinEdge(fullGraph,edge)
        # print("full Graph",fullGraph)
        codeFullGraph = np.array2string(fullGraph)
        for idGraph in self.spaceGraphs[encodeGraph].keys():
            topo = []
            for sub in self.spaceGraphs[encodeGraph][idGraph]:
                subGraph = sub.copy()
                flag = True
                # print("beforeSubgraph",subGraph)
                for i,edge in enumerate(canEdges):
                    # print("edge",edge)
                    # print("subGraphEdge",self.graphs[idGraph][subGraph[edge[0],edge[0]],subGraph[edge[1],edge[1]]])
                    if  self.graphs[idGraph][subGraph[edge[0],edge[0]],subGraph[edge[1],edge[1]]] == edge[2]:
                        subGraph = self.joinEdge(subGraph,edge)
                    else:
                        flag = False
                        break
                # print("subGraph",subGraph,flag)
                if flag:
                    topo.append(subGraph)
            if len(topo) > 0:
                if codeFullGraph not in self.spaceGraphs:
                    self.spaceGraphs[codeFullGraph] = {}
                self.spaceGraphs[codeFullGraph][idGraph] = topo
                # print("codeTopo2",codeFullGraph)
                # print("hasTopo",self.spaceGraphs.keys())
            # print("topo",topo)
            # print("encodeGraph",self.spaceGraphs[encodeGraph])
        # print('inMerge fullGraph',codeFullGraph)
        # print("returnMerge",self.spaceGraphs.keys())
        return codeFullGraph

    def checkLethal(self):
        initialTree = self.matrixAdj.copy()
        for asEdge in self.associativeEdges:
            self.matrixAdj = self.joinEdge(self.matrixAdj,asEdge)
        
        if canonicalForm(initialTree) != canonicalForm(self.matrixAdj):
            return True
        
        # print("spaceGraphsInMerge",self.spaceGraphs.keys())
        # print("matrix",self.matrixAdj)
        # print("ass",self.associativeEdges)
        self.mergeToGraph(initialTree,self.associativeEdges)
        # print("afterMergeMatrix",self.spaceGraphs.keys())
        return False 

        # print("initialTree",initialTree)
        # print("matrix lethal",self.matrixAdj)
        # return True if canonicalForm(initialTree) != canonicalForm(self.matrixAdj) else initialTree

    def eliminateAssEdges(self):
        # self.canEdges = list(set(self.canEdges) - set(self.associativeEdges))
        # print("can edges",self.canEdges)
        # print("associative edges",self.associativeEdges)
        newCans = []
        for edge in self.canEdges:
            if edge not in self.associativeEdges:
                newCans.append(edge)
        self.canEdges = newCans
        # print("after canEdges",self.canEdges)    

    def expand(self):
        # print("canEdges",self.canEdges)
        # print("associative edges",self.associativeEdges)
        # print("isLethal",self.checkLethal())
        # exit(0)

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
        

        if self.checkLethal():
            return {}
        
        # print("before eliminate",self.canEdges)
        self.eliminateAssEdges()
        print("eliminated",self.canEdges)
        # exit(0)

        self.searchGraph(self.matrixAdj,self.canEdges)

        print("space graphs",self.spaceGraphs)
        frequents = {}
        for k,v in self.spaceGraphs.items():
            if len(v.items()) >= self.theta:
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
        



