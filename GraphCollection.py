import numpy as np
from typing import List
import operator
from ExpansionGraph import ExpansionGraph
from algorithm import string2matrix,extendOneNode,canonicalForm
import datetime

class GraphCollection():
    tempTrees = {}
    def __init__(self,graphs_,theta_):
        self.graphs = graphs_
        self.theta = theta_ * len(graphs_)
        self.freqEdges = self.getFrequentEdges(self.graphs,self.theta)
        self.initTempTree()

    def getFrequentEdges(self,graphs : List[np.ndarray],theta):
        frequentEdges = {}
        for idGraph, graph in enumerate(graphs):
            edgesSet = set()
            for i in range(graph.shape[0]):
                indicesEdge = np.where(graph[i,i+1:] > 0)
                for des in indicesEdge[0]:
                    labelNodes = [graph[i,i],graph[i+des+1,i+des+1]]
                    labelNodes = sorted(labelNodes)#,reverse=True)
                    encodeEdges = (labelNodes[0],labelNodes[1],graph[i,i+des+1])
                    if encodeEdges not in edgesSet:
                        if encodeEdges not in frequentEdges:
                            frequentEdges[encodeEdges] = {}
                            frequentEdges[encodeEdges]['freq'] = 1
                            frequentEdges[encodeEdges]['edges'] = {}
                        else:
                            frequentEdges[encodeEdges]['freq'] += 1

                        edgesSet.add(encodeEdges)
                        frequentEdges[encodeEdges]['edges'][idGraph] = [(i,i+des+1) if graph[i,i] == labelNodes[0] else (des + i + 1,i)]
                    else:
                        frequentEdges[encodeEdges]['edges'][idGraph].append((i,i + des + 1) if graph[i,i] == labelNodes[0] else (des + i + 1,i))

        frequents = {}
        for k,v in frequentEdges.items():
            if v['freq'] >= theta:
                frequents[k] = v['edges']
        return frequents


    def initTempTree(self):
        for edge,matches in self.freqEdges.items():
            matrix = np.array([[edge[0],edge[2]],[edge[2],edge[1]]])
            encodeEdge = np.array2string(matrix)
            self.tempTrees[encodeEdge] = {}
            for i in matches.keys():
                topo = []
                for e in matches[i]:
                    m = np.array([e[0],e[1]])
                    topo.append(m)
                self.tempTrees[encodeEdge][i] = topo  

    def joinCase3bFFSM(self,graphX: np.ndarray,graphY: np.ndarray):
        n = graphX.shape[0]
        rowExpand = np.zeros((1,n),dtype=int)
        rowExpand[0] = graphY[n-1]
        rowExpand[0,n-1] = 0
        graphX = np.r_[graphX,rowExpand]
        colExpand = np.concatenate([rowExpand[0],np.array([graphY[n-1,n-1]])])
        colExpand = np.reshape(colExpand,(n+1,1))
        graphX = np.c_[graphX,colExpand]
        return graphX

    def extendCase3bFFSM(self,graphX: np.ndarray):
        n = graphX.shape[0]
        pad = np.zeros((1,n),dtype=int)
        pad[0] = graphX[-1]
        pad[0,n-1] = 0
        graphX = np.r_[graphX,pad]
        pad = np.concatenate([pad[0],np.array([graphX[n-1,n-1]])])
        pad = np.reshape(pad,(n+1,1))
        graphX = np.c_[graphX,pad]
        return graphX

    def extend(self,X: np.ndarray,pad: np.ndarray):
        n = X.shape[0]
        X = np.r_[X,pad[:,:-1]]
        pad = np.reshape(pad,(pad[0].shape[0],1))
        X = np.c_[X,pad]
        return X

    def extendFFSM(self,X: np.ndarray,Y: np.ndarray):
        pos =  np.where(Y.diagonal() == X[-1,-1])[0]
        n = X.shape[0]
        extensions = []
        for p in pos:
            if p<n-1:
                indices = np.where(Y[p,p+1:] > 0)[0]
                for i in indices:
                    pad = np.zeros((1,n + 1),dtype=int)
                    pad[0,-1] = Y[p+i+1,p+i+1]
                    pad[0,-2] = Y[p,p+i+1]
                    extensions.append(self.extend(X,pad))
        return extensions


    def hasNoExternalAssEdge(self,tree : np.ndarray,embeddings : dict):
        numEmb = 0
        for k in embeddings.keys():
            numEmb += len(embeddings[k])
        for i in range(tree.shape[0]):
            externalEdges = {}
            for idGraph in embeddings.keys():
                for subGraph in embeddings[idGraph]:
                    curNode = subGraph[i]
                    indices = np.where(self.graphs[idGraph][curNode] > 0)[0] # neighbor of curNode
                    nodes = list(set(indices) - set(subGraph)) # node not in subgraph of node curNode
                    edges = set()
                    for node in nodes:
                        if node != curNode:
                            edges.add((self.graphs[idGraph][curNode,curNode],self.graphs[idGraph][node,node],self.graphs[idGraph][curNode,node]))
                    for edge in edges:
                        try:
                            externalEdges[edge]  += 1
                        except:
                            externalEdges[edge]  = 1
            for k,v in externalEdges.items():
                if v == numEmb:
                    print("hasExternalAssEdge",tree)
                    return False
        print("hasNoAss")
        return True

    def freqEdges2matrix(self):
        matrices = []
        for edge in self.freqEdges.keys():
            matrices.append(
                np.array([[edge[0],edge[2]],[edge[2],edge[1]]])
            )
        return matrices

    def frequentTrees2(self,tree: np.ndarray,embeddings: dict):
        n = tree.shape[0]
        canFreqTree = {}
        print("frequent trees 2: ",datetime.datetime.now())
        for i in range(n):
            curLabel = tree[i,i]
            canEdges = []
            for edge in self.freqEdges.keys():
                if edge[0] == curLabel:
                    canEdges.append(edge)
                elif edge[1] == curLabel and edge[0] != edge[1]:
                    canEdges.append((edge[1],edge[0],edge[2]))
            for edge in canEdges:
                pad = np.zeros((1,n+1),dtype=int)
                pad[0,-1] = edge[1]
                pad[0,i] = edge[2]
                canTree = self.extend(tree.copy(),pad)
                embededCanTree = np.array2string(canTree)

                for idGraph in embeddings.keys():
                    canNodes = []
                    for subNodes in embeddings[idGraph]:
                        curNode = subNodes[i]
                        indices = np.where(self.graphs[idGraph][curNode] == edge[2])[0] # neighbor of edge
                        for idNeibor in indices:
                            if idNeibor not in subNodes and self.graphs[idGraph][idNeibor,idNeibor] == edge[1]:
                                canNodes.append(np.append(subNodes,np.array([idNeibor])))
                        
                    if len(canNodes) > 0:
                        if embededCanTree not in canFreqTree:
                            canFreqTree[embededCanTree] = {}
                        canFreqTree[embededCanTree][idGraph] = canNodes
        freqTree = {}
        for k,v in canFreqTree.items():
            if len(v.items()) >= self.theta:
                freqTree[k] = v

        return freqTree

    countGen = 0
    def exploreGenericTree(self,C : dict,R : dict):
        print("timestamp",datetime.datetime.now())
        Q = {}
        for idCan,can in enumerate(C.items()):
            X = string2matrix(can[0])
            nextCan = {k:C[k] for k in list(C.keys())[idCan+1:]}

            encodeX = np.array2string(X)
            S = self.frequentTrees2(X,C[encodeX])
            
            # remove duplicate
            # S - R 
            canonicalS = {canonicalForm(string2matrix(k))['code']:k for k in S.keys()}
            for kR in R.keys():
                canonicalKR = canonicalForm(string2matrix(kR))['code']
                if canonicalKR in canonicalS.keys():
                    if canonicalS[canonicalKR] in S:
                        del S[canonicalS[canonicalKR]]

            self.countGen += 1
            print("Loop in generic tree: ",self.countGen)
            U,V = self.exploreGenericTree(S,R.copy())
            for k,v in U.items():
                Q[k] = v

            R[encodeX] = C[encodeX]
            for k,v in V.items():
                R[k] = v
            if self.hasNoExternalAssEdge(X,C[encodeX]):
                print("MaxExpansion",X)
                eg = ExpansionGraph(
                    X.copy(),
                    C[encodeX],
                    self.graphs,
                    self.freqEdges,
                    self.theta
                )
                expansionX = eg.expand()
                for k,v in expansionX.items():
                    Q[k] = v

            
        return Q,R

    def frequentGraph(self):
        graphDemo = np.array([
            [2,11,10,11],
            [11,1,0,11],
            [10,0,1,10],
            [11,11,10,1]
        ])

        graphDemo2 = np.array([
            [2,11,10,15],
            [11,1,0,15],
            [10,0,1,16],
            [15,15,16,3]
        ])
        results,S = self.exploreGenericTree(self.tempTrees,{})
        print("S-final",S)
        print("final result",results.keys())
        if len(results.items()) == 0:
            return []

        with open('result-mico.txt','w') as f:
            for k in results.keys():
                f.write('$' + k.replace('\n',''))
                for g,v in results[k].items():
                    f.write('\n#' + str(g) + '\n')
                    for topo in v:
                        f.write(' '.join(str(x) for x in topo))
                f.write('\n')


        numNodeGraphs = np.array([len(np.where(string2matrix(k) > 0)[0]) for k,v in results.items()])
        indicesFreq = np.where(numNodeGraphs == numNodeGraphs.max())[0]
        return [string2matrix(list(results.keys())[i]) for i in indicesFreq]
                




        






