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
        # print("theta",self.theta)
        print("freqEdges")#,self.freqEdges)
        for edge in self.freqEdges:
            print(edge)
        print("end freq edges")
        # exit()
        # print("init temp Tree",self.tempTrees)

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
            # print("edge",edge)
            # matrix = self.canonicalForm(np.array([[edge[0],edge[2]],[edge[2],edge[1]]]))
            matrix = np.array([[edge[0],edge[2]],[edge[2],edge[1]]])
            # print("matrix",matrix)
            # matrix = canonicalForm(matrix,code=False)
            # print("after matrix",matrix)
            encodeEdge = np.array2string(matrix)
            self.tempTrees[encodeEdge] = {}
            for i in matches.keys():# range(len(self.graphs)):
                topo = []
                # print("match i ",matches[i])
                for e in matches[i]:
                    m = np.array([e[0],e[1]])
                    # print("i",i,"m",m)
                    topo.append(m)
                self.tempTrees[encodeEdge][i] = topo  

    def joinCase3bFFSM(self,graphX: np.ndarray,graphY: np.ndarray):
        n = graphX.shape[0]
        # if not np.array_equal(graphX[:-1,:-1],graphY[:-1,:-1]):
            # return None
        rowExpand = np.zeros((1,n),dtype=int)
        rowExpand[0] = graphY[n-1]
        rowExpand[0,n-1] = 0
        graphX = np.r_[graphX,rowExpand]
        colExpand = np.concatenate([rowExpand[0],np.array([graphY[n-1,n-1]])])
        colExpand = np.reshape(colExpand,(n+1,1))
        graphX = np.c_[graphX,colExpand]
        # print(graphX)
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
                # indices = np.where(Y[p] > 0)[0]
                for i in indices:
                    # if i!=p:
                    pad = np.zeros((1,n + 1),dtype=int)
                        # pad[0,-1] = Y[i,i]
                        # pad[0,-2] = Y[p,i]
                    pad[0,-1] = Y[p+i+1,p+i+1]
                    pad[0,-2] = Y[p,p+i+1]
                    extensions.append(self.extend(X,pad))
        return extensions


    def hasNoExternalAssEdge(self,tree : np.ndarray,embeddings : dict):
        # return True
        # print("TReeInHasNo",tree)
        numEmb = 0
        for k in embeddings.keys():
            numEmb += len(embeddings[k])
        print("numEmbedding",numEmb)
        # print("embeddings",embeddings)
        for i in range(tree.shape[0]):
            externalEdges = {}
            # print("iExternal",externalEdges)
            for idGraph in embeddings.keys():
                # print("idGRaph",idGraph)
                # print("numSubgraph",len(embeddings[idGraph]))
                # print(embeddings[idGraph])
                for subGraph in embeddings[idGraph]:
                    curNode = subGraph[i]
                    # print("curNode",curNode)
                    indices = np.where(self.graphs[idGraph][curNode] > 0)[0] # neighbor of curNode
                    #indexNodes = subGraph#.diagonal() #[subGraph[j,j] for j in np.where(subGraph[i] > 0)[0]]    # in embedding
                    nodes = list(set(indices) - set(subGraph)) # node not in subgraph of node curNode
                    # print("External-idGraph-subGraph",idGraph,subGraph)
                    # print("incides External",indices)
                    # print("indexNodes Ex",subGraph)
                    # print("nodes External",nodes)
                    # print("curNode",curNode)
                    edges = set()
                    for node in nodes:
                        if node != curNode:
                            # print("edgeToadd",(self.graphs[idGraph][curNode,curNode],self.graphs[idGraph][node,node],self.graphs[idGraph][curNode,node]))
                            edges.add((self.graphs[idGraph][curNode,curNode],self.graphs[idGraph][node,node],self.graphs[idGraph][curNode,node]))
                    # if len(nodes) == len(edges):
                    #     print("numNodes",len(nodes))
                    #     print("numEdgesSet",len(edges))
                    #     print("edgesSet",edges)
                    # exit(0)
                    for edge in edges:
                        try:
                            externalEdges[edge]  += 1
                        except:
                            externalEdges[edge]  = 1
            for k,v in externalEdges.items():
                # print("externalEdges key",k)
                # print("externalEdges value",v)
                if v == numEmb:
                    # exit(0)
                    print("hasExternalAssEdge",tree)
                    return False
        # exit(0)
        print("hasNoAss")
        # print("treeAss",tree)
        # print("AssEmbeding",embeddings)
        return True

    def freqEdges2matrix(self):
        matrices = []
        for edge in self.freqEdges.keys():
            matrices.append(
                # self.canonicalForm(np.array([[edge[0],edge[2]],[edge[2],edge[1]]]))
                np.array([[edge[0],edge[2]],[edge[2],edge[1]]])
            )
        return matrices

    def frequentTrees(self,X,restCan,C):
        S = []
        newTempTrees = {}
        print("timestamp frequent trees",datetime.datetime.now())
        # print("X in",X)
        # print("C X in",C)
        for Y in [string2matrix(k) for k in restCan]:
            candidates = []
            # print("X",X)
            # print("Y",Y)
            if np.array_equal(X[:-1,:-1],Y[:-1,:-1]) and not np.array_equal(X[-1],Y[-1]):
                candidates.append(self.joinCase3bFFSM(X.copy(),Y.copy()))
                candidates.append(self.extendCase3bFFSM(X.copy()))
            # print("XX",X)
            # print("YY",Y)
            # print("can join",candidates)
            

            # if X.shape[0] < 3:
                # extensions = extendOneNode(X,Y) #self.extendFFSM(X,Y)
                # if len(extensions) > 0:
                    # candidates.extend(extensions)
            # print("candidates",len(candidates))
            
            extensions = extendOneNode(X,Y)
            if len(extensions) > 0:
                candidates.extend(extensions)
            
            for joinedTree in candidates:
                # print("preJoinedTreeX",X)
                # print("preJoinedTReeY",Y)
                # print("joinedTree",joinedTree)
                indexAddNode = np.where(joinedTree[-1] > 0)[0][0]
                # print("indexAddNode",np.where(joinedTree[-1] > 0)[0])
                # embedJoinedTree = np.array2string(joinedTree)
                cam = canonicalForm(joinedTree)
                embedCamTree = np.array2string(cam['tree'])
                for i in C[np.array2string(X)].keys():
                    topo= []
                    for subNodes in C[np.array2string(X)][i]:
                        linkedNode = subNodes[indexAddNode]
                        for j in np.where(self.graphs[i][linkedNode] > 0)[0]:
                            if self.graphs[i][linkedNode,j] == joinedTree[-1,indexAddNode] and j not in subNodes:
                                # print("prevNodes",subNodes)
                                newSubNodes = np.append(subNodes,np.array([j]))
                                # print("newSubNodes",newSubNodes)
                                # print("camTree",cam["index"])
                                reindicedNodes = np.array([newSubNodes[idNode] for idNode in cam['index']])
                                # print("newReindexNode",reindicedNodes)
                                topo.append(reindicedNodes)
                
                    # for subGraph in C[np.array2string(X)][i]:
                    #     linkedNode = subGraph[indexAddNode,indexAddNode] # node is extended
                    #     for j in np.where(self.graphs[i][linkedNode] > 0)[0]: # get neighbor of linked node
                    #         if self.graphs[i][linkedNode,j] == joinedTree[-1,indexAddNode] and j not in subGraph.diagonal():
                    #             pad = np.zeros((1,subGraph.shape[0]+1),dtype=int)
                    #             pad[0,indexAddNode] = joinedTree[-1,indexAddNode]
                    #             pad[0,-1] = j
                    #             topo.append(self.extend(subGraph,pad))
                    # cam = self.canonicalForm()
                    if len(topo) > 0:
                        # print("topo",topo)
                        if embedCamTree not in newTempTrees:
                            newTempTrees[embedCamTree] = {}
                            # S.append(joinedTree)
                        newTempTrees[embedCamTree][i] = topo
                # print("CamJoined",cam)
                # print("embedCamTRee",embedCamTree)
            
        temp = {}
        
        for k,v in newTempTrees.items():
            if len(v.items()) >= self.theta:
                # temp[np.array2string(canonicalForm(string2matrix(k),code=False))] = v
                temp[k] = v
        # print(temp)
        return temp

    def frequentTrees2(self,tree: np.ndarray,embeddings: dict):
        n = tree.shape[0]
        canFreqTree = {}
        print("timestamp frequent trees 2",datetime.datetime.now())
        for i in range(n):
            curLabel = tree[i,i]
            canEdges = []
            for edge in self.freqEdges.keys():
                if edge[0] == curLabel:
                    canEdges.append(edge)
                elif edge[1] == curLabel and edge[0] != edge[1]:
                    canEdges.append((edge[1],edge[0],edge[2]))
            # print('curLabel',curLabel)
            # print('canEdgeLabel',canEdges)
            for edge in canEdges:
                pad = np.zeros((1,n+1),dtype=int)
                pad[0,-1] = edge[1]
                pad[0,i] = edge[2]
                canTree = self.extend(tree.copy(),pad)
                embededCanTree = np.array2string(canTree)
                # print("curTree",tree)
                # print("extendTRee",canTree)

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


    def exploreGenericTree(self,C : dict,R : dict):
        print("C in\n",C.keys())
        # print("Temptrees len", len(tempTrees.items()))
        print("timestamp",datetime.datetime.now())
        Q = {}
        # print("reprGroup",C)
        for idCan,can in enumerate(C.items()):
            X = string2matrix(can[0])
            print('begin frequent trees')
            nextCan = {k:C[k] for k in list(C.keys())[idCan+1:]}

            encodeX = np.array2string(X)
            S = self.frequentTrees2(X,C[encodeX])
            print("After S-freq",S.keys())
            
            # S - R
            canonicalS = {canonicalForm(string2matrix(k))['code']:k for k in S.keys()}
            # print("SSS before",S.keys())
            # print("RRBefore",R.keys())
            for kR in R.keys():
                canonicalKR = canonicalForm(string2matrix(kR))['code']
                if canonicalKR in canonicalS.keys():
                    if canonicalS[canonicalKR] in S:
                        # print('del duplicate',canonicalKR)
                        del S[canonicalS[canonicalKR]]

            # remove duplicate
            # for kR in R.keys():
                # if kR in S:
                    # del S[kR]
                    

            # print("tempTree after",self.tempTrees)
            # if len(S) != 0:
            # print("afterFrequentTRee size S: ",S.keys())
            U,V = self.exploreGenericTree(S,R.copy())
            # print("S empty",S)
            # print("R empty",R)
            # print("U empty",U)
            # print("V empty",V)
            # print("X ok ex",X)
            # print("encode X",encodeX)
            for k,v in U.items():
                Q[k] = v

            R[encodeX] = C[encodeX]
            for k,v in V.items():
                R[k] = v
            if self.hasNoExternalAssEdge(X,C[encodeX]):
                print("ok expansion",X)
                # print("topoX",C[encodeX])
                eg = ExpansionGraph(
                    X.copy(),
                    C[encodeX],
                    self.graphs,
                    self.freqEdges,
                    self.theta
                )
                expansionX = eg.expand()
                # print("beforeExpansionX",X)
                # print("expansion X",expansionX)
                for k,v in expansionX.items():
                    Q[k] = v
                # exit(0)

            
            # del eg
        # print("Q",len(Q.items))
        # print("R",R)
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
        # print(self.tempTrees)
        # self.exploreFFSM(self.freqEdges2matrix())
        
        # self.exploreGenericTree(self.freqEdges2matrix(),[],self.tempTrees)
        # print("canonical",self.canonicalForm(graphDemo))
        # exit(0)
        results,S = self.exploreGenericTree(self.tempTrees,{})
        # print("final results",results[0])
        # try:
        print("S-final",S)
        print("final result",results.keys())
        if len(results.items()) == 0:
            return []
        numNodeGraphs = np.array([len(np.where(string2matrix(k) > 0)[0]) for k,v in results.items()])
        # print("numNodeGraphs",numNodeGraphs)
        indicesFreq = np.where(numNodeGraphs == numNodeGraphs.max())[0]
        # return results[0]
        # print("indicesFreq",indicesFreq)
        return [string2matrix(list(results.keys())[i]) for i in indicesFreq]
        # except:
            # return []
        # freqGraphs = [for k in freqGraphs]
        # for k,v in results:
            
        
        # print("result extend")
        # print(extendOneNode(graphDemo,graphDemo2))

        # print(np.where(graphDemo == 0))
        # eg = ExpansionGraph(self.tempTrees)
        # print("ca",self.canonicalForm(graphDemo))
        # print(self.extendFFSM(graphDemo,graphDemo2))
        # print(self.joinCase3bFFSM(graphDemo,graphDemo2))
        # print(self.extend(graphDemo,np.array([[0,0,0,0,1]])))
        # canonicalForm(graphDemo)
        return True
                




        






