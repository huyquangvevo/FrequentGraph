import numpy as np
from typing import List
import operator

class GraphCollections():
    tempTrees = {}
    def __init__(self,graphs_,theta_):
        self.graphs = graphs_
        self.theta = theta_
        self.freqEdges = self.getFrequentEdges(self.graphs,self.theta)
        self.initTempTree()

    def getFrequentEdges(self,graphs : List[np.ndarray],theta):
        frequentEdges = {}
        for idGraph,graph in enumerate(graphs):
            visited = [False]*len(graph)
            edgesSet = set()
            queue = []
            start = 0
            queue.append(start)
            visited[start] = True
            while queue:
                s = queue.pop(0)
                for i,v in enumerate(graph[s]):
                    if i != s and v > 0:
                        if visited[i] == False:
                            # evaluating edge
                            labelNodes = [graph[s,s],graph[i,i]]
                            labelNodes = sorted(labelNodes)#,reverse=True)
                            encodeEdges = (labelNodes[0],labelNodes[1],v)
                            if encodeEdges not in edgesSet:
                                if encodeEdges not in frequentEdges:
                                    frequentEdges[encodeEdges] = {}
                                    frequentEdges[encodeEdges]['freq'] = 1
                                    frequentEdges[encodeEdges]['edges'] = {}
                                else:
                                    frequentEdges[encodeEdges]['freq'] += 1
                                # frequentEdges[encodeEdges]['freq'] = 0 if encodeEdges not in frequentEdges else frequentEdges[encodeEdges]['freq'] + 1                      
                                edgesSet.add(encodeEdges)
                                frequentEdges[encodeEdges]['edges'][idGraph] = [(s,i) if graph[s,s] == labelNodes[0] else (i,s)]
                            else:
                                frequentEdges[encodeEdges]['edges'][idGraph].append((s,i) if graph[s,s] == labelNodes[0] else (i,s)) 
                            # end evaluating
                            queue.append(i)
                            visited[i] = True

        # tempFrequents = [{k: v['edges']} for k, v in frequentEdges.items() if v['freq'] >theta*len(graphs)]
        frequents = {}
        for k,v in frequentEdges.items():
            if v['freq'] > theta*len(graphs):
                frequents[k] = v['edges']
        return frequents

    def initTempTree(self):
        for edge,matches in self.freqEdges.items():
            # print("edge",edge)
            # matrix = self.canonicalForm(np.array([[edge[0],edge[2]],[edge[2],edge[1]]]))
            matrix = np.array([[edge[0],edge[2]],[edge[2],edge[1]]])
            # print("matrix",matrix)
            encodeEdge = np.array2string(matrix)
            self.tempTrees[encodeEdge] = {}
            for i in matches.keys():# range(len(self.graphs)):
                topo = []
                for e in matches[i]:
                    m = np.array([[e[0],edge[2]],[edge[2],e[1]]])
                    topo.append(m)
                self.tempTrees[encodeEdge][i] = topo  

    def encodeGraph(self,graph):
        visited = [False]*len(graph)
        queue = []
        queue.append(0)
        visited[0] = True
        code = str(graph[0,0]) + '$'
        while queue:
            s = queue.pop(0)
            levelStr = ''
            for i in np.where(graph[s]>0)[0][s+1:]:
                queue.append(i)
                levelStr += str(graph[s,i]) + "_" + str(graph[i,i]) + "_"
                visited[i] = True 
            if levelStr != '':
                code += levelStr[:-1] +  '$'
        code += '#'

        return code


    def canonicalForm(self,graph: np.ndarray):
        labelNodes = graph.diagonal()
        start = np.zeros((1,1),dtype=int)
        maxNodes = np.where(labelNodes == np.max(labelNodes))[0]
        start[0,0] = np.max(labelNodes)
        canonical = {
            "code" : ''
        }
        for idStart in maxNodes:
            S = {
                "tree" : start,
                "index" : np.array([idStart]),
                "code" : self.encodeGraph(start)
            }

            while (len(S["index"]) < len(labelNodes)):
                # trees = []
                newCandidates = {}
                for i in range(graph.shape[0]):
                    if i in S["index"]:
                        continue
                    Q = []
                    t = S["tree"]
                    for id,j in enumerate(S["index"]):
                        if graph[i,j] == 0:
                            continue
                        rowExpand = np.zeros((1,t.shape[0]),dtype=int)
                        rowExpand[0,id] = graph[i,j]
                        tree = np.r_[t,rowExpand]
                        colExpand = np.zeros((tree.shape[0],1),dtype=int)
                        colExpand[id,0] = graph[i,j]
                        colExpand[tree.shape[0]-1,0] = graph[i,i]
                        tree = np.c_[tree,colExpand]
                        indexTree = np.concatenate([S["index"],np.array([i])])
                        codeTree = self.encodeGraph(tree)
                        newCandidates[codeTree] = {
                            "tree" : tree,
                            "index" : indexTree,
                            "code" : codeTree
                        }

                S = newCandidates[max(newCandidates.keys())]
            canonical = S if canonical["code"] < S["code"] else canonical 
        return canonical["tree"]

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


    def exploreFFSM(self,C):
        print("C\n",C)
        if len(C) == 0:
            return
        # newTempTrees = {}
        Q = []
        for X in C:
            S = []
            newTempTrees = {}
            for Y in C:
                candidates = []
                print("X\n",X)
                print("Y\n",Y)
                if np.array_equal(X[:-1,:-1],Y[:-1,:-1]) and not np.array_equal(X[-1],Y[-1]):
                    candidates.append(self.joinCase3bFFSM(X,Y))
                    # print("join",candidates[0],"x",X,"y",Y)
                extensions = self.extendFFSM(X,Y)
                if len(extensions) > 0:
                    print("extension X\n",X,"\nY",Y )
                    candidates.extend(extensions)
                # print("X",X)
                # print("Y",Y)
                print("candidates",candidates)
                for joinedTree in candidates:
                    # joinedTree = self.joinCase3bFFSM(X,Y)
                    print("joinedTree",joinedTree)
                    indexAddNode = np.where(joinedTree[-1] > 0)[0][0]
                    # print("index add node",indexAddNode)
                    embedJoinedTree = np.array2string(joinedTree)
                    for i in self.tempTrees[np.array2string(X)].keys():
                        topo= []
                        for subGraph in self.tempTrees[np.array2string(X)][i]:
                            # print("subgraph",subGraph)
                            linkedNode = subGraph[indexAddNode,indexAddNode] # node is extended
                            # print("indexNode",indexAddNode,"linkNoded",linkedNode)
                            # print("joinedTree",joinedTree)
                            for j in np.where(self.graphs[i][linkedNode] > 0)[0]: # get neighbor of linked node
                                if self.graphs[i][linkedNode,j] == joinedTree[-1,indexAddNode] and j not in subGraph.diagonal():
                                    pad = np.zeros((1,subGraph.shape[0]+1),dtype=int)
                                    pad[0,indexAddNode] = joinedTree[-1,indexAddNode]
                                    pad[0,-1] = j
                                    topo.append(self.extend(subGraph,pad))
                        # print("i",i,"topo",topo)
                        if len(topo) > 0:
                            # print("embed",embedJoinedTree)
                            if embedJoinedTree not in newTempTrees:
                                newTempTrees[embedJoinedTree] = {}
                                S.append(joinedTree)
                            newTempTrees[embedJoinedTree][i] = topo
 
            
            temp = {}
            nextCans = []
            i = 0
            print("new tempTrees",newTempTrees)
            # print("theta",self.theta," len",len(self.graphs))
            print("threshold",self.theta*len(self.graphs))
            for k,v in newTempTrees.items():
                # print("len items",len(v.items()))
                if len(v.items()) > self.theta*len(self.graphs):
                    temp[k] = v
                    nextCans.append(S[i])
                    # print("S i",S[i])
                i += 1
            if len(temp.items()) == 0:
                return
            self.tempTrees = temp
            # print("new tempTrees",self.tempTrees)
            # print("next candidate",nextCans)
            # Q = Q.extend(self.exploreFFSM(nextCans))
            self.exploreFFSM(nextCans)
            return 0

    def freqEdges2matrix(self):
        matrices = []
        for edge in self.freqEdges.keys():
            matrices.append(
                # self.canonicalForm(np.array([[edge[0],edge[2]],[edge[2],edge[1]]]))
                np.array([[edge[0],edge[2]],[edge[2],edge[1]]])
            )
        return matrices


    def exploreGenericTree(self,C,R):
        print("C\n",C)
        Q = []
        for X in C:
            S = []
            newTempTrees = {}
            for Y in C:
                candidates = []
                if np.array_equal(X[:-1,:-1],Y[:-1,:-1]) and not np.array_equal(X[-1],Y[-1]):
                    candidates.append(self.joinCase3bFFSM(X,Y))
                    # print("join",candidates[0],"x",X,"y",Y)
                extensions = self.extendFFSM(X,Y)
                if len(extensions) > 0:
                    # print("extension X\n",X,"\nY",Y )
                    candidates.extend(extensions)
                # print("X",X)
                # print("Y",Y)
                # print("candidates",candidates)
                for joinedTree in candidates:
                    # joinedTree = self.joinCase3bFFSM(X,Y)
                    indexAddNode = np.where(joinedTree[-1] > 0)[0][0]
                    embedJoinedTree = np.array2string(joinedTree)
                    S.append(joinedTree)
                    for i in self.tempTrees[np.array2string(X)].keys():
                        topo= []
                        for subGraph in self.tempTrees[np.array2string(X)][i]:
                            linkedNode = subGraph[indexAddNode,indexAddNode] # node is extended
                            for j in np.where(self.graphs[i][linkedNode] > 0)[0]: # get neighbor of linked node
                                if self.graphs[i][linkedNode,j] == joinedTree[-1,-1] and j not in subGraph.diagonal():
                                    pad = np.zeros((1,subGraph.shape[0]+1),dtype=int)
                                    pad[0,indexAddNode] = joinedTree[-1,-1]
                                    pad[0,-1] = j
                                    topo.append(self.extend(subGraph,pad))

                        if len(topo) > 0:
                            if embedJoinedTree not in newTempTrees:
                                newTempTrees[embedJoinedTree] = {}
                            newTempTrees[embedJoinedTree][i] = topo 
            
            temp = {}
            nextCans = []
            i = 0
            for k,v in newTempTrees.items():
                if len(v.items()) > self.theta*len(self.graphs):
                    temp[k] = v
                    nextCans.append(S[i])
                i += 1

            # S = S - R
            S = []
            # encodeReds = dict((np.array2string(k),1) for k in R)
            for can in nextCans:
                # if np.array2string(can) not in encodeReds:
                if can not in R:
                    S.append(can)
            U,V = self.exploreGenericTree(S,R)
            # encodeQ = dict((np.array2tring(k),1) for k in Q)
            # Q = Q union U union Expansion(X)
            for u in U:
                if u not in Q:
                    Q.append(u)
            
            # R = R union {X} union V 
            if X not in R:
                R.append(X)
            for v in V:
                if v not in R:
                    R.append(v)

        return Q,R 

                

            # self.tempTrees = temp
            # print(S)
            # print("next candidate",nextCans)
            # Q = Q.extend(self.exploreFFSM(nextCans))
    

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
        self.exploreFFSM(self.freqEdges2matrix())
        # print("ca",self.canonicalForm(graphDemo))
        # print(self.extendFFSM(graphDemo,graphDemo2))
        # print(self.joinCase3bFFSM(graphDemo,graphDemo2))
        # print(self.extend(graphDemo,np.array([[0,0,0,0,1]])))
        # canonicalForm(graphDemo)
        return True
                




        






