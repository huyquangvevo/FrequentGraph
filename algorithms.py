import numpy as np
from typing import List
import operator


def getFrequentEdges(graphs : List[np.ndarray],theta):
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
                        labelNodes = sorted(labelNodes)
                        # encodeEdges = '{}-{}-{}'.format(labelNodes["tree"],labelNodes["index"],v)
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
                            frequentEdges[encodeEdges]['edges'][idGraph] = [(s,i)]
                        else:
                            frequentEdges[encodeEdges]['edges'][idGraph].append((s,i)) 
                        # end evaluating
                        queue.append(i)
                        visited[i] = True

    frequents = [{k: v['edges']} for k, v in frequentEdges.items() if v['freq'] >theta*len(graphs)]
    return frequents

# def JoinCase3bFFSM(x: np.ndarray,y: np.ndarray):

def encodeGraph(graph):
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


def canonicalForm(graph: np.ndarray):
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
            "code" : encodeGraph(start)
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
                    codeTree = encodeGraph(tree)
                    newCandidates[codeTree] = {
                        "tree" : tree,
                        "index" : indexTree,
                        "code" : codeTree
                    }

            S = newCandidates[max(newCandidates.keys())]
        canonical = S if canonical["code"] < S["code"] else canonical 
    print(canonical)            
    return canonical




def FFSM(edges):
    Q = []
    # for edge in edges:
    return False

def frequentGraph(graphs,edges):
    graphDemo = np.array([
        [2,11,10,11],
        [11,1,0,11],
        [10,0,1,10],
        [11,11,10,1]
    ])
    canonicalForm(graphDemo)
    return True


