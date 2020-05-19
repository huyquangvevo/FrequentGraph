import numpy as np
from typing import List

def getFrequentEdges(graphs : List[np.ndarray],theta):
    frequentEdges = {}
    for graph in graphs:
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
                        encodeEdges = '{}-{}-{}'.format(labelNodes[0],labelNodes[1],v)
                        if encodeEdges not in edgesSet:
                            frequentEdges[encodeEdges] = 0 if encodeEdges not in frequentEdges else frequentEdges[encodeEdges] + 1                      
                            edgesSet.add(encodeEdges)
                        # end evaluating
                        queue.append(i)
                        visited[i] = True

    frequents = [tuple([int(x) for x in k.split('-')]) for k, v in frequentEdges.items() if v >theta*len(graphs)]
    return frequents


