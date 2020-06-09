from utils import readGraphs,plotGraph
from GraphCollection import GraphCollection
from algorithm import string2matrix
import numpy as np
import json 

datasets = "mico"
graphs = readGraphs('{}.outx'.format(datasets))

def extractResultGraph(results):
    for k,v in results.items():
        numNodeGraphs = np.array([string2matrix(k).shape[0] for k,v in results.items()])
        indicesFreq = np.where(numNodeGraphs == numNodeGraphs.max())[0]
        return [string2matrix(list(results.keys())[i]) for i in indicesFreq]



if __name__ == "__main__":
    
    graphDB = GraphCollection(graphs,1.0)
    print("Frequent edges",len(graphDB.freqEdges.items()))
    
    freqGraphs = graphDB.frequentGraph()
    print("freqGraphs",freqGraphs)
    for freqGraph in freqGraphs:
        plotGraph(freqGraph,False)
