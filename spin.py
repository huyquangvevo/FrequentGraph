from utils import readGraphs,plotGraph
from GraphCollection import GraphCollection
# import pandas as pd
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
    # frequents = getFrequentEdges(graphs,0.8)
    # print(frequents)
    # frequentGraph(graphs,frequents)
    
    graphDB = GraphCollection(graphs,1.0)
    print("Frequent edges",len(graphDB.freqEdges.items()))
    # plotGraph(graphDB.graphs[0])
    
    freqGraphs = graphDB.frequentGraph()
    print("freqGraphs",freqGraphs)
    # with open('result-{}.json'.format(datasets), 'w') as fp:
        # json.dump(freqGraphs, fp)
    # df = pd.DataFrame({"MaxSubgraph" : [np.array2string(g) for g in freqGraphs]})
    # df.to_csv("result-{}.csv".format(datasets),index=False)
    # exit()
    # freqGraphs = extractResultGraph(freqGraphs)
    # print("freqGraph",freqGraph)
    for freqGraph in freqGraphs:
        plotGraph(freqGraph,False)
   

    # for k,v in freqGraphs:
    #     for kk,vv in v.items():
    #         print("graph",kk)
    #         plotGraph(vv[0],False)
    # print(freqGraphs)
    
