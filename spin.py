from utils import readGraphs,plotGraph
from GraphCollection import GraphCollection
import pandas as pd
import numpy as np


datasets = "mico"
graphs = readGraphs('{}.outx'.format(datasets))


if __name__ == "__main__":
    # frequents = getFrequentEdges(graphs,0.8)
    # print(frequents)
    # frequentGraph(graphs,frequents)
    
    graphDB = GraphCollection(graphs,0.8)
    # plotGraph(graphDB.graphs[0])
    
    freqGraphs = graphDB.frequentGraph()
    print("freqGraphs",freqGraphs)
    df = pd.DataFrame({"MaxSubgraph" : [np.array2string(g) for g in freqGraphs]})
    df.to_csv("result-{}.csv".format(datasets),index=False)
    # exit()
    for freqGraph in freqGraphs:
        plotGraph(freqGraph,False)
   

    # for k,v in freqGraphs:
    #     for kk,vv in v.items():
    #         print("graph",kk)
    #         plotGraph(vv[0],False)
    # print(freqGraphs)
    
