from utils import readGraphs,plotGraph
from algorithms import GraphCollections


graphs = readGraphs('mico.outx')


if __name__ == "__main__":
    # frequents = getFrequentEdges(graphs,0.8)
    # print(frequents)
    # frequentGraph(graphs,frequents)
    graphDB = GraphCollections(graphs,0.8)
    # plotGraph(graphDB.graphs[0])
    freqGraphs = graphDB.frequentGraph()
    # print(graphDB.tempTrees)
    for k,v in graphDB.tempTrees.items():
        for kk,vv in v.items():
            print("graph",kk)
            plotGraph(vv[0],False)
    # print(freqGraphs)
    
