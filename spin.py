from utils import readGraphs,plotGraph
from GraphCollection import GraphCollection


graphs = readGraphs('mico-demo.outx')


if __name__ == "__main__":
    # frequents = getFrequentEdges(graphs,0.8)
    # print(frequents)
    # frequentGraph(graphs,frequents)
    
    graphDB = GraphCollection(graphs,0.8)
    # plotGraph(graphDB.graphs[0])
    
    freqGraphs = graphDB.frequentGraph()
    # print(graphDB.tempTrees)
    exit()
    for k,v in graphDB.tempTrees.items():
        for kk,vv in v.items():
            print("graph",kk)
            plotGraph(vv[0],False)
    # print(freqGraphs)
    
