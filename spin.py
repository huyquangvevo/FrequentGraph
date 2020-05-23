from utils import readGraphs
from algorithms import GraphCollections


graphs = readGraphs('mico.outx')


if __name__ == "__main__":
    # frequents = getFrequentEdges(graphs,0.8)
    # print(frequents)
    # frequentGraph(graphs,frequents)
    graphDB = GraphCollections(graphs,0.8)
    freqGraphs = graphDB.frequentGraph()
    # print(freqGraphs)
    
