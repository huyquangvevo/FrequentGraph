from utils import readGraphs,plotGraph
from GraphCollection import GraphCollection


graphs = readGraphs('mico.outx')


if __name__ == "__main__":
    graphDB = GraphCollection(graphs,0.8)
    print(graphDB.freqEdges)
    # exit(0)
    print("End frequent")
    for graph in graphs:
        plotGraph(graph)
