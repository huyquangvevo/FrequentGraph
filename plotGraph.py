from utils import readGraphs,plotGraph
from algorithms import GraphCollections


graphs = readGraphs('mico.outx')


if __name__ == "__main__":
    graphDB = GraphCollections(graphs,0.8)
    for graph in graphs:
        plotGraph(graph)
