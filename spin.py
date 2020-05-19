from utils import readGraphs
from algorithms import getFrequentEdges


graphs = readGraphs('mico.outx')


if __name__ == "__main__":
    frequents = getFrequentEdges(graphs,0.8)
    print(frequents)