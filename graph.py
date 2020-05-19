import numpy as np

class Graph:
    def __init__(self,numVertices_):
        self.edges = np.zeros((numVertices_,numVertices_))
