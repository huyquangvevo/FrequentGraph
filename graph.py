class Vertex:
    def __init__(self,id_,label_):
        self.id = id_
        self.label = label_

class Edge:
    def __init__(self,v1_,v2_,label_):
        self.v1 = v1_
        self.v2 = v2_
        self.label = label_

class Graph:
    def __init__(self,id_,vertices_,edges_):
        self.id = id_
        self.vertices = vertices_
        self.edges = edges_