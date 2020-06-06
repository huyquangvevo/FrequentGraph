import numpy as np

def string2matrix(st):
    strMatrix = st[2:-2]
    rows = strMatrix.split("]\n [")
    # print(row)
    matrix = []
    for row in rows:
        rowClean = row.replace("\n","")
        matrix.append(np.fromstring(rowClean,dtype=int,sep=' '))
    # print(np.array(matrix))
    return np.array(matrix).copy()

def encodeGraph(graph):
    visited = [False]*len(graph)
    queue = []
    queue.append(0)
    visited[0] = True
    code = str(graph[0,0]) + '$'
    while queue:
        s = queue.pop(0)
        levelStr = ''
        for i in np.where(graph[s]>0)[0][s+1:]:
            queue.append(i)
            levelStr += str(graph[s,i]) + "_" + str(graph[i,i]) + "_"
            visited[i] = True 
        if levelStr != '':
            code += levelStr[:-1] +  '$'
    code += '#'

    return code

def canonicalForm(graph: np.ndarray,embeddings=None):
    labelNodes = graph.diagonal()
    start = np.zeros((1,1),dtype=int)
    maxNodes = np.where(labelNodes == np.max(labelNodes))[0]
    start[0,0] = np.max(labelNodes)
    canonical = {
        "code" : ''
    }
    for idStart in maxNodes:
        S = {
            "tree" : start,
            "index" : np.array([idStart]),
            "code" : encodeGraph(start)
        }

        while (len(S["index"]) < len(labelNodes)):
                # trees = []
            newCandidates = {}
            for i in range(graph.shape[0]):
                if i in S["index"]:
                    continue
                Q = []
                t = S["tree"]
                for id,j in enumerate(S["index"]):
                    if graph[i,j] == 0:
                        continue
                    rowExpand = np.zeros((1,t.shape[0]),dtype=int)
                    rowExpand[0,id] = graph[i,j]
                    tree = np.r_[t,rowExpand]
                    colExpand = np.zeros((tree.shape[0],1),dtype=int)
                    colExpand[id,0] = graph[i,j]
                    colExpand[tree.shape[0]-1,0] = graph[i,i]
                    tree = np.c_[tree,colExpand]
                    indexTree = np.concatenate([S["index"],np.array([i])])
                    codeTree = encodeGraph(tree)
                    newCandidates[codeTree] = {
                        "tree" : tree,
                        "index" : indexTree,
                        "code" : codeTree
                    }

            S = newCandidates[max(newCandidates.keys())]
        canonical = S if canonical["code"] < S["code"] else canonical
        if embeddings is not None:
            for k in embeddings.keys():
                topo = []
                for subNodes in embeddings[k]:
                    reindexedNodes = np.array([subNodes[idNode] for idNode in canonical['index']])
                    topo.append(reindexedNodes)
                embeddings[k] = topo
    return canonical
    #return canonical["code"] if code else canonical["tree"]

def extend(X: np.ndarray,pad: np.ndarray):
    n = X.shape[0]
    X = np.r_[X,pad[:,:-1]]
    pad = np.reshape(pad,(pad[0].shape[0],1))
    X = np.c_[X,pad]
    # print("X extend",X)
    return X

def extendOneNode(X: np.ndarray, Y: np.ndarray):
    # print(X,Y)
    n = X.shape[0]
    xLabelNodes = X.diagonal()#np.unique(X.diagonal())
    extensions = []
    for i,lNode in enumerate(xLabelNodes):
        # print("node",lNode,"i",i)
        indices = np.where(Y.diagonal() == lNode)[0]
        for iY in indices:
            for j in np.where(Y[iY] > 0)[0]:
                if j != iY:
                    pad = np.zeros((1,n+1),dtype=int)
                    pad[0,-1] = Y[j,j]
                    pad[0,i] = Y[iY,j]
                    extensions.append(extend(X.copy(),pad))
    return extensions

def extendByCore(X: np.ndarray, Y: np.ndarray):
    # print(X,Y)
    n = X.shape[0]
    xLabelNodes = X.diagonal()#np.unique(X.diagonal())
    extensions = []
    for i,lNode in enumerate(xLabelNodes):
        # print("node",lNode,"i",i)
        if i != n -1:
            continue
        indices = np.where(Y.diagonal() == lNode)[0]
        for iY in indices:
            for j in np.where(Y[iY] > 0)[0]:
                if j != iY:
                    pad = np.zeros((1,n+1),dtype=int)
                    pad[0,-1] = Y[j,j]
                    pad[0,i] = Y[iY,j]
                    extensions.append(extend(X.copy(),pad))
    return extensions
