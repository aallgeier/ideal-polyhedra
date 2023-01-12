import numpy as np

"""
This function returns the adjacency matrix of the vertices of an octahedron.
"""
def createMatrix():
    #returns the adjacency matrix of the vertices of an octahedronn
    return  [[0, 1, 1, 1, 0, 1],[1, 0, 1, 1, 1, 0],[1, 1, 0, 0, 1, 1],
            [1, 1, 0, 0, 1, 1],[0, 1, 1, 1, 0, 1],[1, 0, 1, 1, 1, 0]]
#w is the number of vertices that bounds a single face
w = int(input("Number of vertices bounding a face? "))
#v is the number of all the vertices on the polyhedra
v = len((createMatrix()[0]))

"""
This function returns the edges of a graph.
Each edge is represented as a set the contains two vertices.
parameters:
    n - number of vertices on the polyhedra
"""
def getEdges(n):
    adjacency = createMatrix()
    edges = []
    for i in range(len(adjacency)):
        for j in range(i+1, n):
            if adjacency[i][j] == 1:
                edges.append(set([i, j]))
    return edges

print("-Number of edges:")
print(len(getEdges(v)))
print("-Edges:")
print(getEdges(v))

"""
This function returns the list of vertices that are adjacent to the chosen vertex.
parameters:
    n - number of all the vertices on the polyhedra
    vertex - vertex of our choice
return:
    list of vertices
"""
def adjacentVertices(n, vertex):
    adjacency = createMatrix()
    adjacentVertices = []
    for i in range(n):
        if adjacency[vertex][i] == 1:
            adjacentVertices.append(i)
    return adjacentVertices


"""
This function returns all the simple cycles that begin and end with a vertex of our choice
while avoiding specified vertices.
parameters:
    n - number of vertices
    startingVertex - vertex to start our cycle
    avoid - list of vertices we want to avoid (an empty list will
    give all the cycles that start and end at startingVertex)
return:
    list of lists of vertices that make up simple cycles
"""
def fixedVertexCycles(n, startingVertex, avoid):
    #completecycles records cycles that start and end at startingVertex that don’t contain verticies in avoid.
    completecycles = []
    #path starts from startingVertex
    #go to vertices adjacent to startingVertex that’s not in avoid
    adjacent = adjacentVertices(n, startingVertex)
    #We want to make the first step outside of the loop because we need
    #a list of lists. We don’t want a list of integers.
    afterFirstStep = []
    for i in adjacent:
        if i not in avoid:
            afterFirstStep.append([startingVertex, i])
    #We want to make the first step outside of the loop in order to avoid
    #a "cycle" with 3 vertices that go 0 -> adjacent vertex -> 0 (a digon like situation).
    afterSecondStep = []
    for i in range(len(afterFirstStep)):
        endpt = afterFirstStep[i][len(afterFirstStep[i])-1]
        adjacent = adjacentVertices(n, endpt)
        for j in adjacent:
            if j != startingVertex and j not in avoid:
                afterSecondStep.append(afterFirstStep[i]+[j])
    #Change variable name
    previouspaths = afterSecondStep
    #We choose n for the range below because we already made 2 steps and
    #a walk cannot take steps longer than the number of all vertices.
    for i in range(n):
        #latestpaths will record the updates of the paths
        latestpaths = []
        for j in range(len(previouspaths)):
            #We will check that adjacent vertices of the endpoints of the paths
            endpoint = previouspaths[j][len(previouspaths[j])-1]
            adjacent = adjacentVertices(n, endpoint)
            for k in adjacent:
                if k not in avoid:
                    if k not in previouspaths[j]:
                        latestpaths.append(previouspaths[j]+[k])
                    elif k == startingVertex:
                        completecycles.append(previouspaths[j]+[k])
        previouspaths = latestpaths.copy()
    #at this point, completecycles contain paths that are merely the reversed version of another cycle.
    #We must remove them.
    indicies = []
    for i in range(len(completecycles)):
        reversedCycle = completecycles[i][::-1]
        for j in range(i+1,len(completecycles)):
            if completecycles[j] == reversedCycle:
                indicies.append(i)
    uniquecycles = list(np.array(completecycles)[indicies])
    #make sure all elements are lists
    uniquecycles = [list(x) for x in uniquecycles]
    return uniquecycles

"""
This function returns the list of all distinct simple cycles on the 1-skeleton of the given polyhedron.
Parameters:
    n - number of vertices in a graph
#Explained to Rainie
#explained to Franco
"""
def allSimpleCycles(n):
    allCycles = []
    avoid = []
    for startingVertex in range(n):
        allCycles = allCycles + fixedVertexCycles(n, startingVertex, avoid)
        avoid.append(startingVertex)
    return allCycles
print("-Number of simple cycles (total):")
print(len(allSimpleCycles(v)))
print("-List of simple cycles (total):")
print(allSimpleCycles(v))

"""
This function assigns variables to each face.
"""
def makeVariables(n):
    edges = getEdges(n)
    variables = []
    for i in range(len(edges)):
        variables.append("x"+str(i))
    return variables


"""
 !!! Caution !!!
 !!!This only works for the platonic solids!!!
 This function returns all distinct cycles that bound a single face.
 parameters:
    n - number of vertices
    numVertexForFace - the number of vertices needed to bound a single face
 *Note: For numVertexForFace, actually put the number of vertices around a face.
 The path will contain an additional vertex
 (hence numVertexForFace+1 in the code).
"""
def boundingFaces(n, numVertexForFace):
    allCycles = allSimpleCycles(n)
    boundingFacelist = []
    for i in range(len(allCycles)):
        if len(allCycles[i]) == numVertexForFace+1:
            boundingFacelist.append(allCycles[i])
    return boundingFacelist


#print("Number of cycles bounding a face:")
#print(len(boundingFaces(v, w)))
print("-Number of simple cycles bounding a face:")
print(len(boundingFaces(v, w)))

"""
This function returns the list of simple cycles that don’t bound a single face.
"""
def notBoundingFace(n, numVertexForFace):
    allSimpleCycleslist = allSimpleCycles(n)
    notBoundingFacelist = []
    for i in range(len(allSimpleCycleslist)):
        if len(allSimpleCycleslist[i]) > numVertexForFace+1:
            notBoundingFacelist.append(allSimpleCycleslist[i])
    return notBoundingFacelist
print("-Number of simple cycles not bounding a face:")
print(len(notBoundingFace(v, w)))
print("-Number of inequalities (including edge weight restrictions:)")
print(len(notBoundingFace(v, w))+len(getEdges(v)))


"""
This function prints equations coming from the bounding face condition.
"""
def getEquationsBoundingFace(n, numVertexForFace):
    edges = getEdges(n)
    variables = makeVariables(n)
    cycles = []
    boundingFace = boundingFaces(n, numVertexForFace)
    #[[0]*4 for x in range(4)]
    matrix = [[0]*len(edges) for x in range(len(boundingFace))]
    for i in range(len(boundingFace)):
        cycle = []
        for j in range(len(boundingFace[i]) - 1):
            edge = set([boundingFace[i][j], boundingFace[i][j + 1]])
            for k in range(len(edges)):
                if edge == edges[k]:
                    cycle.append(variables[k])
                    matrix[i][k] = 1
        cycles.append(sorted(cycle))
    for i in range(len(cycles)):
        print('+ '.join(cycles[i]) + " == 2*np.pi:")
    print("-Left side of equations to matrix:")
    print(matrix)
print("-Equations and inequalities:")
getEquationsBoundingFace(v, w)

"""
This function prints the inequalities for
simple circuits not bounding a face.
"""
def getInequalitiesNotBoundingFace(n, numVertexForFace):
    cycles = []
    edges = getEdges(n)
    variables = makeVariables(n)
    nboundingFace = notBoundingFace(n, numVertexForFace)
    for i in range(len(nboundingFace)):
        cycle = []
        for j in range(len(nboundingFace[i]) - 1):
            edge = set([nboundingFace[i][j], nboundingFace[i][j + 1]])
            for k in range(len(edges)):
                if edge == edges[k]:
                    cycle.append(variables[k])

        cycles.append(sorted(cycle))
    for i in range(len(cycles)):
        print("if " + '+ '.join(cycles[i]) + " > 2*np.pi:")
        print("    count += 1")

getInequalitiesNotBoundingFace(v, w)

"""
This function prints the inequalities satisfied
by each edge weight.
"""
def getedgeWeightRange(n):
    variables = makeVariables(n)
    for i in variables:
        if i != variables[len(variables)-1]:
            print("if " + "0 < "+str(i)+" and " + str(i) + " < np.pi:")
            print("    count += 1")
        else:
            print("if " + "0 < "+str(i)+" and " + str(i) + " < np.pi:")
            print("    count += 1")
            print("{", ', '.join(makeVariables(v)), "}")
getedgeWeightRange(v)