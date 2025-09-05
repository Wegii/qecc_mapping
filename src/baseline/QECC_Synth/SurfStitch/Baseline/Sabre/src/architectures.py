import numpy as np
import json
import random
# from qiskit.test.mock import FakeTokyo

triangle = np.array([[0,1,0], [0,0,1], [1,0,0]])
ibmqx4 = np.array([[0,0,0,0,0],[1,0,0,0,0], [1,1,0,0,0], [0,0,1,0,1],[0,0,1,0,0]])
k5 = np.array([[0,1,1,1,1], [1,0,1,1,1], [1,1,0,1,1], [1,1,1,0,1],[1,1,1,1,0]])
ring5 = np.array([[0,1,0,0,0], [0,0,1,0,0], [0,0,0,1,0], [0,0,0,0,1],[1,0,0,0,0]])
ibmqx5 = np.array([[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0], [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0], [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
                   [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0], [0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0], [0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],
                   [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0], [0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0], [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],
                   [0,0,0,0,0,1,0,0,0,0,0,1,0,1,0,0], [0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0], [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0]])

ibmToronto = np.array([[0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
              [0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
              [0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
              [0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
              [0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
              [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0],
              [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0],
              [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,1],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0],
              [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,1,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0]] )

ibmTokyo = np.array([[0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [1,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0], [0,1,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0], [0,0,1,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0], [0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0], 
                     [1,0,0,0,0,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0], [0,1,1,0,0,1,0,1,0,0,1,1,0,0,0,0,0,0,0,0], [0,1,1,0,0,0,1,0,1,0,0,0,1,1,0,0,0,0,0,0], [0,0,0,1,1,0,0,1,0,1,0,0,1,1,0,0,0,0,0,0], [0,0,0,1,1,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0], 
                     [0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,1,0,0,0,0], [0,0,0,0,0,1,1,0,0,0,1,0,1,0,0,0,1,1,0,0], [0,0,0,0,0,0,0,1,1,0,0,1,0,1,0,0,1,1,0,0], [0,0,0,0,0,0,0,1,1,0,0,0,1,0,1,0,0,0,1,1], [0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,1], 
                     [0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0], [0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,1,0,0], [0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,1,0], [0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,1], [0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0]])

def linearArch(n):
    graph = np.zeros((n,n))
    graph[0][1] = 1
    for i in range(1,n-1):
        graph[i][i+1] = 1
        graph[i][i-1] = 1
    graph[n-1][n-2] = 1
    return graph

def meshArch(n,m):
    graph = np.zeros((n*m,n*m))
    for i in range(n*m):
        for j in range(n*m):
            if neighbors(i,j,m):
                graph[i][j] = 1
    return graph

def neighbors(i,j,n):
    oneRowOver =  abs(j//n - i//n) == 1
    sameColumn = j % n  == i % n
    columnLeft = j % n  == (i % n - 1)
    columnRight =  (j%n)  == (i % n + 1)
    if (i // n) %  2 == 0:
        return oneRowOver and (sameColumn or columnLeft)
    else:   
        return oneRowOver and (sameColumn or columnRight)

def knockoutNQubits(arr, n):
    newCM = arr
    faulty = random.sample(range(len(arr)), n)
    newCM = np.delete(newCM, faulty, axis=0)
    newCM = np.delete(newCM, faulty, axis=1)
    return 

def tokyo_all_diags():
    graph = np.copy(ibmTokyo)
    for i in range(4):
        for j in range(3):
            topleft = j*5+i
            bottomRight = (j+1)*5+i+1
            graph[topleft][bottomRight] = 1
            graph[bottomRight][topleft] = 1
            topRight = j*5+i+1
            bottomLeft = (j+1)*5+i
            graph[topRight][bottomLeft] = 1
            graph[bottomLeft][topRight] = 1
    return graph

def tokyo_no_diags():
    graph = np.copy(ibmTokyo)
    for i in range(4):
        for j in range(3):
            topleft = j*5+i
            bottomRight = (j+1)*5+i+1
            graph[topleft][bottomRight] = 0
            graph[bottomRight][topleft] = 0
            topRight = j*5+i+1
            bottomLeft = (j+1)*5+i
            graph[topRight][bottomLeft] = 0
            graph[bottomLeft][topRight] = 0
    return graph


# tokyo has (0,1), (0,3), (1,0), (1,2), (2, 1), (2,3)
def tokyo_minus(squareList):
    graph = ibmTokyo
    for (j,i) in squareList:
            topleft = j*5+i
            bottomRight = (j+1)*5+i+1
            graph[topleft][bottomRight] = 0
            graph[bottomRight][topleft] = 0
            topRight = j*5+i+1
            bottomLeft = (j+1)*5+i
            graph[topRight][bottomLeft] = 0
            graph[bottomLeft][topRight] = 0
    return graph

def tokyo_plus(squareList):
    graph = ibmTokyo
    for (j,i) in squareList:
            topleft = j*5+i
            bottomRight = (j+1)*5+i+1
            graph[topleft][bottomRight] = 1
            graph[bottomRight][topleft] = 1
            topRight = j*5+i+1
            bottomLeft = (j+1)*5+i
            graph[topRight][bottomLeft] = 1
            graph[bottomLeft][topRight] = 1
    return graph

def tokyo_drop_worst_n(n, err_rates_map):
    graph = np.copy(ibmTokyo)
    worst_n = dict(sorted(err_rates_map.items(), key=lambda item: item[1], reverse=True)[:2*n]).keys()
    for (u,v) in np.argwhere(graph > 0):
        if (u,v) in worst_n:
            graph[u][v] = 0
            graph[v][u] = 0
    return graph 
                
def generateMQTFile(cm, fname):
    with open(fname, "w") as f:
        f.write(str(len(cm)) + "\n" )
        for (u, v) in np.argwhere(cm>0):
            f.write(str(u) + " " + str(v) + "\n" )


def generateEnfFile(cm, fname): 
    obj = {
        "qubits" : len(cm),
        "registers" : [{"name": "q" , "qubits": len(cm)}],
        "adj" : [ [ { "v" : "q[" + str(v) + "]"} for v in range(len(cm)) if cm[u][v]==1 ]  for u in range(len(cm)) ]
    }
    with open(fname, "w") as f:
        json.dump(obj, f)

def fake_linear_error_map():
    vals = [ 0.0120651070, 0.0120651070, 0.0219138264,0.0219138264, 0.0353320709,0.0353320709, 0.0434709196,0.0434709196, 0.0446780968, 0.0446780968]
    arch = linearArch(6)
    edges = [tuple(edge) for edge in np.argwhere(arch>0)]
    return dict(zip(edges, vals))

def fake_linear_error_list():
    return list(fake_linear_error_map().values()) 

def write_triq_files(error_map):
    with open("ibmtokyo_T.rlb", 'w') as f:
        f.write(str(len(error_map))+"\n")
        for (edge, error_rate) in error_map.items():
            f.write(str(edge[0]) + " " + str(edge[1]) + " " + str(1-error_rate) + "\n")
    with open("ibmtokyo_M.rlb", 'w') as f:
        f.write("20\n")
        for i in range(20):
            f.write(str(i) + " 1.0" + "\n")
    with open("ibmtokyo_S.rlb", 'w') as f:
        f.write("20\n")
        for i in range(20):
            f.write(str(i) + " 1.0" + "\n")
    
def square(m = 5, n = 5):
    if not m:
        m = 5
        n = 5
    
    Edges = []
    
    for i in range(m):
        for j in range(n - 1):
            Edges.append((i * n + j, i * n + j + 1))
    
    for i in range(m - 1):
        for j in range(n):
                Edges.append((i * n + j, (i + 1) * n + j))
    
    graph = np.zeros((m * n, m * n), dtype=np.int64)
    for (i, j) in Edges:
        graph[i][j] = 1
        graph[j][i] = 1
    
    index2coor = []
    for i in range(m):
        for j in range(n):
            index2coor.append((i, j))
            
    return graph

def hexagon(m = 3, n = 7):
    if not m:
        m = 3
        n = 7
         
    Edges = []
    
    for i in range(m):
        for j in range(n - 1):
            Edges.append((i * n + j, i * n + j + 1))
    
    for i in range(m - 1):
        for j in range(n):
            if (i + j) % 2 == 1:
                Edges.append((i * n + j, (i + 1) * n + j))
    
    graph = np.zeros((m * n, m * n), dtype=np.int64)
    for (i, j) in Edges:
        graph[i][j] = 1
        graph[j][i] = 1
    
    index2coor = []
    for i in range(m):
        for j in range(n):
            index2coor.append((i, j))
            
    return graph



def octagon(m = 3, n = 9):
    if not m:
        m = 3
        n = 9
    
    Edges = []
    
    for i in range(m):
        for j in range(n - 1):
                Edges.append((i * n + j, i * n + j + 1))
    
    for i in range(m - 1):
        for j in range(n):
            if (i*2 + j) % 4 == 0 or (i*2 + j) % 4 == 1 :
                Edges.append((i * n + j, (i + 1) * n + j))
    
    graph = graph = np.zeros((m * n, m * n), dtype=np.int64)
    for (i, j) in Edges:
        graph[i][j] = 1
        graph[j][i] = 1
        
    index2coor = []
    for i in range(m):
        for j in range(n):
            index2coor.append((i, j))    
    
    return graph

def square_minus():
    Edges = []
    
    for i in range(3):
        for j in range(4):
            Edges.append((i * 5 + j, i * 5 + j + 1))
    
    for i in range(2):
        for j in range(5):
                Edges.append((i * 5 + j, (i + 1) * 5 + j))
    Edges.extend([(15, 17), (16, 17), (18, 17), (16, 1), (17, 2), (18, 3)])
    graph = graph = np.zeros((19, 19), dtype=np.int64)
    for (i, j) in Edges:
        graph[i][j] = 1
        graph[j][i] = 1
    return graph


def tokyo(m, n):
    if not m:
        m = 5
        n = 5
    
    def coor2ind(i, j):
        return i * n + j
    Edges = []
    
    for i in range(m):
        for j in range(n - 1):
            Edges.append((coor2ind(i, j), coor2ind(i, j+1)))
    
    for i in range(m - 1):
        for j in range(n):
            Edges.append((coor2ind(i, j), coor2ind(i+1, j)))
            
    for i in range(m - 1):
        for j in range(n - 1):
            if (i + j) % 2 == 0:
                Edges.append((coor2ind(i, j), coor2ind(i+1, j+1)))
                Edges.append((coor2ind(i+1, j), coor2ind(i, j+1))) 
                    
    graph = graph = np.zeros((m * n, m * n), dtype=np.int64)
    for (i, j) in Edges:
        graph[i][j] = 1
        graph[j][i] = 1
    return graph, []


def tokyo_plus(m, n):
    if not m:
        m = 5
        n = 5
    
    def coor2ind(i, j):
        return i * n + j
    Edges = []
    
    for i in range(m):
        for j in range(n - 1):
            Edges.append((coor2ind(i, j), coor2ind(i, j+1)))
    
    for i in range(m - 1):
        for j in range(n):
            Edges.append((coor2ind(i, j), coor2ind(i+1, j)))
            
    for i in range(m - 1):
        for j in range(n - 1):
            Edges.append((coor2ind(i, j), coor2ind(i+1, j+1)))
            Edges.append((coor2ind(i+1, j), coor2ind(i, j+1))) 
                    
    graph = graph = np.zeros((m * n, m * n), dtype=np.int64)
    for (i, j) in Edges:
        graph[i][j] = 1
        graph[j][i] = 1
    return graph


def heavy_hexagon(m, n):
    Edges = []
    
    for i in range(m):
        for j in range(n - 1):
            Edges.append((i * n + j, i * n + j + 1))
    
    for i in range(m - 1):
        for j in range(n):
            if (i + j) % 2 == 1:
                Edges.append((i * n + j, (i + 1) * n + j))
                
    H_Edges = []
    k = 0
    for (i, j) in Edges:
        H_Edges.append((i, m * n + k))
        H_Edges.append((j, m * n + k))
        k += 1
    
    graph = np.zeros((m * n + k, m * n + k), dtype=np.int64)
    for (i, j) in H_Edges:
        graph[i][j] = 1
        graph[j][i] = 1
        
    idx2coord = []
    for i in range(m):
        for j in range(n):
            idx2coord.append((2 * i + 1, 2 * j + 1))
    for (i, j) in Edges:
        x = (idx2coord[i][0] + idx2coord[j][0]) // 2
        y = (idx2coord[i][1] + idx2coord[j][1]) // 2
        idx2coord.append((x, y))
    
    return graph

def heavy_square(m = 5, n = 5):
    Edges = []
    
    for i in range(m):
        for j in range(n - 1):
            Edges.append((i * n + j, i * n + j + 1))
    
    for i in range(m - 1):
        for j in range(n):
                Edges.append((i * n + j, (i + 1) * n + j))
    
    H_Edges = []
    k = 0
    for (i, j) in Edges:
        H_Edges.append((i, m * n + k))
        H_Edges.append((j, m * n + k))
        k += 1
    
    graph = np.zeros((m * n + k, m * n + k), dtype=np.int64)
    for (i, j) in H_Edges:
        graph[i][j] = 1
        graph[j][i] = 1
        
    index2coor = []
    for i in range(m):
        for j in range(n):
            index2coor.append((2 * i + 1, 2 * j + 1))
    for (i, j) in Edges:
        x = (index2coor[i][0] + index2coor[j][0]) // 2
        y = (index2coor[i][1] + index2coor[j][1]) // 2
        index2coor.append((x, y))
    return graph
