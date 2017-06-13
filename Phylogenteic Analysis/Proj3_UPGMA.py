#!/usr/bin/python
import sys
def UPGMA(M, nodes, height_dic):
    while len(nodes) > 1:
        mindist = float("inf")
        for i in range(len(M)):
            for j in range(len(M[i])):
                if M[i][j] < mindist and not not M[i]:
                    mindist = M[i][j]
                    min_i, min_j = i, j
        height = M[min_i][min_j]
        if min_j < min_i:
            min_i, min_j = min_j, min_i
        nodes[min_i] = "(" + nodes[min_i] + ":" + str(height/2.0 - height_dic.pop(nodes[min_i])) + "," + nodes[min_j] + ":" + str(height/2.0-height_dic.pop(nodes[min_j])) + ")"
        height_dic[nodes[min_i]] = height/2.0
        nodes.pop(min_j)
        mergerow = []
        [mergerow.append((M[min_i][i] + M[min_j][i])/2.0) for i in range(0, min_i)]      
        M[min_i] = mergerow
        [(M[i][min_i]+M[min_j][i])/2.0 for i in range(min_i+1, min_j)]
        for i in range(min_j+1, len(M)):
            M[i][min_i] = (M[i][min_i]+M[i][min_j])/2.0
            M[i].pop(min_j)
        del M[min_j]
    return nodes[0]


inputlist = []
M_nodes = []
M = []
filename = sys.argv[-1]
with open(filename) as fn:
    next(fn)
    for line in fn:
        inputlist.append(line.strip())
[M_nodes.append(inputlist[i][0:10]) for i in range(len(inputlist))]
[M.append(list(row[10:len(inputlist[i])].split())) for row in inputlist]
M = [[float(M[i][j]) for j in range(len(M[i]))] for i in range(len(M))]

M_nodes_height = {}
for node in M_nodes:
    M_nodes_height[node] = 0

print UPGMA(M, M_nodes, M_nodes_height)
