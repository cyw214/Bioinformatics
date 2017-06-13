#!/usr/bin/python
import sys
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
    M_nodes_height[node] = 0.0
namelist = {}
for x in range(len(M_nodes)):
	namelist[M_nodes[x]] = 1
it = 0
Liu = {}
N = len(namelist)
Liu[it] = []
while N > 2:
	it = it + 1
	for row in range(1,N+1):
		Liu[it-1].append((sum(M[row+i][row-1] for i in range(N-row)) + sum(M[row-1][j] for j in range(row-1)))/(N-2))
	mindist = float("inf")
	for i in range(len(M)):
	    for j in range(len(M[i])):
	        if (M[i][j]-Liu.get(it-1)[i]-Liu.get(it-1)[j]) < mindist:
	            mindist = M[i][j]-Liu.get(it-1)[i]-Liu.get(it-1)[j]
	            dij = M[i][j]
	            min_i, min_j = i, j
	print M_nodes[min_i],
	print M_nodes[min_j]
	LNode1ui = M[min_i][min_j]/2.0 + (Liu.get(it-1)[min_i] - Liu.get(it-1)[min_j])/2.0
	LNode2ui = M[min_i][min_j]/2.0 + (Liu.get(it-1)[min_j] - Liu.get(it-1)[min_i])/2.0
	newnode =  "U" + str(it)

	M_nodes_height[newnode + " , " + M_nodes[min_i]] = LNode1ui
	M_nodes_height[newnode + " , " + M_nodes[min_j]] = LNode2ui

	if min_j < min_i:
		min_i, min_j = min_j, min_i
	Liu[it] = []
	M_nodes[min_i] = newnode
	M_nodes.pop(min_j)
	mergerow = []
	[mergerow.append((M[min_i][i] + M[min_j][i]-dij)/2.0) for i in range(0, min_i)]      
	M[min_i] = mergerow
	[(M[i][min_i]+M[min_j][i]-dij)/2.0 for i in range(min_i+1, min_j)]
	for i in range(min_j+1, len(M)):
		M[i][min_i] = (M[i][min_i]+M[i][min_j]-dij)/2.0
		M[i].pop(min_j)
	del M[min_j]
	N = len(M_nodes)
M_nodes_height[M_nodes[0] +" , " + M_nodes[1]] = M[1][0]/2
print M_nodes
print M 
print("%10s "%'Between'),
print("%10s "%'And'),
print('Length')
for key in M_nodes_height:
	value = key.split(',')
	if(len(value) > 1):
		print("%10s "%value[0]),
		print("%10s "%value[1]),
		print str(M_nodes_height[key])

	




