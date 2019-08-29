'''A path in the overlap graph that visits every node once = Hamiltonian path
A path in the de Bruijn graph that visits every edge exactly once = Eulerian Path 
Eulers Theorem: Any Eulerian graph must be balanced (all node in degrees = out degrees) AND strongly connected (all nodes can lead to each other)'''
import sys
from collections import defaultdict

def StringComposition(k,text):
	s_list = []
	n =  len(text) -k +1
	for i in range(n):
		s_list.append(text[i:i+k])
	return s_list

def FirstSymbol(string):
	return string[0]

def Suffix(string):
	'''returns the string without the first letter'''
	return string[1:]

def Prefix(string):
	return string[0:-1]

def PathToGenome(path):
	genome = ""
	l = len(path)
	for i in range(l):
		if i == l-1:
			genome += path[i]
			return genome
		genome += FirstSymbol(path[i])

def OverlapGraph(patterns):
	adj_list = {}
	print "Initializing adjacency list..."
	for i in patterns:
		adj_list[i] = []
	s_list = [Suffix(x) for x in adj_list]
	p_list = [Prefix(x) for x in patterns]
	l_p = len(patterns)
	l_a = len(adj_list)
	print "Checking for matching strands..."
	print "Adding strands to list..."
	for i in range(l_p):
		for j in range(l_a):
			if s_list[i] == p_list[j]:
				adj_list.values()[i].append(patterns[j])
	print "Done"
	return adj_list

def DeBruijns(k, text):
	'''implementing de Bruijns graph, glued nodes'''
	edge_list = []
	n = len(text) - k +1
	for i in range(n):
		edge_list.append(text[i:i+k])
	node_list = [Prefix(edge_list[0])]
	node_list2 = []
	for i in edge_list:
		node_list.append(Suffix(i))
		node_list2.append(Suffix(i))
	adj = zip(node_list,node_list2)
	adj_dict = defaultdict(list)
	for x,y in adj:
		adj_dict[x].append(y)
	return adj_dict

def DeBruijnsFromKmers(patterns):
	'''returns de Bruijns graph from kmer composition of genome, using composition graphing '''
	l = len(patterns)
	node_list = []
	node_list2 = []
	for i in patterns:
		node_list.append(Prefix(i))
		node_list2.append(Suffix(i))
	adj = zip(node_list,node_list2)
	adj_dict = defaultdict(list)
	for x,y in adj:
		adj_dict[x].append(y)
	return adj_dict

def NodeDegree(adj_list, node):
	c = 0
	for i in adj_list.keys():
		for j in adj_list[i]:
			if j == node:
				c += 1
	in_degree = c
	out_degree = len(adj_list[node])
	degree = in_degree + out_degree
	return degree, in_degree, out_degree

def FindCircuit(current_node, adj_list, circuit = [], final_node = ""):
	x=0
	print "cn",current_node
	neighbors = (adj_list[current_node])
	n = (len(neighbors))
	if n == 0:
		print "circuit finished"
		circuit.append(current_node)
		return circuit, adj_list
	else:
		x = adj_list[current_node][0]
		circuit.append(current_node)
		adj_list[current_node].remove(adj_list[current_node][0])
		return FindCircuit(tuple(x),adj_list,circuit)

def EulerianCircuit(adj_list):
	current_node = adj_list.keys()[0]
	circuit, adj_list = FindCircuit(current_node, adj_list) 
	while max([len(adj_list[x]) for x in adj_list]) > 0:
		for i in range(len(circuit)):
			if len(adj_list[circuit[i]]) > 0:
				current_node = circuit[i]
				new_circ, adj_list = FindCircuit(current_node,adj_list, [])
				del circuit[i]
				circuit[i:i] = new_circ
	return circuit

def EulerianTour(adj_list):
	final_node = None
	current_node = None
	for i in adj_list.values():
		for j in i:
			if tuple(j) not in adj_list.keys():
				final_node = tuple(j)
				adj_list[tuple(j)] = []
	print "Final node set to:", final_node
	for i in adj_list.keys():
		if (current_node != None) & (final_node != None):
			break
		degree, in_degree, out_degree = NodeDegree(adj_list, i)
		if out_degree - in_degree == 1:
			current_node = i
		if (in_degree - out_degree == 1) & (final_node == None):
			current_node = i
	circuit, adj_list = FindCircuit(tuple(current_node), adj_list,[], final_node) 
	print "circuit:",circuit
	while max([len(adj_list[x]) for x in adj_list]) > 0:
		for i in range(len(circuit)):
			if len(adj_list[tuple(circuit[i])]) > 0:
				current_node = circuit[i]
				new_circ, adj_list = FindCircuit(current_node,adj_list,[], final_node)
				del circuit[i]
				circuit[i:i] = new_circ
	return circuit

def StringReconstruction(patterns):
	dB = DeBruijnsFromKmers(patterns)
	path = EulerianTour(dB)
	text = PathToGenome(path)
	return text

def GenerateBinaryStrings(k, arr, i, s_list):
	if len(s_list) == 2**k:
		return s_list
	if i == k:
		print "String added:",(arr)
		s_list.append("".join(map(str,arr)))
		return 
	arr[i] = 0
	GenerateBinaryStrings(k,arr, i+1,s_list)
	arr[i] = 1
	GenerateBinaryStrings(k,arr,i+1,s_list)
	if len(s_list) == 2**k:
		return s_list

def BinaryStringCycle(k):
	binary_strings = GenerateBinaryStrings(k, ([None]*k), 0, [])
	dB=DeBruijnsFromKmers(binary_strings)
	cycle = EulerianCircuit(dB)
	string = PathToGenome(cycle)
	return string[:(-k+1)]

def DeBruijnsFromReadPairs(paired_kmers):
	'''returns de Bruijns graph from kmer composition of genome, using composition graphing '''
	l = len(paired_kmers)
	node_list = []
	node_list2 = []
	for i in paired_kmers:
		node_list.append([Prefix(i[0]),Prefix(i[1])])
		node_list2.append([Suffix(i[0]),Suffix(i[1])])
	adj = zip(node_list,node_list2)
	adj_dict = defaultdict(list)
	for x,y in adj:
		adj_dict[tuple(x)].append(y)
	return adj_dict

def ReadPairPathToGenome(path, k,d ):
	genome = ''
	genome+= path[0][0]
	l = len(path)
	for i in range(1,l):
		genome += path[i][0][-1]
		if len(genome) == (k+d):
			break
	genome+=path[0][1]
	for i in range(1,l):
		genome += path[i][1][-1]
	return genome

def StringReconstructionReadPairs(k,d,paired_reads):
	dB = DeBruijnsFromReadPairs(paired_reads)
	path = EulerianTour(dB)
	
	sequence = ReadPairPathToGenome(path,k, d)
	return sequence

def MaximalNonBranchingPaths(graph):
	paths = []
	for i in graph.keys():
		degree,in_degree,out_degree = NodeDegree(graph, i)
		if (in_degree == out_degree == 1) == False:
			if out_degree > 0:
				for j in graph[i]:
					path = FirstSymbol(i)
					x=j
					degree,in_degree,out_degree = NodeDegree(graph, x)
					while in_degree == out_degree == 1:
						path += FirstSymbol(x)
						x = graph[x][0]
						degree,in_degree,out_degree = NodeDegree(graph, x)
					if (in_degree == out_degree == 1) == False:
						path += x
					paths.append(path)
	return paths

def ContigGeneration(patterns):
	'''return list of contigs from input of collecion of strings'''
	contigs = []
	dB = DeBruijnsFromKmers(patterns)
	contigs = MaximalNonBranchingPaths(dB)


	'''c = 0
				for i in len(dB.keys()):
								degree,in_degree,out_degree = NodeDegree(dB, dB.keys()[i])
								if len(dB.values()[i]) == 1
			
				for i in dB.keys():
					print c
					if c > 0:
						c -=1
						continue
					contig = ""
					for j in dB[i]:
			
						degree,in_degree,out_degree = NodeDegree(dB, j)
						print in_degree,out_degree
						if in_degree == out_degree == 1:
							print "here"
							x = i
							contig += x
							x = dB[x][0]
							while len(dB[x]) == 1:
								contig += x[-1]
								c += 1
								x = dB[x][0]
							contig += x
							contigs.append(contig)'''
	return " ".join(contigs)





with open("dataset_205_5.txt") as handle:
	data = [x for x in handle.read().split()]


sys.setrecursionlimit(7000)
print "Recursion limit:",sys.getrecursionlimit()
#print ((data[2].replace("->" ,"" )).replace(","," ")).split()
'''read_pairs = []
for x in range(2, len(data)):
	s = (data[x].replace("|" ," " )).split()
	read_pairs.append(s)  #everything is strings'''


contigs= ContigGeneration(data[:])
print contigs
#print len(text)
#print NodeDegree(adj_list, str(246))

'''for i in graph.items():
	if len(i[1])>0:
		print i[0] + " -> " + str(", ".join(i[1]))'''

#data = data.strip()

#print StringComposition(int(data[0], data[1]))

#answer =  PathToGenome(data)
#print answer
'''for i in answer:
	print i'''

'''txt = '\n'.join(answer)
f = open('data.txt','w')
f.write(txt)
f.close()'''
