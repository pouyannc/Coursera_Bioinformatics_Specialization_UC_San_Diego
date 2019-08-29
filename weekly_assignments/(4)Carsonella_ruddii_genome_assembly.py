import sys
from collections import defaultdict

def Suffix(string):
	'''returns the string without the first letter'''
	return string[1:]

def Prefix(string):
	return string[0:-1]

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


sys.setrecursionlimit(7000)
print "Recursion limit:",sys.getrecursionlimit()

#Processing data from text file into correct formatting
with open("Carsonella_ruddii_genome_read_pairs.txt") as handle:
	data = [x for x in handle.read().split()]
read_pairs = []
for x in range(2, len(data)):
	s = (data[x].replace("|" ," " )).split()
	read_pairs.append(s)  #everything is strings

#print read_pairs
genome = StringReconstructionReadPairs(120,1000,read_pairs)

