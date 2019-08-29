
from collections import Counter

def BetterClumpFinder(genome,k,L,t):
	"""Clump finder that is efficient enough to search the entire genome of E Coli."""
    n = len(genome)-k +1
    m = L-k+1
    pattern_count = {}
    pattern_list = []
    kmers = []

    for i in range(n):
    	pattern_list.append(genome[i:i+k])
    	pattern_count[genome[i:i+k]] = 0

    for i in range(m):
    	pattern_count[genome[i:i+k]] += 1

    for i in pattern_count:
    	if pattern_count[i] == 3:
    		kmers.append(i)

    for i in range(m,n):
    	pattern_count[genome[i:i+k]] += 1
    	pattern = pattern_list[i-m]
    	pattern_count[pattern] -= 1
    	if pattern_count[genome[i:i+k]] == 3:
    		
    		kmers.append(genome[i:i+k])
    		
    return len(list(dict.fromkeys(kmers)))



with open('E_coli (1).txt') as handle:
    data = [x for x in handle.read().split()]


print BetterClumpFinder(data[0],9,500,3)