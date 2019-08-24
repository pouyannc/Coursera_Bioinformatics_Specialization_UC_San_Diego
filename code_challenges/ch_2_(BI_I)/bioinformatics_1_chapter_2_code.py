import math
import random
import numpy as np

#Brute Force Algorithm:
def MotifEnumeration(dna, k, d):
	'''returns all (k,d)-motifs (k-mers with d base changes) which are present in entire list of dna'''
	#Slow algorithm
	patterns = []
	l = len(dna)
	n = len(dna[0]) - k + 1
	neighbor_l = []
	for i in range(n):
		neighbor_l.append(Neighbors(dna[0][i:i+k],d))
	neighbor_l = [item for sublist in neighbor_l for item in sublist]
	for h in neighbor_l:
		count = 0
		for i in range(1,l):
			for j in range(n):
				if HammingDistance(h, dna[i][j:j+k]) <= d:
					count += 1
					if count == l-1:
						patterns.append(h)
					break
	return " ".join(list(dict.fromkeys(patterns)))

def HammingDistance(s, t):
    '''returns hamming distance between two strings'''
    dH = len(s)
    for i in range(dH):
        if s[i] == t[i]:
            dH -= 1
    return dH

def MinHD(dna, pattern):
	k = len(pattern)
	n = len(dna) - k +1
	d = HammingDistance(dna[0:k], pattern)
	for i in range(1,n):
		if HammingDistance(dna[i:i+k], pattern) < d:
			d = HammingDistance(dna[i:i+k], pattern)
	return d


def Suffix(string):
	'''returns the string without the first letter'''
	return string[1:]

def FirstSymbol(string):
	return string[0]


def Neighbors(pattern,d):
	'''return all strings that are d nucleotides different from pattern'''
	if d == 0:
		return pattern
	if len(pattern) == 1:
		return ['A','C','G','T']
	neighborhood = []
	suffix_neighbors = Neighbors(Suffix(pattern),d)
	for text in suffix_neighbors:
		if HammingDistance(Suffix(pattern), text) < d:
			for x in ['A','C','G','T']:
				neighborhood.append(x+text)
		else:
			neighborhood.append(FirstSymbol(pattern)+text)
	return neighborhood

def Entropy(p):
	e = 0
	for i in p:
		if i != 0:
		
			e = e + (i*(math.log(i,2)))
	return e*(-1)

def Count(motifs):
	'''returns a count matrix from a motif list'''
	l = len(motifs)
	n = len(motifs[0])
	count_matrix = {"A":[], "C":[], "G":[], "T":[]}
	for key in count_matrix:
		for i in range(n):
			count_matrix[key].append(float(1)) #****implementing pseudocounts
	for i in range(n):
		for j in range(l):
			count_matrix[motifs[j][i]][i] += 1
	return count_matrix

def Profile(motifs):
	'''returns a profile matrix from a motif list'''
	profile_matrix = Count(motifs)
	for i in ['A','C','G','T']:
		for j in range(len(motifs[0])):
			profile_matrix[i][j] /= (len(motifs) + 4) #****implementing pseudocounts
	return profile_matrix

def NumberToPattern(index, k):
    '''returns the k-mer at index from all possible k-mers list (0-4^k long)'''
    base_values = {0: 'A',1:'C',2:'G',3:'T'}
    if k == 1:
        return base_values[index]
    prefix_index = index/4
    r = index%4
    symbol = base_values[r]
    prefix_pattern = NumberToPattern(prefix_index, k-1)
    return (prefix_pattern+symbol)

def MedianString(k, dna):
	'''Brute Force algorithm: returns the pattern which minimizes haming distance among all strings in dna'''
	all_kmers = []
	median = ''
	for i in range(4**k):
		all_kmers.append(NumberToPattern(i, k))

	distance = float('inf')
	for i in all_kmers:
		if distance > sum([MinHD(x, i) for x in dna]):
			distance = sum([MinHD(x, i) for x in dna])
			median = i
	return median, distance

def Pr(pattern,profile):
	count = 1
	for i in range(len(pattern)):
		count = count*(profile[pattern[i]][i])
	return count

def ProfileMostProbableKmer(text,k,profile):
	n = len(text) - k + 1
	p = -1
	best_kmer = ''
	for i in range(n):
		if p < Pr(text[i:i+k],profile):
			p = Pr(text[i:i+k],profile)
			best_kmer = text[i:i+k]
	return best_kmer

def Score(motifs):
	n = len(motifs[0])
	l = len(motifs)
	profile = Profile(motifs)
	score = 0
	score_list = []
	for i in range(n):
		score_list.append([profile[x][i] for x in ['A','C','G','T']])
	for i in score_list:
		score += (1-max(i))
	return score

def GreedyMotifSearch(k, t, dna):
	''' greedy algorithm: returns a set of best scoring k length motifs in dna based on hamming distance'''
	n = len(dna[0]) -k +1
	best_motifs = [dna[x][0:k] for x in range(t)]
	motif = ''
	for i in range(n):
		motifs = []
		motifs.append(dna[0][i:i+k])
		for j in range(1,t):
			profile = Profile(motifs)
			motif = ProfileMostProbableKmer(dna[j],k, profile)
			motifs.append(motif)
		if Score(motifs) < Score(best_motifs):
			best_motifs = np.copy(motifs)
	return best_motifs

def Motifs(profile, dna):
	motifs = []
	k = len(profile.values()[0])
	for i in range(len(dna)):
		motifs.append(ProfileMostProbableKmer(dna[i], k ,profile))
	return motifs

def Random(p_list):
	'''returns a random number from list of probabilities while accounting for weights = probabilities '''
	sum_p_list = sum(p_list)
	p_list_w = np.copy(p_list)
	for i in range(len(p_list)):
		p_list_w[i] = (p_list_w[i]/sum_p_list)
	choice = np.random.choice(p_list,1,p = p_list_w)
	#print list(p_list).index(choice[0])
	return p_list.index(choice[0])

def RandomizedMotifSearch(k,t,dna):
	n = len(dna[0]) - k +1
	final_motifs = []
	for i in range(t):
		r = random.randint(0,n-1)
		final_motifs.append(dna[i][r:r+k])
	for i in range(100):
		c = 0
		motifs = []
		for i in range(t):
			r = random.randint(0,n-1)
			motifs.append(dna[i][r:r+k])
		best_motifs = np.copy(motifs)
		while c == 0:
			profile = Profile(motifs)
			motifs = Motifs(profile,dna)
			if Score(motifs) < Score(best_motifs):
				best_motifs = np.copy(motifs)
			else:
				c = 1
		if Score(best_motifs) < Score(final_motifs):
			print "Final Motif Set Updated"
			final_motifs = np.copy(best_motifs)
	return final_motifs

def ProfileRandomGeneratedKmer(text,k,profile):
	n = len(text) - k + 1
	p_list = []
	best_kmer = ''
	for i in range(n):
		p_list.append(Pr(text[i:i+k],profile))
	kmer_i = Random(p_list)
	#print kmer_i
	best_kmer = text[kmer_i:kmer_i+k]
	return best_kmer

def GibbsSampler(k,t,N,dna):
	''' a smarter randomized algorithm'''
	c = 0
	n = len(dna[0]) - k +1
	best_motifs = []
	motifs = list(RandomizedMotifSearch(k,t,dna))
	'''for i in range(t):
					r = random.randint(0,n-1)
					motifs.append(dna[i][r:r+k])'''
	best_motifs = np.copy(motifs)
	print Score(best_motifs)
	for i in range(N):
		j = random.randint(0,t-1)
		r_motif = motifs.pop(j)
		profile = Profile(motifs)
		r_motif = ProfileRandomGeneratedKmer(dna[j],k,profile)
		motifs.insert(j,r_motif)
		if Score(motifs) < Score(best_motifs):
			print "Updated Best Motifs - Score: " +str(Score(best_motifs))
			best_motifs = np.copy(motifs)
	return best_motifs





'''motifs = ["AAATTGACGCAT",
"GACGACCACGTT",
"CGTCAGCGCCTG",
"GCTGAGCACCGG",
"AGTTCGGGACAG"]

profile = {

    'A': [],
    'C': [],
    'G': [],
    'T': []
}'''





'''l1 = []
matrix = (Profile(motifs)).values()

for i in range(len(matrix[0])):
	l1.append([x[i] for x in matrix])

total = 0
for x in l1:
	total += Entropy(x)'''



with open('dataset_163_4.txt') as handle:
	data = [x for x in handle.read().split()]

'''profile = {

    'A': [],
    'C': [],
    'G': [],
    'T': []
}

n_list = ['A','C','G','T']

q=0
b=0
for i in data[2:]:
	profile[n_list[b]].append(float(i))
	q+=1
	if q%int(data[1] )== 0:
		b+=1

print profile'''
#print data[1]


#print Random([0.1,0.2,0.7])

answer =  GibbsSampler(int(data[0]),int(data[1]),int(data[2]),data[3:])
for i in answer:
	print i

#print MedianString(7, ["CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC","GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC","GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG"])
#print (MedianString(int(data[0]), data[1:] ))
