import random
import math
import numpy as np

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

with open('DosR.txt') as handle:
	data = [x for x in handle.read().split()]

m =  GibbsSampler(12,len(data),1000,data)
for i in m:
	print i


