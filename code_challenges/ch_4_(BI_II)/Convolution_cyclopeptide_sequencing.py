from itertools import chain
from collections import Counter
import copy
import time

'''This is an alternative implementation of Leaderboard cyclopeptide sequencing assuming each mass from 57-200 is a unique amino acid, hence accounting for non-proteinogenic amino acids'''
#WARNING: takes a long time for high N value (~10 min with N=1000)

def CycloSpectrumNP(string):
	'''given a circular peptide string, return the masses of all subpeptides'''
	spectrum = [0, sum(string)]
	l = len(string)
	for i in xrange(l):
		for j in xrange(1,l):
			spectrum.append(sum((string*2)[i:i+j]))
	return spectrum

def LinSpectrumNP(string):
	'''given a linear peptide string, return the masses of all subpeptides'''
	spectrum = [0, sum(string)]
	l = len(string)
	for i in xrange(l):
		for j in xrange(1,l):
			if i+j <= l:
				spectrum.append(sum(string[i:i+j]))
	return spectrum

def PepScore(theo_spec,spectrum):
	'''return the score of a theoretical spectrum based on values in common with experimental spectrum'''
	count = 0
	for i in set(theo_spec):
		j = round(i,0)
		if spectrum.count(j) > 0:
			count += min(theo_spec.count(i), spectrum.count(j))
	return count

def ExpandNP(leaderboard, alph):
	'''expands a list of peptide mass sequences by adding each mass between 57 and 200 to each sequence'''
	if leaderboard == [['']]:
		return [[aa] for aa in alph]
	return [pep + [aa] for pep in leaderboard for aa in alph]

def TrimNP(leaderboard,spectrum,N):
	'''trims the leaderboard from the leaderboard sequencing algorithm to keep the top N peptides including ties for the last place'''
	
	if len(leaderboard) == 0:
		print "empty leaderboard, nothing to trim.."
		return leaderboard

	#sort the scores in descending order
	sorted_peps = [scores[score] for score in sorted(scores.keys(), reverse=True)]

	leaderboard_trimmed = list()
	c = 0

	#add top scoring peptides to new set as long as the length of the set does not exceed N
	while len(leaderboard_trimmed) < N:
		if (len(leaderboard_trimmed)) > len(leaderboard):
			print "length leader trimmed exceeded leaderboard"
			break
		try:
			for pep in sorted_peps[c]:
				leaderboard_trimmed.append(pep)
			c+=1
		except IndexError:
			print "warning: index error at trimming"
			break
		
	return leaderboard_trimmed

def SpectralConvolution(spectrum):
	convoluted = [(i-j) for i in spectrum for j in spectrum if (200 >= (i-j) >= 57)]
	return convoluted

scores = dict()
def LeaderboardCyclopeptideSequencingNP(spectrum,N, alph):
	'''Branch-and-bound algorithm which sequences a cyclopeptide given its mass spectrum, implements a leaderboard to account for error in the experimental spectrum.'''
	leaderboard = [['']]
	leader_pep = [""]
	while leaderboard != []:
		
		leaderboard = (ExpandNP(leaderboard, alph))
		
		leaderboard_c = copy.copy(leaderboard)

		for pep in leaderboard_c:
			pep_mass = sum(pep)
			#print 'mass',pep_mass
			#use max(spectrum) if spectrum is not sorted in ascending order
			if pep_mass == spectrum[-1]:
				peptide_score = PepScore(CycloSpectrumNP(pep),spectrum)
				leader_score = PepScore(CycloSpectrumNP(leader_pep[0]), spectrum)
				print "pep score", peptide_score
				print "leader score", leader_score
				if peptide_score == leader_score:
					leader_pep.append(pep)
					print "leader list expanded..."
				elif peptide_score > leader_score:
					leader_pep = [pep]
					print "new leader list created"
			elif pep_mass > spectrum[-1]:
				#print "removed"
				leaderboard.remove(pep)

		print "creating scores"
		for pep in leaderboard:
			score = PepScore(LinSpectrumNP(pep), spectrum)
			if score in scores:
				scores[score].append(pep)
			else:
				scores[score] = list([pep])
		print "done"

		leaderboard = TrimNP(leaderboard, spectrum, N)
		print "length after trim:",len(leaderboard)
	print "final seq list length", len(leader_pep)
	return list(set(tuple(x) for x in leader_pep))

def ConvolutionCyclopeptideSequencing(M,N,spectrum):
	'''performs convolution of the spectrum in order to reduce initial amino acid alphabet down to only the ones in the target protein, followed by leaderboard cyclopeptide sequencing'''
	c = 0
	convoluted = Counter(SpectralConvolution(spectrum))
	conv_sorted = sorted(convoluted.keys(), key = lambda x : convoluted[x], reverse= True)
	reduced_alph = []
	for i in conv_sorted:
		if len(reduced_alph) == M:
			c = 1
		if c == 1 and convoluted[i] != convoluted[reduced_alph[-1]]:
			break
		reduced_alph.append(i)
	leader_pep = LeaderboardCyclopeptideSequencingNP(spectrum, N, reduced_alph)


	return leader_pep



start_time = time.time()

with open("dataset_104_8.txt") as handle:
	data = [int(x) for x in handle.read().split()]


output = ConvolutionCyclopeptideSequencing(data[0],data[1],data[2:])
for i in output:
	print "-".join([str(int(round(j,0))) for j in i])


print ("--- %s seconds ---" % (time.time() - start_time))
