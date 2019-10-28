
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
		if spectrum.count(i) > 0:
			count += min(theo_spec.count(i), spectrum.count(i))
	return count

def ExpandNP(leaderboard):
	'''expands a list of peptide mass sequences by adding each mass between 57 and 200 to each sequence'''
	m_l = list(range(57,201))
	if leaderboard == [['']]:
		return [[aa] for aa in m_l]
	return [pep + [aa] for pep in leaderboard for aa in m_l]

def TrimNP(leaderboard,spectrum,N):
	'''trims the leaderboard from the leaderboard sequencing algorithm to keep the top N peptides including ties for the last place'''
	if len(leaderboard) == 0:
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
			c += 1
		except IndexError:
			print "warning: index error at trimming"
			break
		
	return leaderboard_trimmed

scores = dict()
def LeaderboardCyclopeptideSequencingNP(spectrum,N):
	'''Branch-and-bound algorithm which sequences a cyclopeptide given its mass spectrum, implements a leaderboard to account for error in the experimental spectrum.'''
	#we can make this MORE ACCURATE by knowing the amino acid composition of the sequence and only incorporate that alphabet instead of all amino acids!
	leaderboard = [['']]
	leader_pep = [""]
	while leaderboard != []:
		
		leaderboard = (ExpandNP(leaderboard))
		
		leaderboard_c = copy.copy(leaderboard)
		#print leaderboard
		for pep in leaderboard_c:
			pep_mass = sum(pep)

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
				leaderboard.remove(pep)
				#print "peptide exceeded parent mass: REMOVED"

		print "creating scores"
		for pep in leaderboard:
			score = PepScore(LinSpectrumNP(pep), spectrum)
			if score in scores:
				scores[score].append(pep)
			else:
				scores[score] = list([pep])
		print "done"

		leaderboard = TrimNP(leaderboard, spectrum, N)
		print leaderboard
	print "final seq list length", len(leader_pep)
	return list(set(tuple(x) for x in leader_pep))

start_time = time.time()

with open("dataset_103_2.txt") as handle:
	data = [int(x) for x in handle.read().split()]

print data[1:]
output = LeaderboardCyclopeptideSequencingNP(data[1:], 1000)
for i in output:
	print "-".join([str(j) for j in i])
#print ProteinToMass(output)


print ("--- %s seconds ---" % (time.time() - start_time))
