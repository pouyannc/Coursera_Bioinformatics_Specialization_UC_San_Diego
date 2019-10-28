import sys
#import numpy as np
import collections
import itertools
from collections import Counter
import copy
import time

def ProteinTranslation(pattern):
	'''convert an RNA sequence into the translated amino acid sequence'''
	codon_map = """UUU F      CUU L      AUU I      GUU V
	UUC F      CUC L      AUC I      GUC V
	UUA L      CUA L      AUA I      GUA V
	UUG L      CUG L      AUG M      GUG V
	UCU S      CCU P      ACU T      GCU A
	UCC S      CCC P      ACC T      GCC A
	UCA S      CCA P      ACA T      GCA A
	UCG S      CCG P      ACG T      GCG A
	UAU Y      CAU H      AAU N      GAU D
	UAC Y      CAC H      AAC N      GAC D
	UAA Stop   CAA Q      AAA K      GAA E
	UAG Stop   CAG Q      AAG K      GAG E
	UGU C      CGU R      AGU S      GGU G
	UGC C      CGC R      AGC S      GGC G
	UGA Stop   CGA R      AGA R      GGA G
	UGG W      CGG R      AGG R      GGG G"""

	rna = pattern
	aa = ''

	trans_list =  codon_map.split() #Turn string into list
	trans_dict = dict(zip(trans_list[0::2], trans_list[1::2])) #Zip together items from list to create dictionary
	for i in range(0,len(rna)-3, 3): #Iterate through dictionary using input string
		aa += trans_dict[rna[i:i+3]]

	return aa

def AAMass(aa):
	'''returns the mass of a single amino acid'''
	aa_mass = "G 57 A 71 S 87 P 97 V 99 T 101 C 103 I 113 L 113 N 114 D 115 K 128 Q 128 E 129 M 131 H 137 F 147 R 156 Y 163 W 186"
	mass_l = aa_mass.split()
	mass_d = dict(zip(mass_l[0::2], [int(x) for x in (mass_l[1::2])]))
	return mass_d[aa]

def DToR(DNA_string):
	'''returns a DNA sequence with T replaced with U (RNA)'''
	n = DNA_string.replace('T','U')
   	return n

def Complement(string):
	'''returns the reverse complement of a DNA sequence'''
	reverse = string[::-1]
	reverse = reverse.replace("A","Ax")
	reverse = reverse.replace("C","Cx")
	reverse = reverse.replace("G","C")
	reverse = reverse.replace("T","A")
	reverse = reverse.replace("Ax","T")
	reverse = reverse.replace("Cx","G")
	return reverse

def PeptideEncoding(text, peptide):
	'''return all substrings of DNA string "text" that code for "peptide"'''
	c = 0
	l = len(peptide)
	peptide_l = []
	peptide_c = 0
	rna = DToR(text)
	rc_rna = DToR(Complement(text))
	while c < 3:
		for i in range(c,len(text),3):									
			if ProteinTranslation(rna[i:i+3]+"UAA") == peptide[0]:
				if ProteinTranslation(rna[i:i+(3*l)]+"UAA") == peptide:
					#peptide_l.append(text[i:i+(3*l)])
					peptide_c+=1
			if ProteinTranslation(rc_rna[i:i+3]+"UAA") == peptide[0]:
				if ProteinTranslation(rc_rna[i:i+(3*l)]+"UAA") == peptide:
					#peptide_l.append(text[::-1][i:i+(3*l)][::-1])
					peptide_c+=1
		c+=1
	return peptide_c

def NumSubpeptides(n):
	'''given length of cyclic polypeptide, return number of possible subpeptides'''
	return n*(n-1)

def NumSubpeptidesLin(n):
	'''given length of linear polypeptide, return number of possible subpeptides'''
	s = 0
	for i in range(n):
		s += n-i
	return s+1

def CountPeptides(mass_l,mass):
	'''return number of possible sub-peptides, given mass of peptide '''
	if mass in cache:
		return cache[mass]

	if mass == 0:
		return 1
	if mass < 0:
		return  0
	value = sum([CountPeptides(mass_l,mass-x) for x in mass_l])
	cache[mass] = value
	return value

def CyclopeptideSequencing(spectrum):
	'''Branch-and-bound algorithm which sequences a cyclopeptide given its mass spectrum, does not account for errors in experimental spectrum '''
	cand_pep = ['']
	final_pep = []
	while cand_pep != []:
		try:
			cand_pep.remove('')
		except ValueError:
			pass
		cand_pep_c = copy.copy(Expand(cand_pep))
		for pep in cand_pep_c:
			if ProteinMass(pep) == spectrum[-1]:
				if (all([x in spectrum for x in LinSpectrum(pep)]) == True) & ([PepMass(x) for x in pep] not in final_pep):
						final_pep.append([PepMass(x) for x in pep])
						cand_pep.remove(pep)
			elif all([x in spectrum for x in LinSpectrum(pep)]) != True:
				cand_pep.remove(pep)
	return final_pep

def ProteinMass(string):
	'''return protein mass from sequence'''
	aa_mass = "G 57 A 71 S 87 P 97 V 99 T 101 C 103 L 113 I 113 N 114 D 115 K 128 Q 128 E 129 M 131 H 137 F 147 R 156 Y 163 W 186"
	mass_l = aa_mass.split()
	mass_d = dict(zip(mass_l[0::2], [int(x) for x in (mass_l[1::2])]))
	mass = 0
	mass = sum([mass_d[x] for x in string])
	return mass

def LinSpectrum(string):
	'''given a linear peptide string, return the masses of all subpeptides'''
	spectrum = [0, ProteinMass(string)]
	l = len(string)
	for i in xrange(l):
		for j in xrange(1,l):
			if i+j <= l:
				spectrum.append(ProteinMass(string[i:i+j]))
	return spectrum

def CycloSpectrum(string):
	'''given a circular peptide string, return the masses of all subpeptides'''
	spectrum = [0, ProteinMass(string)]
	l = len(string)
	for i in xrange(l):
		for j in xrange(1,l):
			spectrum.append(ProteinMass((string*2)[i:i+j]))
	return spectrum

def PepScore(theo_spec,spectrum):
	'''return the score of a theoretical spectrum based on values in common with experimental spectrum'''
	count = 0
	for i in set(theo_spec):
		if spectrum.count(i) > 0:
			count += min(theo_spec.count(i), spectrum.count(i))
	return count

def Expand(leaderboard):
	'''expands a list of peptides by adding each amino acid to the end of each peptide in the list, multiplying the list length by 18..'''
	aa_l = ['G','A', 'S', 'P', 'V', 'T', 'C', 'I', 'N', 'D', 'K', 'E', 'M', 'H', 'F' ,'R', 'Y', 'W']
	return [pep+aa for pep in leaderboard for aa in aa_l]

def Trim(leaderboard,spectrum,N):
	'''trims the leaderboard from the leaderboard sequencing algorithm to keep the top N peptides including ties for the last place'''
	if len(leaderboard) == 0:
		return leaderboard

	#sort the scores in descending order
	sorted_peps = [scores[score] for score in sorted(scores.keys(), reverse=True)]

	leaderboard_trimmed = set()
	c = 0

	#add top scoring peptides to new set as long as the length of the set does not exceed N
	while len(leaderboard_trimmed) < N:
		if (len(leaderboard_trimmed)) > len(leaderboard):
				print "length leader trimmed exceeded leaderboard"
				break
		try:
			for pep in sorted_peps[c]:
				leaderboard_trimmed.add(pep)
			c += 1
		except IndexError:
			print "warning: index error at trimming"
			break
		
	return leaderboard_trimmed

def ProteinToMass(leader_peps):
		'''outputs a protein sequence into its mass sequence'''
		seq_list = list()
		for pep in leader_peps:
			seq_list.append("-".join([str(ProteinMass(aa)) for aa in pep]))
		return " ".join(seq_list)

scores = dict()
def LeaderboardCyclopeptideSequencing(spectrum,N):
	'''Branch-and-bound algorithm which sequences a cyclopeptide given its mass spectrum, implements a leaderboard to account for error in the experimental spectrum.'''
	#we can make this MORE ACCURATE by knowing the amino acid composition of the sequence and only incorporate that alphabet instead of all amino acids!
	leaderboard = set([""])
	leader_pep = [""]
	while leaderboard != set([]):
		
		leaderboard = set(Expand(leaderboard))
		
		leaderboard_c = copy.copy(leaderboard)

		for pep in leaderboard_c:
			pep_mass = ProteinMass(pep)

			if pep_mass == spectrum[-1]:
				peptide_score = PepScore(CycloSpectrum(pep),spectrum)
				leader_score = PepScore(CycloSpectrum(leader_pep[0]), spectrum)
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
			score = PepScore(LinSpectrum(pep), spectrum)
			if score in scores:
				scores[score].add(pep)
			else:
				scores[score] = set([pep])
		print "done"

		leaderboard = Trim(leaderboard, spectrum, N)

	print "final seq list length", len(list(dict.fromkeys(leader_pep)))
	return list(dict.fromkeys(leader_pep))


start_time = time.time()

sys.setrecursionlimit(10000)
print "Recursion limit:",sys.getrecursionlimit()

with open("dataset_103_2.txt") as handle:
	data = [int(x) for x in handle.read().split()]

print data[1:]
output = LeaderboardCyclopeptideSequencingNP(data[1:], 1000)
for i in output:
	print "-".join([str(j) for j in i])
#print ProteinToMass(output)


print ("--- %s seconds ---" % (time.time() - start_time))
