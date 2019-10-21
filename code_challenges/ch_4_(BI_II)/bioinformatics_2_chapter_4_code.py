
import sys
import numpy as np
import collections

from collections import Counter
import copy



class Peptide(object):

	#cache = TTLCache(maxsize = 10000000, ttl = 300)

	def __init__ (self,string):
		'''initialize Peptide instance'''
		self.string = string
		self.mass = self.ProteinMass(self.string)
		self.cyc_spec = self.CycloSpectrum(string)

	def ProteinMass(self, string):
		'''return protein mass from sequence'''
		aa_mass = "G 57 A 71 S 87 P 97 V 99 T 101 C 103 L 113 I 113 N 114 D 115 K 128 Q 128 E 129 M 131 H 137 F 147 R 156 Y 163 W 186"
		mass_l = aa_mass.split()
		mass_d = dict(zip(mass_l[0::2], [int(x) for x in (mass_l[1::2])]))
		mass = 0
		mass = sum([mass_d[x] for x in string])
		return mass

	def CycloSpectrum(self,string):
		'''given a circular peptide string, return the masses of all subpeptides'''
		spectrum = [0, self.mass]
		l = len(string)
		for i in xrange(l):
			for j in xrange(1,l):
				spectrum.append(self.ProteinMass((string*2)[i:i+j]))
		return (spectrum)

	def ProteinToMass(self):
		seq = ""
		aa_mass = "G 57 A 71 S 87 P 97 V 99 T 101 C 103 L 113 I 113 N 114 D 115 K 128 Q 128 E 129 M 131 H 137 F 147 R 156 Y 163 W 186"
		mass_l = aa_mass.split()
		aa_d = dict(zip(mass_l[0::2], [int(x) for x in (mass_l[1::2])]))
		for i in self.string:
			seq += str(aa_d[i])
			seq += "-"
		return seq[:-1]


def ProteinTranslation(pattern):
	'''convert an RNA sequence into a amino acid sequence'''


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

def PepMass(pep):
	aa_mass = "G 57 A 71 S 87 P 97 V 99 T 101 C 103 I 113 N 114 D 115 K 128 E 129 M 131 H 137 F 147 R 156 Y 163 W 186"
	mass_l = aa_mass.split()
	mass_d = dict(zip(mass_l[0::2], [int(x) for x in (mass_l[1::2])]))
	return mass_d[pep]

def ProteinToMass(prot):
	seq = ""
	aa_mass = "G 57 A 71 S 87 P 97 V 99 T 101 C 103 L 113 I 113 N 114 D 115 K 128 Q 128 E 129 M 131 H 137 F 147 R 156 Y 163 W 186"
	mass_l = aa_mass.split()
	aa_d = dict(zip(mass_l[0::2], [int(x) for x in (mass_l[1::2])]))
	for i in prot:
		seq += str(aa_d[i])
		seq += "-"
	return seq[:-1]

def DToR(DNA_string):
    n = DNA_string.replace('T','U')
    return n

def Complement(string):
    reverse = string[::-1]
    reverse = reverse.replace("A","Ax")
    reverse = reverse.replace("C","Cx")
    reverse = reverse.replace("G","C")
    reverse = reverse.replace("T","A")
    reverse = reverse.replace("Ax","T")
    reverse = reverse.replace("Cx","G")
    return reverse

def PeptideEncoding(text, peptide):
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

def LinSpectrum(pep_string):
	'''given a linear peptide string, return the masses of all subpeptides'''
	spectrum = [0]
	pep_string_extend = ''
	l = len(pep_string)
	r = l * (l+1) //2
	for i in range(r):
		pep_string_extend +=  pep_string
	c = 1
	while l > 0:
		for i in range(l):
			spectrum.append(ProteinMass(pep_string_extend[i:i+c]))
		c +=1
		l -=1
	spectrum.append(ProteinMass(pep_string))
	return (spectrum)

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
		cand_pep_c = np.copy(Expand(cand_pep))
		for pep in cand_pep_c:
			if ProteinMass(pep) == spectrum[-1]:
				if (all([x in spectrum for x in LinSpectrum(pep)]) == True) & ([PepMass(x) for x in pep] not in final_pep):
						final_pep.append([PepMass(x) for x in pep])
						cand_pep.remove(pep)
			elif all([x in spectrum for x in LinSpectrum(pep)]) != True:
				cand_pep.remove(pep)
	return final_pep

def CyclopeptideScoring(string,spectrum):
	'''returns the score of a peptide with regards to an experimental spectrum, based on matching with theoretical linspectrum '''

	#initialize count
	count  = 0
	
	#assign score of -1 to peptides with mass grater than parent mass, to be filtered later
	if string.mass > int(spectrum[-1]):
		print "FILTERED!"
		return -1

	#assign the cyclic theoretical spectrum
	theo_spec = string.cyc_spec
	#assign the score based on every substring in common with spectrum
	count = sum([min(theo_spec.count(i), spectrum.count(str(i))) for i in set(theo_spec)])

	#Similar efficiency method..:
	#count  =  len(filter(lambda X: str(X) in set(spectrum),theo_spec))

	return count

def Expand(leaderboard):
	'''takes a list of peptides and returns a new expanded list: each old peptide with each amino acid attached'''
	aa = "G A S P V T C I N D K E M H F R Y W"

	#assign list of amino acids
	aa_list =  aa.split()
	#initialize list of peptides to be returned as empty list
	expanded_pep = []
	
	#check if this is the first time expand is being called..
	if leaderboard == ['']:
		for aa in aa_list:
			expanded_pep.append(Peptide(aa))
	else:
		#iterate through peptides and add each peptide + each amino acid to new list.
		for peptide in leaderboard:
			for aa in aa_list:
				expanded_pep.append(Peptide(peptide.string+aa))

	return expanded_pep

def Trim (leaderboard,N,scores):
	'''return the top N highest scoring peptides in leaderboard, with regards to spectrum '''

	#check total number of peptides
	#if total is less than N,
	if sum([len(peptide_list) for peptide_list in scores.values()]) < N:
		#assign all top scoring ties to leaderboard to be expanded (gets rid of unnecessary starter peptides)
		leaderboard = scores[max(scores.keys())]
	else:
		#add all the top scoring peptides to leaderboard until it has N values
		while len(leaderboard) < N:
			leaderboard += scores[max(scores.keys())]
			del scores[max(scores.keys())]
	return leaderboard

def LeaderboardCyclopeptideSequencing (spectrum, N):
	'''Leaderboard based algorithm for sequencing a cyclopeptide given an experimental spectrum, able to account for some error. WARNING: takes several minutes to run for large N values (>100).'''

	#initialize leaderboard
	leaderboard = ['']
	#initialize dictionary to store scores
	scores = dict()
	#create initial list 
	cand_pep = [[CyclopeptideScoring(pep, spectrum), pep] for pep in Expand(leaderboard)]

	while cand_pep != []:

		#initialize scores as an empty dictionary
		scores= dict()
		#initialize leaderboard as empty list
		leaderboard = []

		print "candidates:" ,cand_pep

		#iterate through peptides and store scores
		for i in cand_pep:
			if i[0] in scores:
				scores[i[0]].append(i[1])
			else:
				scores[i[0]] = [i[1]]

		#trim the list of peptides based on not N scores
		leaderboard = Trim(leaderboard, N, scores)
		#create a new list of candidate peptides based on expanding leaderboard
		cand_pep = filter(lambda X: X[0] != -1, [[CyclopeptideScoring(pep, spectrum), pep] for pep in Expand(leaderboard)])
		
	return leaderboard[0]

sys.setrecursionlimit(10000)
print "Recursion limit:",sys.getrecursionlimit()
with open("dataset_102_8.txt") as handle:
	data = [x for x in handle.read().split()]

Bacillus_gene = ""
for i in data:
	Bacillus_gene += i
#print data
#print data
#print PeptideEncoding(Bacillus_gene,"VKLFPWFNQY")
#for i in TheoreticalSpectrum(data[0]):
#	print i
m_list =  [57 ,71 , 87 , 97 , 99 , 101 , 103 , 113  ,114 ,115 , 128  ,129 ,131 ,137 ,147 ,156,163, 186]
l = len(m_list)
#print data[0]
#print data[1:]
print LeaderboardCyclopeptideSequencing(data[1:],int(80)).ProteinToMass()
