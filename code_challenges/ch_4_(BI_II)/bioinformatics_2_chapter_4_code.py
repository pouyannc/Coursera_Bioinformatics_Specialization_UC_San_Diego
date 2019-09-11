import sys
import numpy as np
try:
    from functools import lru_cache
except ImportError:
    from backports.functools_lru_cache import lru_cache



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

def ProteinMass(string):
	aa_mass = "G 57 A 71 S 87 P 97 V 99 T 101 C 103 I 113 N 114 D 115 K 128 E 129 M 131 H 137 F 147 R 156 Y 163 W 186"
	mass_l = aa_mass.split()
	mass_d = dict(zip(mass_l[0::2], [int(x) for x in (mass_l[1::2])]))
	mass = 0
	#for i in string:
	mass = sum([mass_d[x] for x in string])
	return mass

def PepMass(pep):
	aa_mass = "G 57 A 71 S 87 P 97 V 99 T 101 C 103 I 113 N 114 D 115 K 128 E 129 M 131 H 137 F 147 R 156 Y 163 W 186"
	mass_l = aa_mass.split()
	mass_d = dict(zip(mass_l[0::2], [int(x) for x in (mass_l[1::2])]))
	return mass_d[pep]

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

def CycloSpectrum(pep_string):
	'''given a circular peptide string, return the masses of all subpeptides'''
	spectrum = [0]
	pep_string_extend = ''
	l = len(pep_string)
	for i in range(l*l):
		pep_string_extend +=  pep_string
	c = 1
	while c < l:
		for i in range(l):
			spectrum.append(ProteinMass(pep_string_extend[i:i+c]))
		c +=1
	spectrum.append(ProteinMass(pep_string))
	return (spectrum)

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

cache = {}

def CountPeptides(mass_l,mass):
	'''return number of possible peptides, given mass of peptide '''
	if mass in cache:
		return cache[mass]

	if mass == 0:
		return 1
	if mass < 0:
		return  0
	value = sum([CountPeptides(mass_l,mass-x) for x in mass_l])
	cache[mass] = value
	return value

def Expand(cand_pep):
	aa = "G A S P V T C I N D K E M H F R Y W"
	aa_list =  aa.split()
	cand_pep_c = np.copy(cand_pep)
	if len(cand_pep_c) == 0:
		for i in aa_list:
			cand_pep.append(i)
	else:
		for i in range(0, len(cand_pep_c)):
			for aa in aa_list:
				cand_pep.append(cand_pep_c[i]+aa)
			cand_pep.remove(cand_pep_c[i])
	return cand_pep



def CyclopeptideSequencing(spectrum):
	'''Branch-and-bound algorithm which sequences a cyclopeptide given its mass spectrum '''
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




sys.setrecursionlimit(10000)
print "Recursion limit:",sys.getrecursionlimit()
with open("dataset_100_6.txt") as handle:
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
#print CountPeptides(m_list,int(data[0]))
print [int(x) for x in data]
for i in CyclopeptideSequencing([0 ,71, 99 ,101 ,103 ,128 ,129 ,199, 200 ,204, 227, 230, 231, 298, 303, 328 ,330 ,332 ,333]):
	print "-".join ([str(x) for x in i])
#print CountPeptides(m_list,1024)**(1.000/1024)
#print NumSubpeptidesLin(20341)
#[0, 113, 128, 186, 241, 299, 314 ,427]
#[int(x) for x in data]
#[0, 71, 97, 99, 103, 113, 113, 114, 115, 131, 137, 196, 200, 202, 208, 214, 226, 227, 228, 240, 245, 299, 311, 311, 316, 327, 337, 339, 340, 341, 358, 408, 414, 424, 429, 436 ,440 ,442, 453, 455, 471, 507, 527, 537, 539, 542, 551, 554, 556, 566, 586, 622, 638, 640, 651, 653 ,657, 664, 669 ,679, 685, 735 ,752, 753, 754, 756, 766, 777, 782, 782, 794 ,848, 853, 865, 866, 867 ,879, 885, 891, 893 ,897 ,956, 962, 978, 979, 980, 980 ,990, 994, 996, 1022, 1093]