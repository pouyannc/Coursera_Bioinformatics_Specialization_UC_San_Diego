from Bio import SeqIO

def MinimumSkewProblem(genome):
    '''return the positions of minimum value(s) in the skew diagram for genome'''
    skew = [0]
    m = 0
    m_i = []
    for i in range(len(genome)):
        if genome[i] == 'C':
            skew.append(skew[i]-1)
        if genome[i] == 'G':
            skew.append(skew[i]+1)
        if genome[i] in ['A','T']:
            skew.append(skew[i])
        if skew[i] < m:
            m = skew[i]
    for i in range(len(skew)): 
        if skew[i] == m:
            m_i.append(i)        
    return m_i

with open('Salmonella_enterica.txt') as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        data = record.seq

def Complement(string):
    '''returns the reverse complement of string'''
    reverse = string[::-1]
    reverse = reverse.replace("A","Ax")
    reverse = reverse.replace("C","Cx")
    reverse = reverse.replace("G","C")
    reverse = reverse.replace("T","A")
    reverse = reverse.replace("Ax","T")
    reverse = reverse.replace("Cx","G")
    return reverse

def HammingDistance(s, t):
    '''returns hamming distance between two strings'''
    dH = len(s)
    for i in range(dH):
        if s[i] == t[i]:
            dH -= 1
    return dH

def ApproxPatternMatching(pattern,text,d):
    '''returns all starting positions where Pattern appears as a substring of Text with at most d mismatches'''
    k = len(pattern)
    n = len(text)-k+1
    pattern_i = []
    for i in range(n):
        if HammingDistance(text[i:i+k], pattern)<= d:
            pattern_i.append(i)
    return pattern_i

def ApproxPatternCount(pattern,text,d):
    pattern_i = ApproxPatternMatching(pattern,text,d)
    return len(pattern_i)

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

def FreqKmersWithMismatch(text, k, d):
    '''returns all most frequent k-mers with up to d mismatches in Text. Updated to include reverse complements as well.'''
    n = len(text)-k+1
    p_list = []
    most_freq = 1
    l = []
    for i in range(4**k):
        l.append(NumberToPattern(i, k))
    for i in range(4**k):
        rcl = Complement(l[i])
        if (ApproxPatternCount(l[i], text,d) + ApproxPatternCount(rcl, text,d))> most_freq:
            most_freq = (ApproxPatternCount(l[i], text,d) + ApproxPatternCount(rcl, text,d))
    for i in range(4**k):
        rcl = Complement(l[i])
        if (ApproxPatternCount(l[i], text,d) + ApproxPatternCount(rcl, text,d)) == most_freq:
            p_list.append(l[i])
    return " ".join(p_list)

window_start = 3764856
window_stop = 3765300

#WARNING: current algorithm is inefficient and takes a very long time to run


#print MinimumSkewProblem(data)
print FreqKmersWithMismatch((data[window_start:window_stop]), 9, 1)
