"""
    Project2
    Name: Gianmaria
    Surname: Lucca
    MAT: 241440
"""
def kmer_search(S, k, freq, x):
    """
    Function that given a string S (of combination of A,C,G,T),
    the integers k,freq and x in [0,k-1] retruns
    2 lists, L1 and L2:
    -) L1 is a list of lists, where for each k-mer
       there is a list of all its positions in S
    -) L2 is a list where at entry i is stored
       the number of point-x mutations (in S) of the k-mer
       in the i-th position    
    """
    L1 = []
    L2 = []

    # We'll move through the whole string, considering the substrings of length k, and we check, using the Rabin Karp Hash Search, the occourencies
    # of the substring
    # if occurencies >= freq ==> We append the whole list of occurencies to L1
    sub_list = [] # support list that tells us which positions has been already evaluated
    for i in range(0,len(S)-k+1): #position of the first element of the substring
      if ((i in sub_list) == False):
         sub = S[i:i+k] # the substring
         occurencies = rabin_karp_search(S,sub)
         #print(f"The occurencies are: {occurencies}")
         sub_list = sub_list+occurencies # we update the already-visited substrings
         # We check if occurencies >= freq
         if (len(occurencies) >= freq):
               # We append the list occurencies to L1
               L1.append(occurencies)
               # Now, we search for the point x-mutation of the last entry of L1, calling the support function x_mutation_counter
               L2 = x_mutation_counter(S, sub, x, L2)

    return L1,L2


def x_mutation_counter(S, sub, x, L2):
    """
    Function that appends the number of x_mutations' positions of the substring sub in the string S
    """
    # We create a list of all the possible nucleotides
    no_mutations = 0
    nuc = ['A', 'C', 'G', 'T']
    nuc.remove(sub[x]) # We remove from the list the nucleotide in sub at position x
    for j in nuc:
        # since we cannot modify directly the string, we actually create a list and then convert it back to string
        mutant_sub_list = list(sub)
        mutant_sub_list[x] = j
        mutant_sub = ''.join(mutant_sub_list) # now we have our mutant substring (type string)
        # We check, using the Rabin Karp algorithm, the occurencies of the mutant string in S
        mutant_occ = rabin_karp_search(S, mutant_sub)
        no_mutations += len(mutant_occ) # update the number of occurencies

    L2.append(no_mutations) #append the total value to L2

    return L2


#########################################
""" implement Rabin Karp algorithm """
base = 4
mod = 100049
def extend_by_one(h, c):
    return (base * h + ord(c)) % mod

def remove_left(h, c, bstar):
    return (h - bstar * ord(c)) % mod

def hash_value(s):
    h = 0
    for c in s: h = extend_by_one(h, c)
    return h

def rabin_karp_search(s, p):
    """ return all the occurances of p in s using the Rabin-Karp algorithm """
    indeces = []
    n, m = len(s), len(p)
    bstar = pow(base, m-1, mod)
    hp = hash_value(p)
    hs = hash_value(s[:m])
    for i in range(n-m+1):
        if hp == hs:
            if s[i:i+m] == p: indeces.append(i)
        if i < n-m:
            hs = remove_left(hs, s[i], bstar)
            hs = extend_by_one(hs, s[i+m])
    return indeces
#########################################


def spet_location(S, k, p):
    """
    function that returns the starting position q
    of a k-mer K in S,
    such that
    -) p-2k <= q < p-k
    -) the k-mer can't appear more than 5 times in S
    -) the proportion in K of "C" and "G" must be between 35% and 65% 
    """
    q = None # Initialization
    q_frequency = 0

    for i in range(p-2*k, p-k):
        sub = S[i:i+k] # the substring we're working with
        if ((0.35*k) <= ((sub.count('C')+sub.count('G'))) <= (0.65*k)): # Proportion check
            frequency = len(rabin_karp_search(S, sub)) # We search the number of times the k-mer appears in S 
            if (frequency <= 5):
                if ((frequency <= q_frequency) or (q is None)):
                    q = i # here we substitute the value, 
                          # since if frequency == q_frequency we take the closest k-mer, which is the one we have in this loop (sub)
                    q_frequency = frequency # we update the value of the frequency

    return q    



def get_my_pqs():
    # First we call the cleaning function clean_Fasta()
    # the file sequence.fasta should be in the same directory as this file
    dna = clean_Fasta()
    k = 40
    list = []
    a = (156, spet_location(dna,k,156))
    list.append(a)
    b = (999, spet_location(dna,k,999))
    list.append(b)
    c = (130800, spet_location(dna,k,130800))
    list.append(c)
    return list

def clean_Fasta():
    """
    Function that cleans the data in the file sequence.fasta and returns the string of the DNA
    """
    string = ""
    with open("sequence.fasta") as file:
        for item in file: 
            stripped = item.rstrip() # We copy all the non empty elements of a given row
            if (stripped != ""):
                if (stripped[0] != ">" and stripped[0] != ";"): # So the row is not the first one, which in general is a comment
                    # We remove all the values in the row that are not nucleotides (A,G,C,T)
                    clean_stripped = ""
                    for q in range(len(stripped)):
                        if (stripped[q] in ['A','G','C','T']) == True:
                            clean_stripped+=stripped[q]
                    string += clean_stripped           
    return string


if __name__ == "__main__":
    from time import time
    s = "GTCGTATTTGAGGTCCACCGCCCCGGCCATGCCGGGTGTTGGCGTCGCCGTCTAGGGGCGCAGTATACGGCTAAGATTTGCCGGATACAGAGTATAACTTTCGAGATTGAGGATGTGCATACTAAATGCACAAACCTGTCACTGAACGACAGACCCGAAACGCGAGACCCTCTAATATCACTGTCCAAAGATCCGTATGTAAATACACGCTTCACACGTCTTTTCTTAATGGGCTTTCGCGCTGCGCATGCAATCGGTTTGCGGGGGTCGGACTAGCCTACTGGTTCGGAAAATTCACGGGTAACACAAATCGGTACCGACCAAGAACATCCCAATATAGATCGCATTCCACACAAGAAGTGTCAGGTTCTATTATGAATAGTGGAGGTAAAGGAAGCATGCCGATTACTCGTGTCCCCGATGGGGACTAGACTGCTCCTGCTCTATCACGGTTTGCCTACCAGGAGTCGTGCCGGTAGAGCGCGCTCTACACTATTCACATGTTCCTATCATCATACGCTCTCTGGATTATTATCCCACCCTCCCGGATGTAGTCACTAGGAACTGGGTTAAAACAGGTCGAACAGCGCTACAATCTTGCCTGATGGTGATGTCAGTTCCCACCATATTACGACACTATTCTAATTCGAGATATATATTTCGCTGTGTGCGGTCACTGGAATTAGAGGTGCCGTCAACGCTCCATATTTCAAAGAACGATCAACTAGGAATGCTTGACATTCGACTAGTGCTGTATCAATAATTCGATAGCGAAGACCGTTCGGGTAAATTGCTGGCCCTTCTGTCGCCGCCCAAGCACAATCACTTGGCTGCATCAAAAAGGTTGGCAACGGAACTGCTTTATTGTAAAGGCCGCCTAACCGATAGGTTGCAAAAGGACATGAGCGACGGGAAAGCGTGACTGGTATCCGGAGTGCCAAGCATAGGGTAAACAGAGAGGTGTTGTCGGAACCGTACACTGCCAGTCCCCACATCCGATGCGGGCGAGCCAAATTAACGGAGCCAAAGCCTGCCGCTTACTTATACTTGCTACAGAAAAAGGGTGTTTCCCAGCGTCCCCGCAGGTGTGTATGTAATTTCTGACTTCTCTGCCCGGAGCGGAGTGGGTGGCCGTGACGTCAGCAGCGTACCAGACTCGGAATGAGGGGTTCTCCTACCATTTTTCGCGGGTGGACGGCTGATCGCTTCTCTGACGTCGATAATGAGCGCTGCGCTGGATTCAAACATCCTCTGATTCCCTTGACCCTGCGCTAGGCGTCGACTGCCCTACTTGCGGTGAGTAGAGTCTGGGGATAGTCTCATGTACCCTCTAAATCTTGATATTAGAGCGGCCGAATACGGTCAAAACGTCTGGGGTGCATCCGTACTCCACTCTCACAGACCGTAAAGCTTAGTGGACCTGACCAGTTGGCTTCACCGGCACTCGTAAGTCGCTGTCCACGGGGACCTTGCAACCACATCTGATCCCGTCCCGCCTCTCGACTACGCCGAGTCACCTTAACCTGGTTCCACCCCCGAACGCGCCTTACGTCTCTTATGGACCGGCATTCTGTTCACGAGACCAAGCCAGGGTGACTAAATAACTAATTACGCACTTCTCAAAAAATGACACCAGGGATTAATAGAAAAGAATCCTCGCAGGCGGCAACTGCTACTACAGCTGAACCCATTTCTACACTCAGTTCTGGCCATCACGTGGTTTAGGTGAGGAGTCAAGTAACGAGCGGGAGGCAGACAGCGGTTTCGCCTTAGAGTTCTCGTGCGCTAGGCCTCGACACCAAGGGTATTGTATTACACTGTTCCTCCGGTGCCCCATTCCCTCGTTCCGTCTCATCCTTCTCTCGGCACGTAGATTGGACTTAACCAGACTCACGAGCGTATGAGTAGTGTACTTTCATGCATAAGCTGACATCAAAATATCAGGCGTTTCCGCGGCCGATAATCAGGAGGTCAGGGGATAGGGGCACTTCGGTCCCGCATTTGTGATCGCCATGTTATCGATCAAATCAGATAATATTATGCGCGTGACTCTAAGTAAGCGAAATTGCAGATAGGCGTTAAGTGGGCAGTAAAACCAAGCCACAGCGTAGTCGAACCCAAGTCCGACCGGTTGTAAGGATAAGAAGTTATACGGGTTGGGACCTACCTATTCCTGCGGCGATGAATGGGGTATGACACACCCAGGTGCTTTCGACTTGGTCCTGCACCGCGCAACGAACGTCCAGCAAAATGTTCATGAAGTTGTCGGGCTCACTAGAACCAGTAATCGGAATTGCCCAGGGTTCTACATGACAGTCGAAGGTCCAACCCTGGTAGCATGGCCCGGAGGGTGTGCAGAAGTTAACCCGCTAGAACGGTGGCGCATGTACATCTGTTAAAATTTATAGTCCAGTAATGAGGATTCTAAGAGGTCAGAAACATTGAAAGAAGTACGACCCTACCTAGTGGGAACCGGTCGGCGTCCGGGCACGTAACTTATCTTAGAGGACGACCATGCCTTACCCTTGCAAGCGATTTGCGGCATGCGGGCCGTACTGATGGGAACGTTCTAGTGTAAATGAGCGGAGACAGCTACTCTAGCAATCCAAAATCTTACTCAGTGTGGGTGGGAAGGCTCCAGCGGGAGAGCAGCCCGGGCTTCAGTTTTGGAACTAATCTACTCGTTTAGTCGAGGGGTTATCCTTGGAGCGTGCCCGCAGGAAACGGTATGATATATGGATATGGGAGTAGGATACCTCGTGTCGGGGCCAAGGGTCCTTTGCTACCAGCGACTCTTCGTATCTGTTACGGCTTCTTCCGGAAGGATATATGATGCACGAGAACCAGGTAGGGTGACGACATAACTAATTTAAAAGTCTACGCACGGATACCACGTCAGTTAGATGTCCAAGGAAAATTTAACACATCAGGCGTTATTCATAATCTTCTGGAGAGGGTGATGAGCGACGGGA"
    q = time()
    L1,L2 = kmer_search(s, 7, 3, 1)
    elapsed = time()-q
    print(L1)
    print(L2)
    print(elapsed)
    print(len(s))
    print(get_my_pqs())
    k = 40
    string = clean_Fasta()
    print(len(string))
    print(spet_location(string,k,156))
    third = spet_location(string, k, 130800)
    print(string[third:third+k])
    for p in range(2*k,len(string)+1):
        val = spet_location(string,k,p)
        if val is None: 
            print(p)
            break
        else: print(f"Ok for p = {p}")
    print(spet_location(string, k, 130800))   