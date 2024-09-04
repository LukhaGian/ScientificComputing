""" Test the three functions from project2.py. """
from project2 import kmer_search, spet_location, get_my_pqs

k, x, freq = 7, 1, 3
with open("dna.txt") as f: S = f.read().rstrip()
L1, L2 = kmer_search(S, k, freq, x)
for i, l in enumerate(L1):
    pos = l[0]
    kmer = S[pos:pos+k]
    print(f"{kmer} found {len(l)} times at locations {l}.")
    print(f"{' '*k} there are {L2[i]} point-{x} mutations of it in S.")
if spet_location("GC"*100+"ATTA", 10, 50) is not None:
    print("Error due to too high GC content.")
if spet_location("AT"*100+"ATTA", 10, 50) is not None:
    print("Error due to too low GC content.")
if spet_location("ACGT"*100+"ATTA", 4, 50) is not None:
    print("Error due to too high frequency.")

for k in [4, 7, 14, 21, 28]:
    q = spet_location(S, k, 50)
    if q is None:
        print(f"For k = {k} no location found.")
    else:
        print(f"For k = {k} best location is {q} with sequence {S[q:q+k]}.")

L = get_my_pqs()
foundNone, foundnotNone = False, False
for p, q in L:
    if type(p) is not int: print("p should be an integer")
    if q is None:
        foundNone = True
    else: 
        if type(q) is not int: print("p should be an integer")
        foundnotNone = True
if not foundNone: print("One q should be None.")
if not foundnotNone: print("One q should be a number.")