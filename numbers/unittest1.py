""" Test the functions from project1.py. """
from project1 import number2base_rep, admissible, count_admissible, count_admissible_width, largest_multi_admissible
from time import time

if number2base_rep(42, 2) != "101010": print("error in number2base_rep")
if admissible(256, 3) or not admissible(1202102012021, 10): 
    print("error in admissible")
if count_admissible(5, 1, 10**4) != 3153: print("error in count_admissible")
if count_admissible_width(3, 3) != 8 : print("error in count_admissible_width")
if largest_multi_admissible([3, 4], 1, 10**4) != 7880:
      print("error in largest_multi_admissible")
print("unittest completed.")
"""
print("First test")
for k in range(1,8):
        print("k = ",k)
        tic = time()
        n = count_admissible(5, 10**k, 10**(k+1))
        elapsed = time() - tic
        print(f"The value is " ,n)
        print(f"Elapsed time = {elapsed:6.4e}.")
input("Finished")
"""
print("Second test")
for k in range(1,14):
        tic = time()
        n = count_admissible_width(3, k)
        elapsed = time() - tic
        print(f"k = {k}, The value is {n}, Elapsed time = {elapsed:6.4e}.", )
input("Finished")
print("Third test")
for k in range(1,8):
        print("k = ", k)
        tic = time()
        n = largest_multi_admissible([3, 5, 7, 10], 1, 10**k)
        elapsed = time() - tic
        print(f"The value is " ,n)
        print(f"Elapsed time = {elapsed:6.4e}.")