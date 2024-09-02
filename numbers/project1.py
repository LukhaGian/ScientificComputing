"""
    Project1
    Name: Gianmaria
    Surname: Lucca
    MAT: 241440
"""
def number2base_rep(n, b):
    """
    Function that takes 2 integers and returns
    the representation in base b of n
    """
    rep = "" # Initialize the string
    while n > 0:  
        rem = n % b
        n = n // b
        rep = str(rem)+rep
    return rep   


def admissible(n, b):
    """
    Function that takes integer n and returns
    if n is b-admissible
    """
    # First, generate the n_b (string)
    n_b = number2base_rep(n, b)
    is_admissible = True
    # Now, we know that a substring z of n_b s.t
    # it exists a neighbouring substring y of z,
    # of course len(z) <= len(n_b) // 2
    upper = len(n_b) // 2
    # Now, we check for every substring of length i <= upper
    i = 1
    while (i <= upper):
        # we create a list that contains all the substrings of length i
        list = []
        head = 0 # head of the substring
        tail = i # tail of the substring
        while tail <= len(n_b):
            list.append(n_b[head:tail])
            head += 1
            tail += 1
        # Given the list of substrings of length i,
        # we check if there are equal and contiguous substrings
        # OBS: contiguous substrings are far away i-1 spaces in the list (meaning between 2 contiguous substrings there are i-1 other substrings)
        j = 0
        jump = i
        while ((j+jump) <= len(list)-1):
            if (list[j] == list[j+jump]): # if the contiguous substrings are the same
                is_admissible = False
                # break both cycles
                j = len(list)
                i = upper+1
            else:
                j += 1                   # otherwise, check the next substrings        
        i += 1
    return is_admissible    
 

def count_admissible(b, start, end):
    """
    Function that takes three integers and
    returns the number of b-admissible numbers n with start <= n < end.
    """
    c = 0 # Initialize the counter
    for k in range(start, end):
        if admissible(k, b) == True:
            c += 1    
    return c # our counter


def count_admissible_width(b, width):
    """
    Function that takes two integers and returns
    the number of b-admissible numbers n 
    whose b-representation has exactly width digits
    """
    # Remark: given b and width, we have [b^(width-1)]*[b-1] possible strings which represent 
    # in base b all the numbers that have exactly width digits
    
    # First of all, we create a list which contains all the possible values that a digit in the b-representation can have (so from 0 to b-1)
    digits_list = []
    for i in range(0,b):
        digits_list.append(str(i))
    # Now, we know that the left digit cannot be 0, so we use a for loop in order to have only the strings with 1,2,...,b-1 on the left
    comb_list = [] # Initialize the list which contains all the possible numbers
    for x in range(1,b):
        left = str(x)
        All_Width_Length(digits_list, width-1, left, comb_list) # width - 1, since the first value is occupied by left
    # Now, given the list comb_list, we check for each element if it's admissible or not
    admissible_counter = 0
    for i in range(len(comb_list)):
        if (admissible_string(comb_list[i])): # We call the function that returns if the number, which representation in base b is 
                                              # comb_list[i], is admissible or not admissible
            admissible_counter += 1
            #print(f"The element {comb_list[i]} is admissible.")
    return admissible_counter

############### Support functions for count_admissible_width(b, width)
###############
def All_Width_Length(set, k, left, comb_list):
    """
    Wrapper method for All_Width_Length_Rec(set, prefix, n, k, left, comb_list)
    """
    n = len(set)
    All_Width_Length_Rec(set, "", n, k, left, comb_list)
    return
 
def All_Width_Length_Rec(set, prefix, n, k, left, comb_list):
    """
    The main recursive method
    to append to comb_list all possible
    strings of length k, 
    given -) a set of elements, 
          -) a (recursive) prefix, 
          -) the most left prefix (left)
          -) the list itself
     """
    # Base case: k is 0,
    # We append the string to comb_list
    if (k == 0) :
        #print(left+prefix)
        comb_list.append(left+prefix)
        return
 
    # One by one add all elements
    # from set and recursively
    # call for k equals to k-1
    for i in range(n):
 
        # Next element of input added
        newPrefix = prefix + set[i]
         
        # k is decreased, because
        # we have added a new element
        All_Width_Length_Rec(set, newPrefix, n, k - 1, left, comb_list)
 
def admissible_string(n_b):
    """
    Function that takes the number n in base b n_b and returns
    if n is b-admissible
    """
    # It's almost identical to admissible(n, b)
    is_admissible = True
    # Now, we know that a substring z of n_b s.t
    # it exists a neighbouring substring y of z,
    # of course len(z) <= len(n_b) // 2
    upper = len(n_b) // 2
    # Now, we check for every substring of length i <= upper
    i = 1
    while (i <= upper):
        # we create a list that contains all the substrings of length i
        list = []
        head = 0 # head of the substring
        tail = i # tail of the substring
        while tail <= len(n_b):
            list.append(n_b[head:tail])
            head += 1
            tail += 1
        # Given the list of substrings of length i,
        # we check if there are equal and contiguous substrings
        # OBS: contiguous substrings are far away i-1 spaces in the list (meaning between 2 contiguous substrings there are i-1 other substrings)
        j = 0
        jump = i
        while ((j+jump) <= len(list)-1):
            if (list[j] == list[j+jump]): # if the contiguous substrings are the same
                is_admissible = False
                # break both cycles
                j = len(list)
                i = upper+1
            else:
                j += 1                   # otherwise, check the next substrings        
        i += 1
    return is_admissible    
###############
###############


def largest_multi_admissible(L, start, end):
    """
    Function that takes a list of bases b L and
    returns the biggest number n s.t. start <= n < end
    which is b-admissible for all the bases b in L
    It returns None if no such number exist
    """
    largest_admissible = None
    for i in range(start, end): # We check for all the numbers n
        for j in range(len(L)): # We check all the bases in list L
            if (admissible(i, L[j]) == False): 
                break # We break the cycle if for the j-th base i is not b-admissible
            elif j == len(L)-1:
                largest_admissible = i # otherwise, if we're at the final base, we update the value
    return largest_admissible            


"""
if __name__ == "__main__":
    #count_admissible_width(5,10)
    #print(count_admissible_width(3,4))
    for k in range(1,5):
        print(count_admissible(5, 10**k, 10**(k+1)))
        print(count_admissible_width(3, k))
        print(largest_multi_admissible([3, 5, 7, 10], 1, 10**k))
        print("Done")

    if number2base_rep(42, 2) != "101010": print("error in number2base_rep")
    if admissible(256, 3) or not admissible(1202102012021, 10): 
        print("error in admissible")
    if count_admissible(5, 1, 10**4) != 3153: print("error in count_admissible")
    if count_admissible_width(3, 3) != 8 : print("error in count_admissible_width")
    if largest_multi_admissible([3, 4], 1, 10**4) != 7880:
        print("error in largest_multi_admissible")
    print("unittest completed.")    
"""  

