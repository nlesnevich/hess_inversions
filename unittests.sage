"""
This is the testing function file for
 - bexpansion.sage
"""

import random

#bexpansion
attach("bexpansion.sage")

"""
Code to compute the b_k(q) expansion of h-inversion polynomials
and check whether it is stongly q-log concave
i.e., whether b_k^2(q) - b_{k-1}(q)b_{k+1}(q) has all on-negative coefficients
"""

R.<q> = PolynomialRing(QQ) #Sets the polynomial ring globally

def test_makePoset(hfam, S):
    """
    hessenberg func  -hfam
    list of inversions -S

    returns the poset P associated to these two
    """
    n = hfam[-1]
    ints_n = range(1,n+1) #[1,...,n]
    rels = [(i,j) for i in ints_n for j in range(i+1,hfam[i-1]+1)] #edges in digraph
    direls = [(j,i) if (i,j) in S else (i,j) for (i,j) in rels] #flips those in S
    return Poset((ints_n,direls)) #returns poset on [1,...,n]

def test_perms(hfam,S):
    """
    returns the set of permutations whose hfam-inversion set is S

    these are the inverses of the linear extensions of the associated poset
    """
    P = makePoset(hfam,S)
    return [Permutation(L).inverse() for L in P.linear_extensions()]

def test_bfam_dict(hfam, S):
    """
    returns a dictionary where 
    k -> (list of perms w in S_{hfam(m)} where w(hfam(m)) = k)
    other words:
    k -> B_k(S,hfam(m))
    """
    ret_dict = dict()
    if len(S) == 0:
        return ret_dict
    m = max([inv[0] for inv in S])
    hm = hfam[m-1]
    hfam = [min([hm,hfam[i]]) for i in range(hm)]
    ##Now everything is the right size
    AllPerms = perms(hfam,S)
    for i in range(1,hm+1):
        Bkset = [p for p in AllPerms if p[-1] == i]
        if len(Bkset) > 0:
            ret_dict[i] = Bkset
    return ret_dict

def test_bfam_count(hfam,S):
    """
    Returns a list of tuples (k,b_k)
    """
    BD = bfam_dict(hfam,S)
    return [(k,len(BD[k])) for k in BD.keys()]

def test_length_genfunc(perms):
    """
    returns the length generating function for a set of permutations
    """
    retpoly = 0
    if len(perms) == 0:     #necessary because we pull "n" from perms[0]
        return retpoly
    n = len(perms[0])
    for p in perms:
        ell = p.length()
        retpoly += q^(ell)      #q is defined globally
    return retpoly

def test_bfam_poly(hfam,S):
    """
    returns dictionary 
    k -> b_k(q)
    """
    BD = bfam_dict(hfam,S)
    retdict = dict()
    for i in BD.keys():
        retdict[i] = length_genfunc(BD[i])
    return retdict

def test_strongly_logconcave(hfam,S):
    """
    T/F check if the b_k(q) sequence is strongly log concave
    """
    BD = bfam_poly(hfam,S)
    coefs = sorted(BD.keys())
    polylist = []
    for i in range(1,len(coefs)-1):
        polylist += [(BD[coefs[i]])^2 - BD[coefs[i-1]]*BD[coefs[i+1]]]
    #for i in range(len(polylist)):
        #print(polylist[i])
    for poly in polylist:
        for co in list(poly):
            if co < 0:
                return False
    return True 
        
def test_all_admit(hfam):
    """
    For a given hfam, returns a list of all admissable inversion sets
    """ 
    print()
    n = hfam[-1]
    retlist = []
    ints_n = range(1,n+1) #[1,...,n]
    rels = [(i,j) for i in ints_n for j in range(i+1,hfam[i-1]+1)]
    G = Graph(rels)
    acyclics = G.acyclic_orientations()  #much quicker than looping over all permutations
    for P in acyclics:
        S = [(e[1],e[0]) for e in P.edges() if e[0] > e[1]]
        retlist += [S]
    return retlist

def test_conjecture():
    """
    For a particular hfam, checks the log-concavity condition for all inversion sets
    """
    return None

def test_conj():
    """
    --Prints results
    Generates all hfam on(n), then checks the
    log-concavity conditions for each inversion set
    for each of those hfam
    """ 
    return None

def test_all_hfam():
    """
    returns all (connected) hess functions for a particular n
    """
    print("The function all_hfam(n) returns all (connected) hessenberg functions of size n, for example:")
    print("all_hfam(3) = " + str(all_hfam(3)) )
    print("should be 233 and 333")
    print("and")
    print("The function all_hfam(n) returns all (connected) hessenberg functions of size n, for example:")
    print("all_hfam(4) = " + str(all_hfam(4)) + "\n")
    print("should be 2344 and 2444 and 3344 and 3444 and 4444")
    return None


def test_dyck_to_hess(word):
    """
    turns a dyck word in to the hessenberg function (in the obvious way)
    """
    print("dyck_to_hess([1,1,1,1,0,0,0,0]) = " + str(dyck_to_hess([1,1,1,1,0,0,0,0])))
    print("should be 4444")
    print("dyck_to_hess([1,1,1,0,1,0,0,0]) = " + str(dyck_to_hess([1,1,1,0,1,0,0,0])))
    print("should be 3444")
    print("dyck_to_hess([1,1,1,0,0,0,1,0]) = " + str(dyck_to_hess([1,1,1,0,0,0,1,0])))
    print("should be 3334")
    print("dyck_to_hess([1,0,1,0,1,0,1,0]) = " + str(dyck_to_hess([1,0,1,0,1,0,1,0])))
    print("should be 2344")
