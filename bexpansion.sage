"""
Code to compute the b_k(q) expansion of h-inversion polynomials
and check whether it is stongly q-log concave
i.e., whether b_k^2(q) - b_{k-1}(q)b_{k+1}(q) has all on-negative coefficients
"""

R.<q> = PolynomialRing(QQ) #Sets the polynomial ring globally

def makePoset(hfam, S):
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

def perms(hfam,S):
    """
    returns the set of permutations whose hfam-inversion set is S

    these are the inverses of the linear extensions of the associated poset
    """
    P = makePoset(hfam,S)
    return [Permutation(L).inverse() for L in P.linear_extensions()]

def bfam_dict(hfam, S):
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

def bfam_count(hfam,S):
    """
    Returns a list of tuples (k,b_k)
    """
    BD = bfam_dict(hfam,S)
    return [(k,len(BD[k])) for k in BD.keys()]

def length_genfunc(perms):
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

def bfam_poly(hfam,S):
    """
    returns dictionary 
    k -> b_k(q)
    """
    BD = bfam_dict(hfam,S)
    retdict = dict()
    for i in BD.keys():
        retdict[i] = length_genfunc(BD[i])
    return retdict

def strongly_logconcave(hfam,S):
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
        
def all_admit(hfam):
    """
    For a given hfam, returns a list of all admissable inversion sets
    """ 
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

def conjecture(hfam):
    """
    For a particular hfam, checks the log-concavity condition for all inversion sets
    """
    allS = all_admit(hfam)
    for S in allS:
        if not strongly_logconcave(hfam,S):
            print(S)
            return False
    return True

def conj(n):
    """
    --Prints results
    Generates all hfam on (n), then checks the
    log-concavity conditions for each inversion set
    for each of those hfam
    """
    for hess in all_hfam(n):
        print(str(hess) + " --> " + str(conjecture(hess)))
        if not conjecture(hess):
            print("!!!!!COUNTEREXAMPLE!!!!!!")

def all_hfam(n):
    """
    returns all (connected) hess functions for a particular n
    """
    Dwords = [DyckWord([1] + list(d) + [0]) for d in DyckWords(n-1)]
    return [dyck_to_hess(d) for d in Dwords]

def dyck_to_hess(word):
    """
    turns a dyck word in to the hessenberg function (in the obvious way)
    """
    area = list(word.to_area_sequence())
    area.reverse()
    for i in range(len(area)):
        area[i] += i+1
    return area
