"""
Code to compute the b_k expansion of h-inversion polynomials
"""

R.<q> = PolynomialRing(QQ)

def makePoset(hfam, S):
    """
    hessenberg func hfam
    list of inversions S

    returns the poset P associated to this
    """
    n = hfam[-1]
    ints_n = range(1,n+1) #[1,...,n]
    rels = [(i,j) for i in ints_n for j in range(i+1,hfam[i-1]+1)] #edges in digraph
    direls = [(j,i) if (i,j) in S else (i,j) for (i,j) in rels] #flips those in S
    return Poset((ints_n,direls)) #returns poset on [1,...,n]

def perms(hfam,S):
    """
    returns the permutations whose hfam-inversion set is S

    these are the inverses of the linear extensions of the associated poset
    """
    P = makePoset(hfam,S)
    return [Permutation(L).inverse() for L in P.linear_extensions()]

def bfam_dict(hfam, S):
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
    BD = bfam_dict(hfam,S)
    return [(k,len(BD[k])) for k in BD.keys()]

def length_genfunc(perms):
    retpoly = 0
    if len(perms) == 0:
        return []
    else:
        n = len(perms[0])
        for i in range((n*(n-1))//2+1):
            lenset = [p for p in perms if p.length() == i]
            if len(lenset) > 0:
                retpoly += len(lenset)  * q^(i)
    return retpoly

def bfam_poly(hfam,S):
    BD = bfam_dict(hfam,S)
    retdict = dict()
    for i in BD.keys():
        retdict[i] = length_genfunc(BD[i])
    return retdict

def strongly_logconcave(hfam,S):
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
    n = hfam[-1]
    retlist = []
    ints_n = range(1,n+1) #[1,...,n]
    rels = [(i,j) for i in ints_n for j in range(i+1,hfam[i-1]+1)]
    G = Graph(rels)
    acyclics = G.acyclic_orientations()
    for P in acyclics:
        S = [(e[1],e[0]) for e in P.edges() if e[0] > e[1]]
        retlist += [S]
    return retlist

def conjecture(hfam):
    allS = all_admit(hfam)
    for S in allS:
        if not strongly_logconcave(hfam,S):
            print(S)
            return False
    return True
