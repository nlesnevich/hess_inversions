"""
Code to compute 
h_sequences(n) - all sequences h(n) such that h(i) > i
noninc_hsequences(n) - all non-increasing h(n) such that h(i) > i
"""


def h_sequences(n):
    """ 
    yields all sequences h of length n with h(i) > i
    note this is a generator object so if you want a list do list(sequences_h(n))
    """
    def rec(i, current):
        if i == n-2:
            yield current+[n,n]
        for v in range(i+2, n+1):   # allowed {i+1,...,n}
            yield from rec(i+1, current + [v])
    yield from rec(0, [])

def noninc_hsequences(n):
    """
    yields all sequences h of length n with h(i) > i that aren't hessenbergs
    """
    for hseq in h_sequences(n):
        if hseq != sorted(hseq):
            yield hseq
            

def all_admit(hfam):
    """
    For a given hfam, returns a list "S" of all admissable inversion sets
    and a list of all permutations "hS_perms" of inversions associated to the hfam for I_hfam(n,J(S))
    in a tuple: (S,P)
    tosses all of the ones with j(S) != n-1 since they'll be handled by others.
    we keep the "n" around because we need to check edges with j(S) + 1 = n in the expansion formula later.
    """ 
    n = hfam[-1]
    retlist = []
    ints_n = range(1,n+1) #[1,...,n]
    rels = [(i,j) for i in ints_n for j in range(i+1,hfam[i-1]+1)]
    G = Graph(rels)
    acyclics = G.acyclic_orientations()  #much quicker than looping over all permutations
    for P in acyclics:
        if P.out_degree(n-1) > 0 and P.out_degree(n) == 0:#makes sure j(S) = n-1 
            #first the list of S
            S = [(e[1],e[0]) for e in P.edges() if e[0] > e[1]]
            #now the list of P
            #we only want permutations for I_H(n,J(S)), so we have to get rid on n
            P.delete_vertex(n)
            the_poset = Poset(P)
            hS_perms = [Permutation(L).inverse() for L in the_poset.linear_extensions()]
            yield (S,hS_perms)


def t_sequences(n):
    """
    returns all t-sequences for the relevant n (remember n = j(S) for these ones)
    the t(sigma) value is really all that we need, so this will be off of 
    possible t(sigmas) = [n,n-1,n-2,...,1]
    and returns
    [#{perms w/ t(sigma)=n},...,#{perms w/ t(sigma)=1}]     --note: many will be 0 but I'll keep the full length for polynomial conversions
    corresponding to the poly basis (in var x)
    [(x-n choose 0),(x-n+1 choose 1),...,(x-1 choose n-1)]
    if it's not log concave then it'll print that out in detail
    """
    print("This is going to check everything where j(S) = "+str(n-1))
    HSeq = noninc_hsequences(n)
    m = 0 #maxnontriv_achieved
    for hfam in HSeq:
        print(hfam)
        #now we have to handle t(sigma), which requires which indices k have (k,J(s)+1 = n) as an edge.
        checkset = [i for i in range(n-1) if hfam[i] == n]
        for S,Perms in all_admit(hfam):
            coef_sequence = [0]*(n-1)
            for sigma in Perms:
                #print("sigma is" + str(sigma))
                #print("check is " + str(checkset))
                t_sigma = max([sigma[i] for i in checkset])
                #print("tsig is " + str(t_sigma))
                coef_sequence[-1*t_sigma] += 1
            m = max([m,len([i for i in coef_sequence if i>0])])
            if False:
                print("S = "+str(S))
                print(coef_sequence)
                #print(Perms)
            if not is_log_concave(coef_sequence):
                print("TIME TO CHECK THIS THING:")
                print("    S = " + str(S))
                print("    Perms = " + str(Perms))
    return m

def tseq(hfam,perms):
    """
    lives under the assumption that j(S) = n-1
    """
    n = hfam[-1]
    checkset = [i for i in range(n-1) if hfam[i] == n]
    coef_sequence = [0]*(n-1)
    for sigma in perms:
        t_sigma = max([sigma[i] for i in checkset])
        coef_sequence[-1*t_sigma] += 1
    return coef_sequence

def t_a(tseq,hfam,m):
    n = hfam[-1]
    hm = hfam[m-1]
    aseq = [0]*(m+4)
    for k in range(m+1):
        ak = 0
        #print(k)
        for i in range(1,hfam[-1]):
            #print(str(hm-i-1) + " choose " + str(n-1-i-k))
            ak += binomial(hm-i-1,n-1-i-k)*tseq[-i]
        aseq[k] = ak
    return aseq




def is_log_concave(intlist):
    """
    takes in a list of integers
    returns if they are log concave
    """
    for i in range(1,len(intlist)-1):
        check_at_i = (intlist[i])*(intlist[i]) - (intlist[i-1])*(intlist[i+1])
        if check_at_i < 0:
            return False
    return True


