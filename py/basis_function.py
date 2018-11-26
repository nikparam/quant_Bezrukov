import numpy as np
np.set_printoptions(precision=12)

def fact2(n):
     if n <= 0:
         return 1
     else:
         return n * fact2(n-2) 

class BasisFunction(object):
    ''' A class that contains all our basis function data
        Attributes:
        origin: array/list containing the coordinates of the Gaussian origin
        shell:  tuple of angular momentum
        exps:   list of primitive Gaussian exponents
        coefs:  list of primitive Gaussian coefficients
        norm:   list of normalization factors for Gaussian primitives
    '''
    def __init__(self,origin=[0.0,0.0,0.0],shell=(0,0,0),exps=[],coefs=[]):
        self.origin = np.asarray(origin)
        self.shell = shell
        self.exps  = exps
        self.coefs = coefs
        self.norm = None
        self.normalize()

    def normalize(self):
        ''' Routine to normalize the basis functions, in case they
            do not integrate to unity.
        '''
        l,m,n = self.shell
        L = l+m+n
        # self.norm is a list of length equal to number primitives
        # normalize primitives first (PGBFs)
        self.norm = np.sqrt(np.power(2,2*(l+m+n)+1.5)*
                        np.power(self.exps,l+m+n+1.5)/
                        fact2(2*l-1)/fact2(2*m-1)/
                        fact2(2*n-1)/np.power(np.pi,1.5))

        # now normalize the contracted basis functions (CGBFs)
        # Eq. 1.44 of Valeev integral whitepaper
        prefactor = np.power(np.pi,1.5)*\
            fact2(2*l - 1)*fact2(2*m - 1)*fact2(2*n - 1)/np.power(2.0,L)

        N = 0.0
        num_exps = len(self.exps)
        for ia in range(num_exps):
            for ib in range(num_exps):
                N += self.norm[ia]*self.norm[ib]*self.coefs[ia]*self.coefs[ib]/\
                         np.power(self.exps[ia] + self.exps[ib],L+1.5)

        N *= prefactor
        N = np.power(N,-0.5)
        for ia in range(num_exps):
            self.coefs[ia] *= N

def E(i,j,t,Qx,a,b):
    ''' Recursive definition of Hermite Gaussian coefficients.
        Returns a float.
        a: orbital exponent on Gaussian 'a' (e.g. alpha in the text)
        b: orbital exponent on Gaussian 'b' (e.g. beta in the text)
        i,j: orbital angular momentum number on Gaussian 'a' and 'b'
        t: number nodes in Hermite (depends on type of integral, 
           e.g. always zero for overlap integrals)
        Qx: distance between origins of Gaussian 'a' and 'b'
    '''
    p = a + b
    q = a*b/p
    if (t < 0) or (t > (i + j)):
        # out of bounds for t  
        return 0.0
    elif i == j == t == 0:
        # base case
        return np.exp(-q*Qx*Qx) # K_AB
    elif j == 0:
        # decrement index i
        return (1/(2*p))*E(i-1,j,t-1,Qx,a,b) - \
               (q*Qx/a)*E(i-1,j,t,Qx,a,b)    + \
               (t+1)*E(i-1,j,t+1,Qx,a,b)
    else:
        # decrement index j
        return (1/(2*p))*E(i,j-1,t-1,Qx,a,b) + \
               (q*Qx/b)*E(i,j-1,t,Qx,a,b)    + \
               (t+1)*E(i,j-1,t+1,Qx,a,b)

def overlap(a,lmn1,A,b,lmn2,B):
    ''' Evaluates overlap integral between two Gaussians
        Returns a float.
        a:    orbital exponent on Gaussian 'a' (e.g. alpha in the text)
        b:    orbital exponent on Gaussian 'b' (e.g. beta in the text)
        lmn1: int tuple containing orbital angular momentum (e.g. (1,0,0))
              for Gaussian 'a'
        lmn2: int tuple containing orbital angular momentum for Gaussian 'b'
        A:    list containing origin of Gaussian 'a', e.g. [1.0, 2.0, 0.0]
        B:    list containing origin of Gaussian 'b'
    '''
    l1,m1,n1 = lmn1 # shell angular momentum on Gaussian 'a'
    l2,m2,n2 = lmn2 # shell angular momentum on Gaussian 'b'
    S1 = E(l1,l2,0,A[0]-B[0],a,b) # X
    S2 = E(m1,m2,0,A[1]-B[1],a,b) # Y
    S3 = E(n1,n2,0,A[2]-B[2],a,b) # Z
    return S1*S2*S3*np.power(np.pi/(a+b),1.5)

def S(a,b):
    '''Evaluates overlap between two contracted Gaussians
       Returns float.
       Arguments:
       a: contracted Gaussian 'a', BasisFunction object
       b: contracted Gaussian 'b', BasisFunction object
    '''
    s = 0.0
    for ia, ca in enumerate(a.coefs):
        for ib, cb in enumerate(b.coefs):
            s += a.norm[ia]*b.norm[ib]*ca*cb*\
                     overlap(a.exps[ia],a.shell,a.origin,
                     b.exps[ib],b.shell,b.origin)
    return s

#HYDROGEN
#S   3
  #1      3.42525091             0.15432897
  #2      0.62391373             0.53532814
  #3      0.16885540             0.44463454
Hexps   = [3.42525091, 0.62391373, 0.16885540] 
Hcoefs  = [0.15432897, 0.53532814, 0.44463454]

#OXYGEN
#S   3
  #1    130.7093200              0.15432897
  #2     23.8088610              0.53532814
  #3      6.4436083              0.44463454
#L   3
  #1      5.0331513             -0.09996723             0.15591627
  #2      1.1695961              0.39951283             0.60768372
  #3      0.3803890              0.70011547             0.39195739
Oexps1 = [130.7093200, 23.808610, 6.4436083]
Oexps2 = [5.0331513, 1.1695961, 0.3803890]

Ocoefs1 = [0.15432897, 0.53532814, 0.44463454]
Ocoefs2 = [-0.09996723, 0.39951283, 0.70011547]
Ocoefs3 = [0.15591627, 0.60768372, 0.39195739]

#water = """
#0 1
#O  0.000000000000 -0.143225816552 0.000000000
#H  1.638036840407  1.136548822547 0.000000000
#H -1.638036840407  1.136548822547 0.000000000
#"""
oxygen_s1 = BasisFunction(origin=[0.0, -0.143225816552, 0.0],shell=(0,0,0),exps=Oexps1,coefs=Ocoefs1) # separate s-function 
oxygen_ls = BasisFunction(origin=[0.0, -0.143225816552, 0.0],shell=(0,0,0),exps=Oexps2,coefs=Ocoefs2) # s-function in L
oxygen_lp1 = BasisFunction(origin=[0.0, -0.143225816552, 0.0],shell=(1,0,0),exps=Oexps2,coefs=Ocoefs3) # p-function in L
oxygen_lp2 = BasisFunction(origin=[0.0, -0.143225816552, 0.0],shell=(0,1,0),exps=Oexps2,coefs=Ocoefs3) # p-function in L
oxygen_lp3 = BasisFunction(origin=[0.0, -0.143225816552, 0.0],shell=(0,0,1),exps=Oexps2,coefs=Ocoefs3) # p-function in L

hydrogen1 = BasisFunction(origin=[ 1.638036840407, 1.136548822547, 0.0], shell=(0,0,0), exps=Hexps, coefs=Hcoefs)
hydrogen2 = BasisFunction(origin=[-1.638036840407, 1.136548822547, 0.0], shell=(0,0,0), exps=Hexps, coefs=Hcoefs)

bfs = [oxygen_s1, oxygen_ls, oxygen_lp1, oxygen_lp2, oxygen_lp3, hydrogen1, hydrogen2 ]

tmp = []
for ia in range(len(bfs)):
    for ib in range(len(bfs)):
        tmp.append( S(bfs[ia], bfs[ib]) )

S = np.array(tmp).reshape((7, 7))
print("S:\n{0}".format(S))
