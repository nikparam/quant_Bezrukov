from __future__ import print_function
import scipy.special as sp
import numpy as np

def fact2(n):
     if n <= 0:
         return 1
     else:
         return n * fact2(n-2)

class GTO(object):
    def __init__(self, origin = [0.0,0.0,0.0], shell = (0,0,0), exp = None, coeff = None):
        self.origin = np.asarray(origin)
        self.shell = shell
        self.exp = exp
        self.coeff = coeff

class STO(object):
    def __init__(self, orbs = [], ngauss=None):
        self.orbs = orbs
        self.ngauss = len(orbs)

def od_overlap(i, j, alpha, beta, xa, xb):
    
    eta = alpha + beta
    xpb = alpha / eta * (xa - xb)
    xpa = -beta / eta * (xa - xb)

    if ( i < 0 or j < 0 ):
        return 0.0

    elif ( i == 0 and j == 0 ):
        return 1.0

    elif ( i == 0 ):
        return xpb * od_overlap(i, j - 1, alpha, beta, xa, xb) + 0.5 / eta * (j - 1) * od_overlap(i, j-2, alpha, beta, xa, xb)
    
    else:
        return xpa * od_overlap(i-1, j, alpa, beta, xa, xb) + 0.5 / eta * ((i - 1) * od_overlap(i-2, j, alpha, beta, xa, xb) + j*od_overlap(i-1, j-1, alpha, beta, xa, xb))

def compute_gto_overlap(g1, g2):
    exp1 = g1.exp
    exp2 = g2.exp
    gamma = exp1 + exp2

    P = (exp1 * g1.origin + exp2 * g2.origin) / gamma
    PA = P - g1.origin
    PB = P - g2.origin

    fx_arr, fy_arr, fz_arr = compute_f_arr(g1, g2)
    
    sumx = fx_arr[0]
    for i in range( (g1.shell[0] + g2.shell[0]) / 2 + 1 ):
        sumx += fx_arr[2*i] * fact2(2*i-1) / (2*gamma)**i

    sumy = fy_arr[0]
    for j in range( (g1.shell[1] + g2.shell[1]) / 2 + 1 ):
        sumy += fy_arr[2*j] * fact2(2*j-1) / (2*gamma)**j

    sumz = fz_arr[0]
    for k in range( (g1.shell[2] + g2.shell[2]) / 2 + 1 ):
        sumz += fz_arr[2*k] * fact2(2*k-1) / (2*gamma)**k

    rsq = np.linalg.norm(g1.origin - g2.origin)**2 # distance between two orbitals 

    return sumx * sumy * sumz * np.exp(-exp1*exp2*rsq / gamma) * (np.pi / gamma)**1.5

def compute_f_arr(g1, g2):
    fx_arr, fy_arr, fz_arr = [], [], []
    l1, m1, n1 = g1.shell
    l2, m2, n2 = g2.shell

    for i in range(l1 + l2 + 1):
        fx_arr.append( f(i, g1.shell[0], g2.shell[0], g1.origin[0], g2.origin[0]) )
    
    for i in range(m1 + m2 + 1):
        fy_arr.append( f(i, g1.shell[1], g2.shell[1], g1.origin[1], g2.origin[1]) )

    for i in range(n1 + n2 + 1):
        fz_arr.append( f(i, g1.shell[2], g2.shell[2], g1.origin[2], g2.origin[2]) )

    return fx_arr, fy_arr, fz_arr

def f(j, l, m, a, b):
    res = 0.0

    Ninterm = l + m - j

    for k in range(Ninterm + 1):
        if ( Ninterm - k > m ):
            continue
        if ( k > l ):
            break

        res += sp.binom(l, k) * sp.binom(m, Ninterm-k) * a**k * b**(Ninterm - k)

    return res

def compute_sto_overlap( s1, s2 ):
    norm = 0.0

    for g1 in s1.orbs:
        for g2 in s2.orbs:
            norm += g1.exp * g2.exp * compute_gto_overlap(g1, g2)

    return norm


origin = [0.0, 0.0, 0.0]
shell = (0, 0, 0)
Hexps   = [3.42525091, 0.62391373, 0.16885540] 
Hcoeffs  = [0.15432897, 0.53532814, 0.44463454]

s = STO(orbs = [GTO(origin, shell, Hexps[0], Hcoeffs[0]),
                GTO(origin, shell, Hexps[1], Hcoeffs[1]),
                GTO(origin, shell, Hexps[2], Hcoeffs[2])])

print(compute_sto_overlap(s, s))
