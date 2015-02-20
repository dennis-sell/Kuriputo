#!/usr/bin/sage

from sage.all import *
import struct
import re
from Crypto.PublicKey import RSA
from Crypto.Cipher import AES
from Crypto import Random
from operator import mul

def satisfies_conditions(N, beta, bounds):
    m = len(bounds)
    geo_average_noise = reduce(mul, bounds, 1) ** (1.0 / m)

    beta_bound = bool(beta > 1 / log(N, 2) ** (1/2))
    noise_bound = bool(geo_average_noise < pow(N, pow(beta, m+1/m)))
    return beta_bound and noise_bound


def simple_howgrave_graham(N, a1, X1):
    dim = 2
    L = matrix(ZZ, dim, dim,
                  [[-1*a1, X1],
                   [    N,  0]] )
    B = L.LLL()

    x = PolynomialRing(RationalField(), 'x').gen()
    f = sum(B[0][i] * x**i / X1**i for i in range(dim))
    return f, B


def make_hg_rows(N, a, X, k, t, matrix_form=False):
    L = []
    for i in range(t+1):
        current_row = []

        for j in range(t+1):
            coefficient = binomial(i, j) * pow(-1 * a, i-j) * pow(X, j)
            coefficient *= pow(N, k-i) if i <= k else 1 # might be < k
            current_row.append(coefficient)
        L.append(current_row)

    if matrix_form:
        return matrix(L)
    return L

"""
To get formula

R.<N,a,x> = PolynomialRing(QQ, 3)
matrix_formula = ACD.make_hg_rows(N,a,x, <k>, <t>, matrix_form=True)
"""

def howgrave_graham(N, a1, X1, k, t):
    print "Howgrave graham for k=%i, t=%i" % (k,t)
    rows = make_hg_rows(N, a1, X1, k, t)
    L = matrix(ZZ, t+1, t+1, rows)
    B = L.LLL()

    x = PolynomialRing(RationalField(), 'x').gen()
    f = sum(B[0][i] * (x**i) / (X1**i) for i in range(t+1))
    return f, B


p = 10313213783924
q =  4793213432432
q1 = 3732132143232
r1 = -1323

X1 = 2000
N = p * q
a1 = q1 * p + r1

#print satisfies_conditions(N, log(p, N), [X1])
fs, Bs = simple_howgrave_graham(N, a1, X1)


#print howgrave_graham_formula(2, 3)
f, B = howgrave_graham(N, a1, X1, 2, 3)

