#!/usr/bin/sage

from sage.all import *
import struct
import re
from Crypto.PublicKey import RSA
from Crypto.Cipher import AES
from Crypto import Random
from operator import mul

p = 131
q = 41
q1 = 37
r1 = 5
X1 = 8

N = p * q
a1 = q1 * p + r1

def satisfies_conditions(N, beta, bounds):
    m = float(len(bounds))
    geo_average_noise = reduce(mul, bounds, 1) ** (1.0 / m)
    beta_bound = beta > 1 / sqrt(log(N, 2))
    noise_bound = geo_average_noise < pow(N, pow(beta, m+1/m))
    return bool(beta_bound) and bool(noise_bound)


def simple_howgrave_graham(N, a1, X1):
    L = matrix(ZZ, 2, 2,
                  [[ -1*X1, a1],
                   [    0,  N]])
    B = L.LLL()

    x = PolynomialRing(RationalField(), 'x').gen()
    f = sum(B[0][1-i] * x**i / float(X1)**i for i in range(2))
    m = f.roots()
    # How to get to a value of r?
    return m[0][0]


def make_hg_rows(N, a, x, k, t):
    L = []
    for i in range(t):
        current_row = []

        for j in range(t):
            coefficient = binomial(i, j) * pow(-1 * a1, i-j) * pow(X1, j)
            coefficient *= pow(N, k-i) if i <= k else 1 # might be < k
            current_row.insert(0, coefficient)
        L.append(current_row)
    return L


def howgrave_graham_formula(k, t):
    var("a,X,N")
    L = make_hg_rows(N, a, X, k, t)
    L = matrix(t, t, L)
    return L


def howgrave_graham(N, a1, X1, k, t):
    print "Howgrave graham for k=%i, t=%i" % (k,t)
    L = make_hg_rows(N, a1, X1, k, t)
    L = matrix(ZZ, t, t, L)
    print L

    B = L.LLL()
    print B
    x = PolynomialRing(RationalField(), 'x').gen()
    f = sum(B[0][1-i] * x**i / float(X1)**i for i in range(2))
    return f.roots()


print satisfies_conditions(N, log(p, N), [X1])
r = simple_howgrave_graham(N, a1, X1)


#print howgrave_graham_formula(2, 3)
print howgrave_graham(N, a1, X1, 2, 3)

