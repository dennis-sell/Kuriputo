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
r1 = 10
X1 = 12

N = p * q
a1 = q1 * p + r1

def satisfies_conditions(N, beta, bounds):
    m = float(len(bounds))
    geo_average_noise = reduce(mul, bounds, 1) ** (1.0 / m)
    beta_bound = beta > 1 / sqrt(log(N, 2))
    noise_bound = geo_average_noise < pow(N, pow(beta, m+1/m))
    return bool(beta_bound) and bool(noise_bound)


def howgrave_graham(N, beta, a1, X1):
    L = matrix(ZZ, 2, 2,
                  [ -1*X1, a1,
                        0,  N])
    B = L.LLL()
    print B

    x = PolynomialRing(RationalField(), 'x').gen()
    f = sum(B[0][1-i] * x**i / float(X1)**i for i in range(2))
    m = f.roots()
    print m

    # How to get to a value of r?
    return m[0][0]



print satisfies_conditions(N, log(p, N), [X1])
r = howgrave_graham(N, log(p, N), a1, X1)
print r - a1

