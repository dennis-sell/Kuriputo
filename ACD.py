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
r1 = 11

N = p * q
a1 = q1 * p + r1

def satisfies_conditions(N, beta, bounds):
    m = float(len(bounds))
    geo_average_noise = reduce(mul, bounds, 1) ** (1.0 / m)
    beta_bound = beta > 1 / sqrt(log(N, 2))
    noise_bound = geo_average_noise < pow(N, pow(beta, m+1/m))
    return bool(beta_bound and noise_bound)

print satisfies_conditions(N, log(p, N), [r1 + 1])

