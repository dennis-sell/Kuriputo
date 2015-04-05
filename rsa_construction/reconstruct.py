#!/usr/bin/sage

from sage.all import *
from Crypto.PublicKey import RSA

rsa_key = RSA.generate(1024)

def recover_from_ds(key):
    d_p = mod(key.d, key.p - 1)
    d_q = mod(key.d, key.q - 1)

    a = ZZ.random_element(0, d_p)

    p = gcd(pow(a, ZZ(key.e) * ZZ(d_p) - 1, key.n) - 1, key.n)
    q = gcd(pow(a, ZZ(key.e) * ZZ(d_q) - 1, key.n) - 1, key.n)
    return p, q

