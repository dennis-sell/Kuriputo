#!/usr/bin/sage

from sage.all import *
from Crypto.PublicKey import RSA
from Crypto.Cipher import AES
from Crypto import Random

import itertools
from operator import mul

class ACD_solver:
    def __init__(self, m, lenn, lenp, lenr):
        self.m = m
        self.lenn = lenn
        self.lenp = lenp
        self.lenr = lenr
        self.X = 2 ** lenr
        # Should these be prime numbers?
        self.p = ZZ.random_element(2**(lenp-1), 2**lenp)
        self.q = ZZ.random_element(2**(lenn-lenp-1), 2**(lenn-lenp))
        self.N = self.p*self.q
        self.r_list = [ZZ.random_element(-self.X, self.X) for i in range(m)]

        q_list = [ZZ.random_element(2**(lenn-lenp)) for i in range(m)]
        self.a_list = [self.p*q + r for q, r in zip(q_list, self.r_list)]
        self.R = PolynomialRing(QQ,m,'x',order='lex')


    def gen_lattice(self, t, k):
        variables = self.R.gens()
        f_list = [a - self.X * x for a, x in zip(self.a_list, variables)]

        indices = [list(ind) + [max(0, k - sum(ind))]
                    for ind in itertools.product(range(t+1))
                    if sum(ind) <= t] # sum i_j <= t

        dim = len(indices)
        functions = f_list + [self.N]
        pindex = [prod(map(operator.pow, functions, exponents)) for exponents in indices]
        monomial_list = [prod(map(operator.pow, variables, exponents[:-1]))
                            for exponents in indices]


        L = matrix(ZZ, dim, dim)
        for i in range(dim):
            for j in range(dim):
                L[i, j] = pindex[i].monomial_coefficient(monomial_list[j])

        def getf(M,i):
            return sum(self.R(b/self.X**monomial.degree()) * monomial
                        for b, monomial in zip(M[i],monomial_list))

        return L, getf


    def approx_factor(self, B, t, k):
        dim = B.nrows()
        det = (self.X**(self.m * dim * t / (self.m+1))
              * self.n**(binomial(k+self.m,self.m) * k/(self.m+1)))
        v_m_norm = norm(B[self.m-1], 2)
        # Not sure I'll get back to this later
        return RR((v_m_norm / (det ** (1/dim))) ** (1/dim))


    def solve(self, t, k): # Use magma later?
        print "Generating lattice"
        A,getf = self.gen_lattice(t,k)
        print "Running LLL",
        start = time.clock()
        B = A.LLL()
        print time.clock()-start,
        return B, getf


    def groebner(self, B, getf, basis_size=0):
        if not basis_size:
            basis_size = B.ncols()-1
            """ Are we sure??? How about
            basis_size = B.n
            """
        R = self.R
        algorithm = 'libsingular:groebner'
        #if use_magma:
        #    rp = random_prime(2**(self.lenr+2),lbound=2**(self.lenr+1))
        #    R = PolynomialRing(GF(rp),self.m,'x',order='lex')
        #    algorithm = 'magma:GroebnerBasis'
        I = (tuple(getf(B,i) for i in range(basis_size)))*R
        #print "groebner basis:",
        start = time.clock()
        J = I.groebner_basis(algorithm)
        print time.clock()-start

        root_count = 0
        for b in J:
            if len(b.variables()) == 1:
                factorization = b.factor(proof=False)
                for poly,j in factorization:
                    if poly.degree() == 1:
                        root1 = -ZZ(poly.monomial_coefficient(poly.variables()[0]**0))
                        if use_magma:
                            rp = R.characteristic()
                            if root1 < -2**self.lenr:
                                root1 = root1 + rp
                            else:
                                if root1 > 2**self.lenr:
                                    root1 = root1 - rp
                        if root1 in self.r_list:
                            root_count += 1
#        print "Found", root_count, "correct roots."


def satisfies_conditions(N, beta, bounds):
    m = len(bounds)
    geo_average_noise = prod(bounds) ** (1.0 / m)

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


p = 10313213783923
q =  4793213432437
q1 = 3732132143239
r1 = -1323

X1 = 2000
N = p * q
a1 = q1 * p + r1

#print satisfies_conditions(N, log(p, N), [X1])
fs, Bs = simple_howgrave_graham(N, a1, X1)


#print howgrave_graham_formula(2, 3)
f, B = howgrave_graham(N, a1, X1, 2, 3)

