#!/usr/bin/sage

from sage.all import *
from Crypto.PublicKey import RSA
from Crypto.Cipher import AES
from Crypto import Random

import itertools
import time


def gen_index(t,k,m):
    if t==1:
        return [tuple([0]*(m-i)+[1]+[0]*i) for i in range(m+1)]
    else:
        # Must optimize this shit
        """
        rows = [tuple([0]*(m+1))]
        last_row = list(rows)
        new_rows = []
        # At each step, we generate all indices with a sum which is one higher.
        for row_sum in range(t):
            for index in last_row:
                for x in range(m+1):
                    # Increments index at position x by 1
                    inc_index = index[0:x-1] + (index[x]+1,) + index[x+1:m+1]
                    rows.append(inc_index)
                    new_rows.append(inc_index)
            last_rows = new_rows
            new_rows = []

        return rows
        """
        rows = []
        for indices in itertools.product(range(t+1),repeat=m):
            if sum(indices) < t+1:
                rows.append(tuple(list(indices)+[max(0,k-sum(indices))]))
        return rows


class ACD_solver:
    def __init__(self, m, lenn, lenp, lenr, verbose=False):
        self.verbose = verbose
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


    def find_roots(self):
        t, k, dim = self.find_tk(rangelim=40, dimlim=2000)
        if not t and not k:
            if self.verbose:
                print "Not solvable or dimension too large"
            return []
        if self.verbose:
            print "Lattice (t, k, dim) =", t, k, dim
        B, getf = self.solve(t, k)
        roots = self.groebner(B, getf)
        return roots


    def solve(self, t, k, use_magma=False, return_times=False):
        if self.verbose:
            print "Generating lattice",
        start = time.time()
        A,getf = self.gen_lattice(t,k)
        if use_magma:
            A = magma(A)
        generating_time = time.time()-start
        if self.verbose:
            print generating_time

        if self.verbose:
            print "Running LLL",
        start = time.time()
        B = A.LLL()
        if use_magma:
            B = B.sage()
        LLL_time = time.time()-start
        if self.verbose:
            print LLL_time

        self.check(B, getf)
        if return_times:
            return B, getf, (generating_time, LLL_time)
        else:
            return B, getf


    def gen_lattice(self, t, k):
        variables = self.R.gens()
        f_list = [a - self.X * x for a, x in zip(self.a_list, variables)]

        indices = gen_index(t,k,self.m)
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


    def find_tk(self, rangelim=20, dimlim = 200, lllfactor=log(1.01)/log(2), lenr=0):
        t,k=0,0
        if lenr == 0:
            lenr = self.lenr
        for test_k, test_t in itertools.product(range(1, rangelim), repeat=2):
            if test_k > test_t:
                continue
            dim = binomial(test_t+self.m,self.m)
            if dim > dimlim:
                continue
            if self.check_tk(test_t,test_k,lllfactor):
                k = test_k
                t = test_t
                dimlim = dim
        return t,k,dimlim


    def check_tk(self, test_t, test_k, lllfactor=log(1.01)/log(2)):
        dim = binomial(test_t+self.m,self.m)
        veclen = (log(dim)/(2*log(2)) +
                  dim*lllfactor +
                  (self.lenr*dim*test_t*self.m/(self.m+1)
                    + self.lenn*binomial(test_k+self.m,self.m)*test_k/(self.m+1))
                        /dim)
        if self.verbose:
            print dim, float(veclen), self.lenp*test_k, bool(veclen < self.lenp*test_k)
        return bool(veclen < self.lenp*test_k)


    """ Checks that each r is a root of the polynomial """
    def check(self, B, getf):
        if any([apply(getf(B, i), self.r_list) for i in range(self.m)]):
            if self.verbose:
                print "Failed: Polynomials do not vanish at roots."
            return False
        return True


    def groebner(self, B, getf, basis_size=0, use_magma=False, return_time=False):
        if not basis_size:
            basis_size = B.ncols()-1
            """ Are we sure??? How about
            basis_size = B.n
            """
        R = self.R
        algorithm = 'libsingular:groebner'
        if use_magma:
           rp = random_prime(2**(self.lenr+2),lbound=2**(self.lenr+1))
           R = PolynomialRing(GF(rp),self.m,'x',order='lex')
           algorithm = 'magma:GroebnerBasis'
        I = (tuple(getf(B,i) for i in range(basis_size)))*R
        if self.verbose:
            print "groebner basis:",
        start = time.time()
        J = I.groebner_basis(algorithm)
        groebner_time = time.time()-start
        if self.verbose:
            print groebner_time

        roots = []
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
                        roots.append(root1)
        extended_groebner_time = time.time() - start
        if self.verbose:
            print "Found", len(roots), "correct roots."
            print "Same as original r's", set(roots) == set(self.r_list)
        if return_time:
            return roots, groebner_time
        else:
            return roots



class univariate_acd_solver:

    def __init__(self, lenn, lenp, lenr):
        self.lenn = lenn
        self.lenp = lenp
        self.lenr = lenr

        lenq = lenn - lenp

        self.p = ZZ.random_element(2**(lenp-1), 2**lenp)
        self.q = ZZ.random_element(2**(lenq-1), 2**lenq)
        self.q1 = ZZ.random_element(2**(lenq-1), 2**lenq)
        self.r = ZZ.random_element(-2**lenr, 2**lenr)

        self.X = 2**lenr
        self.a1 = self.p * self.q1 + self.r
        self.N = self.p * self.q

    def satisfies_conditions(self, N, beta, bounds):
        m = len(bounds)
        geo_average_noise = prod(bounds) ** (1.0 / m)

        beta_bound = bool(beta > 1 / log(N, 2) ** (1/2))
        noise_bound = bool(geo_average_noise < pow(N, pow(beta, m+1/m)))
        return beta_bound and noise_bound


    def simple_howgrave_graham(self):
        a1 = self.a1; X = self.X; N = self.N;
        dim = 2
        L = matrix(ZZ, dim, dim,
                      [[-1*a1, X],
                       [    N, 0]] )
        B = L.LLL()

        x = PolynomialRing(RationalField(), 'x').gen()
        f = sum(B[0][i] * x**i / X**i for i in range(dim))
        return f, B


    def make_hg_rows(self, t, k, matrix_form=False):
        a1 = self.a1; X = self.X; N = self.N;
        L = []
        for i in range(t+1):
            current_row = []

            for j in range(t+1):
                coefficient = binomial(i, j) * pow(-1 * a1, i-j) * pow(X, j)
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

    def howgrave_graham(self, t, k):
        X = self.X
        rows = self.make_hg_rows(t, k)
        L = matrix(ZZ, t+1, t+1, rows)
        B = L.LLL()

        x = PolynomialRing(RationalField(), 'x').gen()
        f = sum(B[0][i] * (x**i) / (X**i) for i in range(t+1))
        roots = f.roots()
        return roots[0][0] if roots else None
