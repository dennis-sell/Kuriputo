import itertools
import time

class gcd_solver:
    def __init__(self,m,lenn,lenp,lenr):
        self.m = m
        self.lenn = lenn
        self.lenp = lenp
        self.lenr = lenr
        self.p = ZZ.random_element(2^(lenp-1),2^lenp)
        self.q = ZZ.random_element(2^(lenn-lenp-1),2^(lenn-lenp))
        self.n = self.p*self.q
        self.r_list = [ZZ.random_element(-2^self.lenr,2^self.lenr) for i in range(m)]
        self.n_list = [self.p*ZZ.random_element(2^(lenn-lenp)) + r for r in self.r_list]
        self.R = PolynomialRing(QQ,m,'x',order='lex')

    def gen_lattice(self,t,k):
        variables = self.R.gens()
        X = 2^self.lenr
        f_list = [a-X*x for a,x in zip(self.n_list,variables)]

        def gen_index():
            if t==1:
                for i in range(self.m+1):
                    yield tuple([0]*(self.m-i)+[1]+[0]*i)
            else:    
                for indices in itertools.product(range(t+1),repeat=self.m):
                    if sum(indices) < t+1:
                        yield tuple(list(indices)+[max(0,k-sum(indices))])

        indices = list(gen_index())
        dim = len(indices)
        functions = f_list + [self.n]#+ZZ.random_element(-2^self.lenr,2^self.lenr)]
        pindex = [prod(map(operator.pow,functions,exponents)) for exponents in indices]
        monomial_list = [prod(map(operator.pow,variables,exponents[:-1])) for exponents in indices]

        A = MatrixSpace(IntegerRing(),dim,dim)(0)
        for i in range(dim):
            for j in range(dim):
                A[i,j] = pindex[i].monomial_coefficient(monomial_list[j])

        def getf(M,i):
            return sum(self.R(b/X^monomial.degree())*monomial for b,monomial in zip(M[i],monomial_list))

#        print "dimension:", dim
        return A,getf

    def approx_factor(self,B,t,k,getf):
        def l2(i):
            return sqrt(sum(b^2 for b in B[i]))
        def l1(i):
            return sum(abs(b) for b in B[i])

        X = 2^self.lenr
        dim = B.nrows()
        det = X^(self.m*dim*t/(self.m+1))*self.n^(binomial(k+self.m,self.m)*k/(self.m+1))
        return RR((l2(self.m-1)/(det^(1/(dim))))^(1/dim))

    def check(self,B,t,k,getf):
#        def l2(i):
#            return sqrt(sum(b^2 for b in B[i]))
#        def l1(i):
#            return sum(abs(b) for b in B[i])

#        X = 2^self.lenr
#        dim = B.nrows()
#        det = X^(self.m*dim*t/(self.m+1))*self.n^(binomial(k+self.m,self.m)*k/(self.m+1))
#        print "approximation factor of 1st vector:", RR((l2(0)/(det^(1/(dim))))^(1/dim))
#        print "approximation factor of mth vector:", RR((l2(self.m-1)/(det^(1/(dim))))^(1/dim))
#        print "approximation factor of dim-1th vector:", RR((l2(dim-1)/(det^(1/(dim))))^(1/dim))
        if any([apply(getf(B,i),self.r_list) for i in range(self.m)]):
            print "Failed: Polynomials do not vanish at roots."#, [apply(getf(B,i),self.r_list) for i in range(self.m)]
            return False
        return True

    def groebner(self,B,getf,basis_size=0,use_magma=False):
        if not basis_size:
            basis_size = B.ncols()-1
        R = self.R
        algorithm = 'libsingular:groebner'
        if use_magma:
            rp = random_prime(2^(self.lenr+2),lbound=2^(self.lenr+1))
            R = PolynomialRing(GF(rp),self.m,'x',order='lex')
            algorithm = 'magma:GroebnerBasis'
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
                        root1 = -ZZ(poly.monomial_coefficient(poly.variables()[0]^0))
                        if use_magma:
                            rp = R.characteristic()
                            if root1 < -2^self.lenr:
                                root1 = root1 + rp
                            else:
                                if root1 > 2^self.lenr:
                                    root1 = root1 - rp
                        if root1 in self.r_list:
                            root_count += 1
#        print "Found", root_count, "correct roots."


    def solve(self,t,k,use_magma=False):
        print "Generating lattice"
        A,getf = self.gen_lattice(t,k)
        if use_magma:
            A = magma(A)
        print "Running LLL",
        start = time.clock()
        B = A.LLL()
        print time.clock()-start,
        if use_magma:
            B = B.sage()
        self.check(B,t,k,getf)
        return B,getf

    def find_tk(self,rangelim=20,dimlim = 200,lllfactor=log(1.01)/log(2),lenr=0):
        t,k=0,0
        if lenr == 0:
            lenr = self.lenr
        for test_k,test_t in itertools.product(range(1,rangelim),repeat=2):
            if test_k > test_t:
                continue
            dim = binomial(test_t+self.m,self.m)
            if dim > dimlim:
                continue
            if self.check_tk(test_t,test_k,lllfactor,verbose=False):
                k = test_k
                t = test_t
                dimlim = dim
        return t,k,dimlim

    def check_tk(self,test_t,test_k,lllfactor=log(1.01)/log(2),verbose=True):
        dim = binomial(test_t+self.m,self.m)
        veclen = log(dim)/(2*log(2)) + dim*lllfactor + (self.lenr*dim*test_t*self.m/(self.m+1)+self.lenn*binomial(test_k+self.m,self.m)*test_k/(self.m+1))/(dim)
        if verbose:
            print dim, float(veclen), self.lenp*test_k, bool(veclen < self.lenp*test_k)
        return bool(veclen < self.lenp*test_k)

def check_params(m,lenn,lenp,lenr,t,k,lllfactor):
        dim = binomial(t+m,m)
        veclen = log(dim)/(2*log(2)) + dim*lllfactor + (lenr*dim*t*m/(m+1)+lenn*binomial(k+m,m)*k/(m+1))/(dim)
        print dim, float(veclen), lenp*k, bool(veclen < lenp*k)
        return bool(veclen < lenp*k)    

def run_tests():
    tests = [
        (1,1000,200,36,41,8),
        (1,1000,400,154,40,16),
        (1,1000,400,156,82,33),
        (2,1000,200,72,9,4),
        (2,1000,400,232,10,6),
        (3,1000,200,87,5,3),
        (3,1000,400,255,4,3),
        (3,1000,400,268,7,5),
        (4,1000,200,94,3,2),
        (4,1000,400,279,4,3),
        (5,1000,200,108,3,2),
        (5,1000,200,110,4,3),
        (5,1000,400,278,3,2),             
        (6,1000,200,115,3,2),
        (6,1000,400,297,3,2),    
        (7,1000,200,122,3,2),
        (7,1000,400,311,3,2),  
        (12,1000,400,347,1,1),
        (18,1000,400,364,1,1),
        (24,1000,400,372,1,1),
        (48,1000,400,383,1,1),
        (96,1000,400,387,1,1)]
             
    for m,logn,logp,logr,t,k in tests:
        print m,"&",logn,"&",logp,"&",logr,"&",
        hulk = gcd_solver(m,logn,logp,logr)
        print t,"&",k,"&",binomial(t+m,m),",&",
        B,getf = hulk.solve(t,k,False)
#        B,getf = hulk.solve(t,k,True)
        print "&",
        hulk.groebner(B,getf,0,False)
#        hulk.groebner(B,getf,0,True)
