from ACD import *

u = univariate_acd_solver(100, 65, 10)
r = u.howgrave_graham(3, 2)



def test_tk(k, t, lenn, lenp, lenr, trials):
    successes = 0
    for i in range(trials):
        #acd = acd_solver(m, lenn, lenp, lenr)
        #if acd.find_roots(k, t) == acd.r_list:
        #    successes += 1
        u = univariate_acd_solver(lenn, lenp, lenr)
        if u.howgrave_graham(k, t) == u.r:
            successes += 1
    success_ratio = float(successes) / trials
    print "success rate is %i out of %i" % (successes, trials)
    print success_ratio
    return success_ratio


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

    time_results = []
    for m,logn,logp,logr,t,k in tests:
        print m, ", ", logn, ", ", logp, ", ", logr
        hulk = ACD_solver(m,logn,logp,logr)
        print t, ", ", k, ", ", binomial(t+m,m)
        start = time()
        B,getf,(generating_time, LLL_time) = hulk.solve(t, k,
                                                    use_magma=True,
                                                    return_times=True)

        _, gtime = hulk.groebner(B,getf, 0, use_magma=True, return_time=True)
        time_results.append((generating_time, LLL_time, gtime))
    return time_results
