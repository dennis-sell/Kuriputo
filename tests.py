import csv
import sys
import time

from ACD import *

def printe(x):
    sys.stderr.write(str(x) + "\n")


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


test_params = [
    (1,1000,200,36,41,8),
    (1,1000,400,154,40,16),
    (1,1000,400,156,82,33),

    (2,1000,200,72,9,4),
    (2,1000,400,232,10,6),
    (2,1999,400,238,15,9),

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

    (7,1000,200,120,3,2),
    (7,1000,400,311,3,2),

    (12,1000,400,347,1,1),
    (18,1000,400,364,1,1),
    (24,1000,400,372,1,1),
    (48,1000,400,383,1,1),
    (96,1000,400,387,1,1)]


def run_tests_2(filename, test_params=test_params, tk_limit=10):
    # Find minimal k,t at which this will work
    # Run that lattice
    # Compare it to higher values.
    # Question: Is it optimal for the minimal t for the minimal k for which it will work?

    start_time = time.time()
    fieldnames = ["m","logn","logp","logr","t","k","dim",
                    "generating_time","LLL_time","gtime",
                    "success","use_rounding"]
    csv_file = open(filename, "wb")
    writer = csv.writer(csv_file)
    for m,logn,logp,logr,_,_ in test_params:
        hulk = ACD_solver(m,logn,logp,logr)

        tks = hulk.find_all_tk(rangelim=50)
        tks = sorted(tks, key=lambda (t,k,d): d)[:tk_limit]
        printe(tks)

        for t, k, dim in tks:
            for rounding_const in [None, 2, 8]:
                B,getf,(generating_time, LLL_time) = \
                        hulk.solve(t, k,
                                   use_magma=True,
                                   return_times=True,
                                   rounding_const=rounding_const)
                roots, gtime = hulk.groebner(B, getf, 2, use_magma=True, return_time=True)
                success = hulk.correct_roots(roots)

                writer.writerow((m, logn, logp, logr, t, k, dim,
                                 generating_time, LLL_time, gtime,
                                 success, rounding_const))
    print "time:", time.time() - start_time
    csv_file.close()


def CLT_toy():
    params = (10,290000,988,26)
    acd = ACD_solver(*params, verbose=True)
    acd.find_roots()

def run_tests():
    time_results = []
    for m,logn,logp,logr,t,k in test_params:
        print m, logn, logp, logr
        hulk = ACD_solver(m,logn,logp,logr)
        print t, k, binomial(t+m,m)
        B,getf,(generating_time, LLL_time) = hulk.solve(t, k,
                                                    use_magma=True,
                                                    return_times=True)

        _, gtime = hulk.groebner(B,getf, 0, use_magma=True, return_time=True)
        print generating_time, LLL_time, gtime
        print
        time_results.append((generating_time, LLL_time, gtime))
    return time_results

if __name__=="__main__":
    print run_tests()
    #print CLT_toy()
