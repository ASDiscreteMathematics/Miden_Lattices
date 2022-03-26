from estimator import *

tb4 = LWE.Parameters(n=1024, q=previous_prime(2^16), Xs=ND.Uniform(-5,5), Xe=ND.Uniform(-5,5))

LWE.estimate.rough(tb4) # does not crash

LWE.estimate(tb4) # does crash


