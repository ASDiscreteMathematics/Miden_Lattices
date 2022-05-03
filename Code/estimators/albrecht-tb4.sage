# albrecht-tb4.sage
# Put this file in the lattice-estimator/ directory, obtained from
# checking out from github:
# $> git checkout https://github.com/malb/lattice-estimator.git
# Then run
# $> sage albrecht-tb4.sage

from estimator import *
import numpy

sec_lvls = [128, 192, 256]
dimensions = [3*256,  4*256, 6*256]

for i in range(3):
    q = 2^16
    sec_lvl = sec_lvls[i]
    dimension = dimensions[i]

    tb4 = LWE.Parameters(n=dimension, q=q, Xs=ND.CenteredBinomial(8), Xe=ND.CenteredBinomial(8))
    print(tb4)
    print("estimating ...")
    LWE.estimate.rough(tb4)
    print("")
    LWE.estimate(tb4) # might crash
    print("")

