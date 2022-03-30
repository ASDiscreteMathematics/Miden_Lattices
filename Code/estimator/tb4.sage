from estimator import *
import numpy

sec_lvls = [128, 160, 192, 256]
dimensions = [3*256, 3*256, 4*256, 6*256]
noise_bounds = [4, 5, 6, 8]

for i in range(4):
    q = 2^16
    sec_lvl = sec_lvls[i]
    dimension = dimensions[i]
    noise_bound = noise_bounds[i]
    tb4 = LWE.Parameters(n=dimension, q=q, Xs=ND.Uniform(-noise_bound, noise_bound), Xe=ND.Uniform(-noise_bound, noise_bound))
    print(tb4)
    LWE.estimate.rough(tb4)
    sigma = float(numpy.std(list(range(-noise_bound, noise_bound+1))))
    print("sigma:", sigma)

    max_noise = 2*noise_bound^2 + 2*noise_bound
    print("decryption failure possible?", max_noise >= q/4)

    #q/4 <= noise
    #q <= noise * 4
    #noise == 2*n*noise_bound^2 + (n+1)*noise_bound
    #      == (2*noise_bound^2+noise_bound)*n + noise_bound
    #q - noise_bound <= (2*noise_bound^2+noise_bound)*n
    #(q - noise_bound) / (2*noise_bound^2+noise_bound) <= n
    n = ceil((q - noise_bound) / (2*noise_bound^2+noise_bound))
    print("homomorphic additions before possible failure:", n)
    print("sane?", q/4 <= 2*n*noise_bound^2 + (n+1)*noise_bound)

    print("")
    # LWE.estimate(tb)


