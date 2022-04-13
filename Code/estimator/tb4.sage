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

    # 1. Multiplying two Fp field elements generates 8 noise terms on one coefficient.
    # To see this, represent two field elements as polynomials in the ring Z[X]/<X^4-X^2+1>.
    # Then (a + b*X + c*X^2 + d*X^3) * (e + f*X + g*X^2 + h*X^3) =
    #   (a*e - b*h - c*g - d*f - d*h) +
    #   (a*f + b*e - c*h - d*g)*X +
    #   (a*g + b*f + c*e - d*h + b*h + c*g + d*f + d*h)*X^2 +  <-- this is the coefficient with 8 noise terns
    #   (d*e + c*f + b*g + a*h + c*h + d*g)*X^3.
    tb_noise_factor = 8
    
    # 2. Multiplying two ring elements in Fp[Y]/<Y^d+1> generates d noise terms on all coefficients.

    # 3. Performing noise Diffie-Hellman generates a noise term of the form  a*d - c*b + b - d.

    max_noise = 2 * tb_noise_factor * dimension * noise_bound^2 + 2 * noise_bound
    print("decryption failure possible?", max_noise >= q/4)
    print("max noise:", max_noise, "versus q/4:", floor(q/4))

    # for decryption failure:
    # q/4 <= noise
    # q <= noise * 4
    # noise == n*2*tbnf*noise_bound^2 + (n+1)*noise_bound
    #       == (2*tbnf*noise_bound^2+noise_bound)*n + noise_bound
    # q - noise_bound <= (2*tbnf*noise_bound^2+noise_bound)*n
    # (q - noise_bound) / (2*tbnf*noise_bound^2+noise_bound) <= n
    n = ceil((q - noise_bound) / (2*tb_noise_factor*noise_bound^2+noise_bound))
    print("homomorphic additions before possible failure:", n)
    print("sane?", q/4 > n * tb_noise_factor * dimension * 2 * noise_bound^2 + (n+1) * noise_bound)

    print("")
    # LWE.estimate(tb)


