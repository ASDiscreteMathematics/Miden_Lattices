from estimator import *
import numpy

sec_lvls = [128, 192, 256]
dimensions = [3*256,  4*256, 6*256]
#noise_bounds = [4, 6, 8]
sigmas = [1.5, 1.5, 1.5]

for i in range(3):
    q = 2^16
    sec_lvl = sec_lvls[i]
    dimension = dimensions[i]
    #noise_bound = noise_bounds[i]
    sigma = sigmas[i]

    tb4 = LWE.Parameters(n=dimension, q=q, Xs=ND.DiscreteGaussian(sigma), Xe=ND.DiscreteGaussian(sigma))
    print(tb4)
    LWE.estimate.rough(tb4)

    print("sigma:", sigma)

    # ## Compute sigma of noise terms

    # 0. Given two noise terms, a and b with std(a)=std(b)=std, the standard deviation of a*b is std^2.
    # The standard deviation of a+b is sqrt(2)*std. Let's use variances instead.
    base_var = sigma^2
    product_var = base_var^2

    # 1. Multiplying two Fp field elements generates 8 noise terms on one coefficient.
    # To see this, represent two field elements as polynomials in the ring Z[X]/<X^4-X^2+1>.
    # Then (a + b*X + c*X^2 + d*X^3) * (e + f*X + g*X^2 + h*X^3) =
    #   (a*e - b*h - c*g - d*f - d*h) +
    #   (a*f + b*e - c*h - d*g)*X +
    #   (a*g + b*f + c*e - d*h + b*h + c*g + d*f + d*h)*X^2 +  <-- this is the coefficient with 8 noise terns
    #   (d*e + c*f + b*g + a*h + c*h + d*g)*X^3.
    tb4_noise_factor = 8
    largest_term_variance = 8 * product_var
    
    # 2. Multiplying two ring elements in Fp[Y]/<Y^d+1> generates d noise terms on all coefficients.
    outer_ring_noise_factor = dimension/4
    outer_ring_var = outer_ring_noise_factor * largest_term_variance

    # 3. Performing Diffie-Hellman generates a noise term of the form  a*d - c*b + b - d.
    diffie_hellman_var = outer_ring_var * 2 + 2 * base_var

    # ## Compute probability of flipping a bit in one coefficient
    # with enough precision
    RR = Reals(3000)
    diffie_hellman_sigma = sqrt(RR(diffie_hellman_var))

    # The bit in one coefficient flips if the noise lies outside [-q/4; q/4]
    # We take the left side of the CDF twice.
    cutoff_point = RR(-q/4 + 0.5)
    cdf = RR(0.5) * (RR(1) + erf(cutoff_point / sqrt(RR(2) * diffie_hellman_sigma)))
    error_probability = RR(2) * cdf
    print("error probability approx 2^", log(error_probability, RR(2)))

    #max_noise = 2 * tb_noise_factor * dimension * noise_bound^2 + 2 * noise_bound
    #print("decryption failure possible?", max_noise >= q/4)
    #print("max noise:", max_noise, "versus q/4:", floor(q/4))

    # for decryption failure:
    # q/4 <= noise
    # q <= noise * 4
    # noise == n*2*tbnf*noise_bound^2 + (n+1)*noise_bound
    #       == (2*tbnf*noise_bound^2+noise_bound)*n + noise_bound
    # q - noise_bound <= (2*tbnf*noise_bound^2+noise_bound)*n
    # (q - noise_bound) / (2*tbnf*noise_bound^2+noise_bound) <= n
    #n = ceil((q - noise_bound) / (2*tb_noise_factor*noise_bound^2+noise_bound))
    #print("homomorphic additions before possible failure:", n)
    #print("sane?", q/4 > n * tb_noise_factor * dimension * 2 * noise_bound^2 + (n+1) * noise_bound)

    print("")
    # LWE.estimate(tb4) # crashes
    # print("")


