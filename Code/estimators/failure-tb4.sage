# Bound the probability of observing a decryption or decapsulation failure.
# Usage:
# $> sage failure-tb4.sage

load("noise_distribution.sage")

sec_lvls = [128, 192, 256]
dimensions = [3*256,  4*256, 6*256]

for i in range(3):
    precision = 1000
    cutoff_point = 2^14-1
    RR = Reals(precision)
    q = 2^16
    sec_lvl = sec_lvls[i]
    dimension = dimensions[i]

    print("TB4  q = 2^16  Xs = Xe = CenteredBinomial(8)  sec lvl =", sec_lvl, " dim =", dimension)
    print("computing failure probability in a single (inner) coefficient ...")

    # ## Compute statistics of noise terms

    noise_term_origin = NoiseDistributionWithOverflow.CenteredBinomial(cutoff_point, precision, ncrumbs=8)
    sigma_origin = 2
    #print("origin sigma:", noise_term_origin.std())
    #print("expected:", sigma_origin)

    # 0. Some noise terms are actually the product of two original noise terms
    noise_term_product = noise_term_origin * noise_term_origin
    sigma_product = sigma_origin * sigma_origin
    #print("product sigma:", noise_term_product.std())
    #print("expected:", sigma_product)


    # 1. Multiplying two Fp field elements generates 8 square noise terms on one coefficient.
    # To see this, represent two field elements as polynomials in the ring Z[X]/<X^4-X^2+1>.
    # Then (a + b*X + c*X^2 + d*X^3) * (e + f*X + g*X^2 + h*X^3) =
    #   (a*e - b*h - c*g - d*f - d*h) +
    #   (a*f + b*e - c*h - d*g)*X +
    #   (a*g + b*f + c*e - d*h + b*h + c*g + d*f + d*h)*X^2 +  <-- this is the coefficient with 8 square noise terms
    #   (d*e + c*f + b*g + a*h + c*h + d*g)*X^3.
    tb4_noise_factor = 8
    noise_term_largest = noise_term_product * tb4_noise_factor
    sigma_largest = sigma_product * sqrt(RR(tb4_noise_factor))
    #print("largest sigma:", noise_term_largest.std())
    #print("expected:", sigma_largest)
    
    # 2. Multiplying two ring elements in Fp[Y]/<Y^d+1> generates d noise terms on all coefficients.
    # We model them as all independent
    outer_ring_noise_factor = dimension//4
    #outer_ring_var = outer_ring_noise_factor * largest_term_variance
    noise_outer_ring = noise_term_largest * outer_ring_noise_factor
    sigma_outer_ring = sigma_largest * sqrt(RR(outer_ring_noise_factor))
    #print("outer sigma:", noise_outer_ring.std())
    #print("expected:", sigma_outer_ring)

	# 3. We ue the above ring in a module, giving rise to a inner product
    module_dimension = dimension//256
    noise_module = noise_outer_ring * module_dimension
    sigma_module = sigma_outer_ring * sqrt(RR(module_dimension))

    # 4. Performing Diffie-Hellman generates a noise term of the form  a*d - c*b + b - d.
    #diffie_hellman_var = outer_ring_var * 2 + 2 * origin_var
    noise_diffie_hellman = noise_module * 2 + noise_term_origin * 2
    sigma_diffie_hellman = sqrt(RR(2)*sigma_module^2 + RR(2)*RR(sigma_origin)^2)
    #print("dh sigma:", noise_diffie_hellman.std())
    #print("expected:", sigma_diffie_hellman)

    # sanitize: make all probabilities are positive
    # (the FFT can make some negative, but on the order of precision^(1/3))
    for n, p in noise_diffie_hellman.dictionary.items():
        if p < 0:
            #print("Problem! Have negative probability:", n, "->", p)
            p = RR(0)

    # 5. Compute probability of flipping a bit in one coefficient
    # with enough precision

    # The bit in one coefficient flips if the noise lies outside [-q/4; q/4]
    # recall that cutoff_point is approx RR(-q/4 + 0.5)
    #print("error probability approx 2^", log(error_probability, RR(2)))
    print("error probability less than 2^", log(noise_diffie_hellman.overflow(), RR(2)))
    # calculate error probability again but using sigma
    #cdf = RR(0.5) * (RR(1) + erf(RR(-2^14-1) / (sqrt(RR(2)) * sigma_diffie_hellman)))
    #error_probability = RR(2) * cdf
    #print("assuming Gaussian, less than 2^", log(error_probability, RR(2)))
    #print("statistical distance:", NoiseDistributionWithOverflow.DiscreteGaussian(cutoff_point, precision, sigma_diffie_hellman).statistical_distance(noise_diffie_hellman))

    #max_noise = 2 * tb_noise_factor * dimension * noise_bound^2 + 2 * noise_bound
    #print("decryption failure possible?", max_noise >= q/4)
    #print("max noise:", max_noise, "versus q/4:", floor(q/4))

    print("")


