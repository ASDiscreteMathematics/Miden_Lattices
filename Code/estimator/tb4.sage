# LWEParameters(n=768, q=65536, Xs=D(σ=2.00), Xe=D(σ=2.00), m=+Infinity, tag=None)
# usvp                 :: rop: ≈2^148.9, red: ≈2^148.9, δ: 1.003356, β: 510, d: 1556, tag: usvp
# dual_hybrid          :: rop: ≈2^152.4, mem: ≈2^148.8, m: 839, β: 521, d: 1597, ↻: 1, ζ: 10, tag: dual_hybrid
# arora-gb             :: rop: ≈2^229.8, dreg: 17, mem: ≈2^229.8, t: 8, m: ≈2^114.9
# error probability less than 2^ -313.167109985835258324752741633145349068546678022645984924962722614594336564000176550675011470910
# 
# LWEParameters(n=1024, q=65536, Xs=D(σ=2.00), Xe=D(σ=2.00), m=+Infinity, tag=None)
# usvp                 :: rop: ≈2^211.7, red: ≈2^211.7, δ: 1.002599, β: 725, d: 2060, tag: usvp
# dual_hybrid          :: rop: ≈2^215.5, mem: ≈2^206.2, m: 1081, β: 738, d: 2092, ↻: 1, ζ: 13, tag: dual_hybrid
# arora-gb             :: rop: ≈2^243.8, dreg: 17, mem: ≈2^243.8, t: 8, m: ≈2^121.9
# error probability less than 2^ -317.415037499278843818546261056052183491240185592307518939544247345458901772205641437477719525082
# 
# LWEParameters(n=1536, q=65536, Xs=D(σ=2.00), Xe=D(σ=2.00), m=+Infinity, tag=None)
# usvp                 :: rop: ≈2^342.5, red: ≈2^342.5, δ: 1.001809, β: 1173, d: 3030, tag: usvp
# dual_hybrid          :: rop: ≈2^346.6, mem: ≈2^343.2, m: 1547, β: 1186, d: 3059, ↻: 1, ζ: 24, tag: dual_hybrid
# arora-gb             :: rop: ≈2^263.5, dreg: 17, mem: ≈2^263.5, t: 8, m: ≈2^131.7
# error probability less than 2^ -312.320519900494553958963709262987701123441500221445234388494985490636651285696057791506703357685

from estimator import *
import numpy

load("noise_distribution.sage")

sec_lvls = [128, 192, 256]
dimensions = [3*256,  4*256, 6*256]
#noise_bounds = [4, 6, 8]
#sigmas = [2, 2, 2]

for i in range(3):
    precision = 320
    #cutoff_point = 2^8-1 # faster
    cutoff_point = 2^14-1
    RR = Reals(precision)
    q = 2^16
    sec_lvl = sec_lvls[i]
    dimension = dimensions[i]
    #noise_bound = noise_bounds[i]
    #sigma = sigmas[i]

    #tb4 = LWE.Parameters(n=dimension, q=q, Xs=ND.DiscreteGaussian(sigma), Xe=ND.DiscreteGaussian(sigma))
    tb4 = LWE.Parameters(n=dimension, q=q, Xs=ND.CenteredBinomial(8), Xe=ND.CenteredBinomial(8))
    print(tb4)
    #print("estimating ...")
    LWE.estimate.rough(tb4)

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

    # 3. Performing Diffie-Hellman generates a noise term of the form  a*d - c*b + b - d.
    #diffie_hellman_var = outer_ring_var * 2 + 2 * origin_var
    noise_diffie_hellman = noise_outer_ring * 2 + noise_term_origin * 2
    sigma_diffie_hellman = sqrt(RR(2)*sigma_outer_ring^2 + RR(2)*RR(sigma_origin)^2)
    #print("dh sigma:", noise_diffie_hellman.std())
    #print("expected:", sigma_diffie_hellman)

    # ## Compute probability of flipping a bit in one coefficient
    # with enough precision

    # The bit in one coefficient flips if the noise lies outside [-q/4; q/4]
    # We take the left side of the CDF twice.
    #cutoff_point = RR(-q/4 + 0.5)
    #cdf = RR(0.5) * (RR(1) + erf(cutoff_point / (sqrt(RR(2)) * diffie_hellman_sigma)))
    #error_probability = RR(2) * cdf
    #print("error probability approx 2^", log(error_probability, RR(2)))
    print("error probability less than 2^", log(noise_diffie_hellman.overflow(), noise_diffie_hellman.overflow().parent()(2)))
    # calculate error probability again but using sigma
    #cdf = RR(0.5) * (RR(1) + erf(RR(-2^14-1) / (sqrt(RR(2)) * sigma_diffie_hellman)))
    #error_probability = RR(2) * cdf
    #print("assuming Gaussian, less than 2^", log(error_probability, RR(2)))
    #print("statistical distance:", NoiseDistributionWithOverflow.DiscreteGaussian(cutoff_point, precision, sigma_diffie_hellman).statistical_distance(noise_diffie_hellman))

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


