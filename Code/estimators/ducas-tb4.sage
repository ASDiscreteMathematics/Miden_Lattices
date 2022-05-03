# ducas-tb4.sage
# Estimate security of the TB4 variants using Ducas et al.'s estimator.
# Usage:
#  - clone the leaky estimator using e.g.
#  $> git clone https://github.com/lducas/leaky-LWE-Estimator.git
#  - create a directory in leaky-LWE-Estimator/ e.g. tb4/ using e.g.
#  $> mkdir tb4
#  - put this file in the director leaky-LWE-Estimator/tb4/
#  - run:
#  $> sage ducas-tb4.sage

def main():

    # core SVP cost models
    cost_model_c = lambda beta: 0.292 * beta
    cost_model_q = lambda beta: 0.265 * beta

    # TB variant parameter sets
    sec_lvls = [128, 192, 256]
    dimensions = [3*256,  4*256, 6*256]
    cutoff_point = 2^14-1

    # for all variants
    for i in range(3):
        sec_lvl = sec_lvls[i]
        dimension = dimensions[i]
        module_dimension = dimension // 256

        print("TB4/", module_dimension, " sec_lvl =", sec_lvl)
 
        load("../framework/instance_gen.sage")

        # estimate security using leaky estimator
        load("../framework/proba_utils.sage")
        D = build_centered_binomial_law(8)

        n = dimension
        q = 2^16
        m = dimension
        A, b, dbdd = initialize_from_LWE_instance(DBDD_predict, n, q, m, D, D)
        _ = dbdd.integrate_q_vectors(q, report_every=20)
        beta, delta = dbdd.estimate_attack()
       
        print("beta:", beta)
        print("delta:", delta) 
        print("classical:", cost_model_c(beta))
        print("quantum:", cost_model_q(beta))

        print("---")
        print("")

if __name__ == '__main__':
    main()
