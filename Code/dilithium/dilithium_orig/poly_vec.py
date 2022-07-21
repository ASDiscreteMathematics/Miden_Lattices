from params import *
from poly import *

class polyvecl:
    def __init__(self):
        self.vec = [poly() for i in range(L)]
        
class polyveck:
    def __init__(self):
        self.vec = [poly() for i in range(K)]


 

#/*************************************************
#* Name:        polyvecl_chknorm
#*
#* Description: Check infinity norm of polynomials in vector of length L.
#*              Assumes input polyvecl to be reduced by polyvecl_reduce().
#*
#* Arguments:   - const polyvecl *v: pointer to vector
#*              - int32_t B: norm bound
#*
#* Returns 0 if norm of all polynomials is strictly smaller than B <= (Q-1)/8
#* and 1 otherwise.
#**************************************************/

def polyvecl_chknorm(v, bound):
	for i in range(L):
		if(poly_chknorm(v.vec[i], bound)):
			return 1
	
	return 0
	

#/*************************************************
#* Name:        expand_mat
#*
#* Description: Implementation of ExpandA. Generates matrix A with uniformly
#*              random coefficients a_{i,j} by performing rejection
#*              sampling on the output stream of SHAKE128(rho|j|i)
#*              or AES256CTR(rho,j|i).
#*
#* Arguments:   - polyvecl mat[K]: output matrix
#*              - const uint8_t rho[]: byte array containing seed rho
#**************************************************/

def polyvec_matrix_expand(mat, rho): 
  
	for i in range(K):
		for j in range(L):
			poly_uniform(mat[i].vec[j], rho, (i << 8) + j)


#/*************************************************
#* Name:        polyvecl_ntt
#*
#* Description: Forward NTT of all polynomials in vector of length L. Output
#*              coefficients can be up to 16*Q larger than input coefficients.
#*
#* Arguments:   - polyvecl *v: pointer to input/output vector
#**************************************************/

def polyvecl_ntt(v):

	for i in range(L):
		poly_ntt(v.vec[i])
		
#/*************************************************
#* Name:        polyveck_ntt
#*
#* Description: Forward NTT of all polynomials in vector of length K. Output
#*              coefficients can be up to 16*Q larger than input coefficients.
#*
#* Arguments:   - polyveck *v: pointer to input/output vector
#**************************************************/

def polyveck_ntt(v):

	for i in range(K):
		poly_ntt(v.vec[i])
		
		
#/*************************************************
#* Name:        polyveck_shiftl
#*
#* Description: Multiply vector of polynomials of Length K by 2^D without modular
#*              reduction. Assumes input coefficients to be less than 2^{31-D}.
#*
#* Arguments:   - polyveck *v: pointer to input/output vector
#**************************************************/

def polyveck_shiftl(v):

	for i in range(K):
		poly_shiftl(v.vec[i])

		
		
#/*************************************************
#* Name:        polyvecl_pointwise_acc_montgomery
#*
#* Description: Pointwise multiply vectors of polynomials of length L, multiply
#*              resulting vector by 2^{-32} and add (accumulate) polynomials
#*              in it. Input/output vectors are in NTT domain representation.
#*
#* Arguments:   - poly *w: output polynomial
#*              - const polyvecl *u: pointer to first input vector
#*              - const polyvecl *v: pointer to second input vector
#**************************************************/

def polyvecl_pointwise_acc(w, u, v):

	t = poly()
	poly_pointwise(w, u.vec[0], v.vec[0])
	for i in range(1,L):
		poly_pointwise(t, u.vec[i], v.vec[i])
		poly_add(w, w, t)
		

def polyvec_matrix_pointwise(t, mat, v): 
  
	for i in range(K):
		polyvecl_pointwise_acc(t.vec[i], mat[i], v)
				
def polyveck_pointwise_poly(r, a, v): 

	for i in range(K):
		poly_pointwise(r.vec[i], a, v.vec[i])

#£/*************************************************
#* Name:        polyveck_sub
#*
#* Description: Subtract vectors of polynomials of length K.
#*              No modular reduction is performed.
#*
#* Arguments:   - polyveck *w: pointer to output vector
#*              - const polyveck *u: pointer to first input vector
#*              - const polyveck *v: pointer to second input vector to be
#*                                   subtracted from first input vector
#**************************************************/

def polyveck_sub(w, u, v): 
	for i in range(K):
		poly_sub(w.vec[i], u.vec[i], v.vec[i])
		
		
#/*************************************************
#* Name:        polyveck_invntt
#*
#* Description: Inverse NTT of polynomials
#*              in vector of length K. Input coefficients need to be less
#*              than 2*Q.
#*
#* Arguments:   - polyveck *v: pointer to input/output vector
#**************************************************/

def polyveck_invntt(v):
	
	for i in range(K):
		poly_invntt(v.vec[i])

#/*************************************************
#* Name:        polyveck_use_hint
#*
#* Description: Use hint vector to correct the high bits of input vector.
#*
#* Arguments:   - polyveck *w: pointer to output vector of polynomials with
#*                             corrected high bits
#*              - const polyveck *u: pointer to input vector
#*              - const polyveck *h: pointer to input hint vector
#**************************************************/

def polyveck_use_hint(w, u, h):

	for i in range(K):
		poly_use_hint(w.vec[i], u.vec[i], h.vec[i])


def polyveck_pack_w1(r, w1): 

	for i in range(K):
		polyw1_pack(r, i*POLYW1_PACKEDBYTES, w1.vec[i])
		
def polyvecl_print(v):
	
	print("[", end='')
	for i in range(L):
		poly_print(v.vec[i])
		if (i < L-1): print(",", end='')
	print("]")	
	
def polyveck_print(v):
	
	print("[", end='')
	for i in range(K):
		poly_print(v.vec[i])
		if (i < K-1): print(",", end='')
	print("]")			