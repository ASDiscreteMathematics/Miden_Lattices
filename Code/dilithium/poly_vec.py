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

