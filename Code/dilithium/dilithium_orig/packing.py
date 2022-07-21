from params import *
from poly import *
from poly_vec import *

#/*************************************************
#* Name:        unpack_pk
#*
#* Description: Unpack public key pk = (rho, t1).
#*
#* Arguments:   - const uint8_t rho[]: output byte array for rho
#*              - const polyveck *t1: pointer to output vector t1
#*              - uint8_t pk[]: byte array containing bit-packed pk
#**************************************************/

def unpack_pk(rho, t1, pk):
	for i in range(SEEDBYTES):
		rho[i] = pk[i]
 
	for i in range(K):
		polyt1_unpack(t1.vec[i], pk[SEEDBYTES + i*POLYT1_PACKEDBYTES:])


#/*************************************************
#* Name:        unpack_sig
#*
#* Description: Unpack signature sig = (c, z, h).
#*
#* Arguments:   - uint8_t *c: pointer to output challenge hash
#*              - polyvecl *z: pointer to output vector z
#*              - polyveck *h: pointer to output hint vector h
#*              - const uint8_t sig[]: byte array containing
#*                bit-packed signature
#*
#* Returns 1 in case of malformed signature; otherwise 0.
#**************************************************/

def unpack_sig(c, z, h, sig):

	for i in range(SEEDBYTES):
		c[i] = sig[i]
	
	for i in range(L):
		polyz_unpack(z.vec[i], sig[SEEDBYTES + i*POLYZ_PACKEDBYTES:])
  
	offset = SEEDBYTES + L*POLYZ_PACKEDBYTES

	#/* Decode h */
	k = 0	
	for i in range(K):
		for j in range(N):
			h.vec[i].coeffs[j] = 0

		if(sig[offset + OMEGA + i] < k or sig[offset + OMEGA + i] > OMEGA):
			return 1
			
		for j in range(k, sig[offset + OMEGA + i]):
			#/* Coefficients are ordered for strong unforgeability */
			if (j > k and sig[offset + j] <= sig[offset + j-1]): 
				return 1
			h.vec[i].coeffs[sig[offset + j]] = 1
	
		k = sig[offset + OMEGA + i]

	#/* Extra indices are zero for strong unforgeability */
	for j in range(k, OMEGA):
		if(sig[offset + j] != 0):
			return 1
	
	return 0
		