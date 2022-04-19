from params import *

class poly:
    def __init__(self):
        self.coeffs = [0 for i in range(N)]

#/*************************************************
#* Name:        poly_chknorm
#*
#* Description: Check infinity norm of polynomial against given bound.
#*              Assumes input coefficients were reduced by reduce32().
#*
#* Arguments:   - const poly *a: pointer to polynomial
#*              - int32_t B: norm bound
#*
#* Returns 0 if norm is strictly smaller than B <= (Q-1)/8 and 1 otherwise.
#**************************************************/

def poly_chknorm(a, B):
	if(B > (Q-1)/8):
		return 1

	#/* It is ok to leak which coefficient violates the bound since
	#	the probability for each coefficient is independent of secret
	#	data but we must not leak the sign of the centralized representative. */
	
	for i in range(N):
		#/* Absolute value */
		t = a.coeffs[i] >> 31;
		t = a.coeffs[i] - (t & 2*a.coeffs[i]);

		if(t >= B):
			return 1
	
	return 0


#/*************************************************
#* Name:        polyt1_unpack
#*
#* Description: Unpack polynomial t1 with 10-bit coefficients.
#*              Output coefficients are standard representatives.
#*
#* Arguments:   - poly *r: pointer to output polynomial
#*              - const uint8_t *a: byte array with bit-packed polynomial
#**************************************************/

def polyt1_unpack(r, a):
	for i in range (N >> 2):
		r.coeffs[4*i+0] = ((a[5*i+0] >> 0) | (a[5*i+1] << 8)) & 0x3FF
		r.coeffs[4*i+1] = ((a[5*i+1] >> 2) | (a[5*i+2] << 6)) & 0x3FF
		r.coeffs[4*i+2] = ((a[5*i+2] >> 4) | (a[5*i+3] << 4)) & 0x3FF
		r.coeffs[4*i+3] = ((a[5*i+3] >> 6) | (a[5*i+4] << 2)) & 0x3FF

#/*************************************************
#* Name:        polyz_unpack
#*
#* Description: Unpack polynomial z with coefficients
#*              in [-(GAMMA1 - 1), GAMMA1].
#*
#* Arguments:   - poly *r: pointer to output polynomial
#*              - const uint8_t *a: byte array with bit-packed polynomial
#**************************************************/

def polyz_unpack(r, a):
	if GAMMA1 == (1 << 17):
		for i in range(N >> 2):
			r.coeffs[4*i+0]  = a[9*i+0];
			r.coeffs[4*i+0] |= a[9*i+1] << 8;
			r.coeffs[4*i+0] |= a[9*i+2] << 16;
			r.coeffs[4*i+0] &= 0x3FFFF;

			r.coeffs[4*i+1]  = a[9*i+2] >> 2;
			r.coeffs[4*i+1] |= a[9*i+3] << 6;
			r.coeffs[4*i+1] |= a[9*i+4] << 14;
			r.coeffs[4*i+1] &= 0x3FFFF;

			r.coeffs[4*i+2]  = a[9*i+4] >> 4;
			r.coeffs[4*i+2] |= a[9*i+5] << 4;
			r.coeffs[4*i+2] |= a[9*i+6] << 12;
			r.coeffs[4*i+2] &= 0x3FFFF;

			r.coeffs[4*i+3]  = a[9*i+6] >> 6;
			r.coeffs[4*i+3] |= a[9*i+7] << 2;
			r.coeffs[4*i+3] |= a[9*i+8] << 10;
			r.coeffs[4*i+3] &= 0x3FFFF;

			r.coeffs[4*i+0] = GAMMA1 - r.coeffs[4*i+0];
			r.coeffs[4*i+1] = GAMMA1 - r.coeffs[4*i+1];
			r.coeffs[4*i+2] = GAMMA1 - r.coeffs[4*i+2];
			r.coeffs[4*i+3] = GAMMA1 - r.coeffs[4*i+3];
			
	elif GAMMA1 == (1 << 19):
		for i in range(N >> 1):
			r.coeffs[2*i+0]  = a[5*i+0];
			r.coeffs[2*i+0] |= a[5*i+1] << 8;
			r.coeffs[2*i+0] |= a[5*i+2] << 16;
			r.coeffs[2*i+0] &= 0xFFFFF;

			r.coeffs[2*i+1]  = a[5*i+2] >> 4;
			r.coeffs[2*i+1] |= a[5*i+3] << 4;
			r.coeffs[2*i+1] |= a[5*i+4] << 12;
			r.coeffs[2*i+0] &= 0xFFFFF;

			r.coeffs[2*i+0] = GAMMA1 - r.coeffs[2*i+0];
			r.coeffs[2*i+1] = GAMMA1 - r.coeffs[2*i+1];

