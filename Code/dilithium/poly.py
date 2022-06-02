from params import *
from rounding import *
from ntt import *
from Crypto.Hash import SHAKE256
from Crypto.Hash import SHAKE128

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


#/*************************************************
#* Name:        challenge
#*
#* Description: Implementation of H. Samples polynomial with TAU nonzero
#*              coefficients in {-1,1} using the output stream of
#*              SHAKE256(seed).
#*
#* Arguments:   - poly *c: pointer to output polynomial
#*              - const uint8_t mu[]: byte array containing seed of length SEEDBYTES
#**************************************************/

def poly_challenge(c, seed):

	shake = SHAKE256.new()
	shake.update(seed[0:SEEDBYTES])
	buf = shake.read(SHAKE256_RATE)
	
	signs = 0
	for i in range(8):
		signs |= buf[i] << 8*i
	pos = 8;

	for i in range(N):
		c.coeffs[i] = 0;

	for i in range(N-TAU, N):
		if(pos >= SHAKE256_RATE):
			buf = shake.read(SHAKE256_RATE)
			pos = 0
		b = buf[pos]
		pos += 1
		
		while(b > i):
			if(pos >= SHAKE256_RATE):
				buf = shake.read(SHAKE256_RATE)
				pos = 0
			b = buf[pos]
			pos += 1

		c.coeffs[i] = c.coeffs[b]
		c.coeffs[b] = 1 - 2*(signs & 1)
		signs >>= 1
  
#/*************************************************
#* Name:        rej_uniform
#*
#* Description: Sample uniformly random coefficients in [0, Q-1] by
#*              performing rejection sampling on array of random bytes.
#*
#* Arguments:   - int32_t *a: pointer to output array (allocated)
#*              - unsigned int len: number of coefficients to be sampled
#*              - const uint8_t *buf: array of random bytes
#*              - unsigned int buflen: length of array of random bytes
#*
#* Returns number of sampled coefficients. Can be smaller than len if not enough
#* random bytes were given.
#**************************************************/

def rej_uniform(a, start, len, buf, buflen):

	ctr = pos = 0
	while(ctr < len and pos + 3 <= buflen):
		t  = buf[pos]
		pos += 1
		t |= buf[pos] << 8
		pos += 1
		t |= buf[pos] << 16
		pos += 1
		t &= 0x7FFFFF

		if(t < Q):
			a[start + ctr] = t
			ctr += 1
  
	return ctr


#/*************************************************
#* Name:        poly_uniform
#*
#* Description: Sample polynomial with uniformly random coefficients
#*              in [0,Q-1] by performing rejection sampling on the
#*              output stream of SHAKE256(seed|nonce)
#*
#* Arguments:   - poly *a: pointer to output polynomial
#*              - const uint8_t seed[]: byte array with seed of length SEEDBYTES
#*              - uint16_t nonce: 2-byte nonce
#**************************************************/


POLY_UNIFORM_NBLOCKS = ((768 + STREAM128_BLOCKBYTES - 1)//STREAM128_BLOCKBYTES)

def poly_uniform(a, seed, nonce):

	print(POLY_UNIFORM_NBLOCKS)

	buflen = POLY_UNIFORM_NBLOCKS*STREAM128_BLOCKBYTES;

	t = nonce.to_bytes(2, 'little')
	shake = SHAKE128.new()
	shake.update(seed[0:SEEDBYTES])
	shake.update(t)

	buf = bytearray(POLY_UNIFORM_NBLOCKS*STREAM128_BLOCKBYTES + 2)
	for i in range(POLY_UNIFORM_NBLOCKS):
		buf[i*SHAKE128_RATE:(i+1)*SHAKE128_RATE]= bytearray(shake.read(SHAKE128_RATE))

	ctr = rej_uniform(a.coeffs, 0, N, buf, buflen)

	while(ctr < N):
		off = buflen % 3
		for i in range(off):
			buf[i] = buf[buflen - off + i]

		buf[off:off+SHAKE128_RATE]= bytearray(shake.read(SHAKE128_RATE))
		buflen = STREAM128_BLOCKBYTES + off
		ctr += rej_uniform(a.coeffs, ctr, N - ctr, buf, buflen)
  

#/*************************************************
#* Name:        poly_pointwise_montgomery
#*
#* Description: Pointwise multiplication of polynomials in NTT domain
#*              representation and multiplication of resulting polynomial
#*              by 2^{-32}.
#*
#* Arguments:   - poly *c: pointer to output polynomial
#*              - const poly *a: pointer to first input polynomial
#*              - const poly *b: pointer to second input polynomial
#**************************************************/

def poly_pointwise(c, a, b):
	
	for i in range(N):
		c.coeffs[i] = (a.coeffs[i] * b.coeffs[i]) %q


#*************************************************
# Name:        poly_ntt
#
# Description: Inplace forward NTT. Coefficients can grow by
#              8*Q in absolute value.
#*
#* Arguments:   - poly *a: pointer to input/output polynomial
#**************************************************/

def poly_ntt(a): 
	NTT_256(a.coeffs)
	
#/*************************************************
#* Name:        poly_sub
#*
#* Description: Subtract polynomials. No modular reduction is
#*              performed.
#*
#* Arguments:   - poly *c: pointer to output polynomial
#*              - const poly *a: pointer to first input polynomial
#*              - const poly *b: pointer to second input polynomial to be
#*                               subtraced from first input polynomial
#**************************************************/

def poly_sub(c, a, b): 
	
	for i in range(N):
		c.coeffs[i] = (a.coeffs[i] - b.coeffs[i])%q

	
#/*************************************************
#* Name:        poly_shiftl
#*
#* Description: Multiply polynomial by 2^D without modular reduction. Assumes
#*              input coefficients to be less than 2^{31-D} in absolute value.
#*
#* Arguments:   - poly *a: pointer to input/output polynomial
#**************************************************/

def poly_shiftl(a): 
  
  for i in range(N):
	a.coeffs[i] <<= D


#/*************************************************
#* Name:        poly_invntt
#*
#* Description: Inplace inverse NTT 
#*              Input coefficients need to be less than Q in absolute
#*              value and output coefficients are again bounded by Q.
#*
#* Arguments:   - poly *a: pointer to input/output polynomial
#**************************************************/

def poly_invntt(a):
		
	iNTT_256(a.coeffs)

	
#/*************************************************
#* Name:        poly_use_hint
#*
#* Description: Use hint polynomial to correct the high bits of a polynomial.
#*
#* Arguments:   - poly *b: pointer to output polynomial with corrected high bits
#*              - const poly *a: pointer to input polynomial
#*              - const poly *h: pointer to input hint polynomial
#**************************************************/

def poly_use_hint(b, a, h): 

	for i in range(N):
		b.coeffs[i] = use_hint(a.coeffs[i], h.coeffs[i])

#/*************************************************
#* Name:        polyw1_pack
#*
#* Description: Bit-pack polynomial w1 with coefficients in [0,15] or [0,43].
#*              Input coefficients are assumed to be standard representatives.
#*
#* Arguments:   - uint8_t *r: pointer to output byte array with at least
#*                            POLYW1_PACKEDBYTES bytes
#*              - const poly *a: pointer to input polynomial
#**************************************************/

def polyw1_pack(r, a): 

	if (GAMMA2 == (Q-1)/88):
		for i in range(N >> 2):
			r[3*i+0]  = a.coeffs[4*i+0];
			r[3*i+0] |= a.coeffs[4*i+1] << 6
			r[3*i+1]  = a.coeffs[4*i+1] >> 2
			r[3*i+1] |= a.coeffs[4*i+2] << 4
			r[3*i+2]  = a.coeffs[4*i+2] >> 4
			r[3*i+2] |= a.coeffs[4*i+3] << 2
	elif (GAMMA2 == (Q-1)/32):
		for in range(N >> 1):
			r[i] = a.coeffs[2*i+0] | (a->coeffs[2*i+1] << 4)