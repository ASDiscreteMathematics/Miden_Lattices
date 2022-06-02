from params import *
from poly import *
from poly_vec import *
from packing import *
from Crypto.Hash import SHAKE256

#/*************************************************
#* Name:        crypto_sign_verify
#*
#* Description: Verifies signature.
#*
#* Arguments:   - uint8_t *m: pointer to input signature
#*              - size_t siglen: length of signature
#*              - const uint8_t *m: pointer to message
#*              - size_t mlen: length of message
#*              - const uint8_t *pk: pointer to bit-packed public key
#*
#* Returns 0 if signature could be verified correctly and -1 otherwise
#**************************************************/

def crypto_sign_verify(sig, siglen, m, mlen, pk):

  unsigned int i;
  uint8_t buf[K*POLYW1_PACKEDBYTES];
  uint8_t rho[SEEDBYTES];
  uint8_t mu[CRHBYTES];
  uint8_t c[SEEDBYTES];
  uint8_t c2[SEEDBYTES];
  poly cp;
  polyvecl mat[K], z;
  polyveck t1, w1, h;
  keccak_state state;

	if(siglen != CRYPTO_BYTES):
		return -1

	# unpackig public key

	rho = bytearray(SEEDBYTES)
	t1 = polyveck()
	unpack_pk(rho, t1, pk)
  
	# unpacking signature
  
    c = bytearray(SEEDBYTES)
	z = polyvecl()
	h = polyveck()
  
	if(unpack_sig(c, z, h, sig) == 1):
		return -1
		
	if(polyvecl_chknorm(&z, GAMMA1 - BETA)):
		return -1

	shake = SHAKE256.new()
    shake.update(pk)
    mu = shake.read(SEEDBYTES)

	shake = SHAKE256.new()	
	shake.update(mu)
	shake.update(m)
	mu = shake.read(CRHBYTES)

	cp = poly()
	poly_challenge(cp, c)
	polyvec_matrix_expand(mat, rho)
	
	polyvecl_ntt(z)
	w1 = polyveck()
	mat = [polyvecl() for j in range(K)]
	polyvec_matrix_pointwise(w1, mat, z)

	poly_ntt(cp)
	
	polyveck_shiftl(t1)
	polyveck_ntt(t1)
	polyveck_pointwise_poly(t1, cp, t1)

	polyveck_sub(w1, w1, t1)
	polyveck_invntt(w1)

	/* Reconstruct w1 */
	
	polyveck_use_hint(w1, w1, h)
	polyveck_pack_w1(buf, w1)



	/* Call random oracle and verify challenge */
	
	shake = SHAKE256.new()
	shake.update(mu[0:CRHBYTES])
	shake.update(buf[0:K*POLYW1_PACKEDBYTES]
	
	c2 = shake.read(SEEDBYTES)
	
	for i in range(SEEDBYTES):
		if(c[i] != c2[i])
			return -1;

	return 0;
