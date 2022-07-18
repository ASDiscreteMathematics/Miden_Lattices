"""
This file implements tests for various parts of the Falcon.py library.

Test the code with:
> make test
"""
from falcon import SecretKey, PublicKey, Params
from falcon import SALT_LEN, HEAD_LEN, SHAKE256
from scripts.sign_KAT import sign_KAT

def test_sign_KAT():
	message = b"data1"
	for n in sign_KAT:
		print("Testing n = ", n)
		sign_KAT_n = sign_KAT[n]
		for D in sign_KAT_n:
			f = D["f"]
			g = D["g"]
			F = D["F"]
			G = D["G"]
			sk = SecretKey(n, [f, g, F, G])
			pk = PublicKey(sk)
			sig = bytes.fromhex(D["sig"])
			if pk.verify(message, sig) is False:
				return False
	return True


# Run all the tests
if (__name__ == "__main__"):
    print("Test Sig KATs")
    print("Test is OK" if (test_sign_KAT() is True) else "Not OK")