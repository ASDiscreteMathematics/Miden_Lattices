from lattice_operations import *
import os

class Pke:
    def __init__(self, security_level, G=None):
        self.security_level = security_level
        if security_level <= 128:
            self.module_dimension = 3
        elif security_level <= 192:
            self.module_dimension = 4
        else:
            self.module_dimension = 6
        if G != None:
            self.G = G
        else:
            M = self.module_dimension
            randomness = os.urandom(9*64*M*M)
            randomness_chunks = [randomness[i*9*64:(i+1)*9*64] for i in range(M*M)]
            self.G = [[sample_random_polynomial(randomness_chunks[j*M+i]) for i in range(M)] for j in range(M)]

    def KeyGen(self):
        # get randomness
        randomness = os.urandom(2*self.module_dimension*64 * 2)
        randomness_chunks = [randomness[i*2*64:(i+1)*2*64] for i in range(2*self.module_dimension)]

        # sample secret key, which consists of two short elements
        s = [sample_short_polynomial(randomness_chunks[i]) for i in range(self.module_dimension)]
        e = [sample_short_polynomial(randomness_chunks[len(randomness_chunks)//2 + i]) for i in range(self.module_dimension)]

        # compute matching public key
        A = map_right(self.G, s, e)

        return (s,e), A

    def Enc(self, pk, msg):
        # get randomness
        randomness = os.urandom(2*self.module_dimension*64 * 2)
        randomness_chunks = [randomness[i*2*64:(i+1)*2*64] for i in range(2*self.module_dimension)]


        # sample ephemeral key, which consists of two short elements
        c = [sample_short_polynomial(randomness_chunks[i]) for i in range(self.module_dimension)]
        d = [sample_short_polynomial(randomness_chunks[len(randomness_chunks)//2 + i]) for i in range(self.module_dimension)]

        # compute matching Diffie-Hellman contribution
        B = map_left(self.G, c, d)

        # compute embedding of message
        assert(len(msg) == 256)
        embedding = embed_msg(msg)

        # compute noisy Diffie-Hellman value
        K = map_inner_product(pk, c)

        # add embedding of message
        for i in range(64):
            K[i] += embedding[i]

        return B, K

    def Dec(self, sk, ctxt):
        s, e = sk
        B, C = ctxt

        # compute noisy Diffie-Hellman value
        K = map_inner_product(s, B)
       
        # cancel it from the second ciphertext coordinate 
        for i in range(64):
            K[i] = C[i] - K[i]

        # extract message
        return extract_msg(K)

    def HomomorphicAdd(self, ctxt_left, ctxt_right):
        B_left, K_left = ctxt_left
        B_right, K_right = ctxt_right
        B = [[(B_left[i][j] + B_right[i][j]) % ((1 << 64) - (1 << 32) + 1) for j in range(64)] for i in range(self.module_dimension)]
        K = [(K_left[j] + K_right[j]) % ((1 << 64) - (1 << 32) + 1) for j in range(64)]
        return B, K

def test_pke():
    sec_lvls = [128, 192, 256]
    for sec_lvl in sec_lvls:
        scheme = Pke(sec_lvl)
        sk, pk = scheme.KeyGen()
        msg = [int(os.urandom(1)[0]) % 2 for i in range(256)]
        ctxt = scheme.Enc(pk, msg)
        dec = scheme.Dec(sk, ctxt)
        assert(dec == msg), "Fail"
    return True

def test_hom_add():
    sec_lvls = [128, 192, 256]
    for sec_lvl in sec_lvls:
        scheme = Pke(sec_lvl)
        sk, pk = scheme.KeyGen()
        msg1 = [int(os.urandom(1)[0]) % 2 for i in range(256)]
        msg2 = [int(os.urandom(1)[0]) % 2 for i in range(256)]
        msg = [(msg1[i] + msg2[i]) % 2 for i in range(256)]
        ctxt1 = scheme.Enc(pk, msg1)
        ctxt2 = scheme.Enc(pk, msg2)
        ctxt = scheme.HomomorphicAdd(ctxt1, ctxt2)
        dec = scheme.Dec(sk, ctxt)
        assert(dec == msg), "Fail"
    return True

