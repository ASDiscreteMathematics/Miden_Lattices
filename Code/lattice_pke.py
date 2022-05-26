from lattice_operations import *

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
            self.G = [[sample_random_polynomial() for i in range(self.module_dimension)] for j in range(self.module_dimension)]

    def KeyGen(self):
        # sample secret key, which consists of two short elements
        s = [sample_short_polynomial() for i in range(self.module_dimension)]
        e = [sample_short_polynomial() for i in range(self.module_dimension)]

        # compute matching public key
        A = map_right(self.G, s, e)

        return (s,e), A

    def Enc(self, pk, msg):
        # sample ephemeral key which consists of two short elements
        c = [sample_short_polynomial() for i in range(self.module_dimension)]
        d = [sample_short_polynomial() for i in range(self.module_dimension)]

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

def test_pke():
    sec_lvls = [128]#, 192, 256]
    for sec_lvl in sec_lvls:
        scheme = Pke(sec_lvl)
        sk, pk = scheme.KeyGen()
        msg = [int(os.urandom(1)[0]) % 2 for i in range(256)]
        ctxt = scheme.Enc(pk, msg)
        dec = scheme.Dec(sk, ctxt)
        assert(dec == msg), "Fail"
    return True

