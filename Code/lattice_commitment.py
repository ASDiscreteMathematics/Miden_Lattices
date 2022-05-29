from lattice_operations import *
from CompactFIPS202 import SHAKE256 as XOF

# Re-randomizable Commitment Scheme
class Rcs:
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

    def Commit(self, message):
        # seed PRNG with message and expand into enough bytes
        randomness = XOF(message, 2*64*self.module_dimension*2)
        randomness_chunks = [randomness[2*64*i:2*64*(i+1)] for i in range(2*self.module_dimension)]

        # sample secrets pseudorandomly
        a = [sample_short_polynomial(randomness_chunks[i]) for i in range(self.module_dimension)]
        b = [sample_short_polynomial(randomness_chunks[i]) for i in range(self.module_dimension, 2*self.module_dimension)]

        # G*a + b
        A = map_right(self.G, a, b)

        return A

    def VerifyRaw(self, commitment, message):
        # seed PRNG with message and expand into enough bytes
        randomness = XOF(message, 2*64*self.module_dimension*2)
        randomness_chunks = [randomness[i*2*64:(i+1)*2*64] for i in range(2*self.module_dimension)]

        # sample secrets pseudorandomly
        a = [sample_short_polynomial(randomness_chunks[i]) for i in range(self.module_dimension)]
        b = [sample_short_polynomial(randomness_chunks[i]) for i in range(self.module_dimension, 2*self.module_dimension)]
        
        # G*a + b
        A = map_right(self.G, a, b)

        return A == commitment

    def Rerandomize(self, commitment):
        # get random bytes
        randomness = os.urandom(2*64*self.module_dimension*2)
        randomness_chunks = [randomness[i*2*64:(i+1)*2*64] for i in range(2*self.module_dimension)]

        # sample secrets
        c = [sample_short_polynomial(randomness_chunks[i]) for i in range(self.module_dimension)]
        d = [sample_short_polynomial(randomness_chunks[i]) for i in range(self.module_dimension, 2*self.module_dimension)]

        # B = c*G + d
        B = map_left(self.G, c, d)
        # C = c*(G*a + b)
        C = map_inner_product(c, commitment)

        return B, C

    def VerifyRerandomized(self, rerandomized_commitment, message):
        # seed PRNG with message and expand into enough bytes
        randomness = XOF(message, 2*64*self.module_dimension*2)
        randomness_chunks = [randomness[i*2*64:(i+1)*2*64] for i in range(2*self.module_dimension)]

        # sample secret (just the one) pseudorandomly
        a = [sample_short_polynomial(randomness_chunks[i]) for i in range(self.module_dimension)]
        a = [sample_short_polynomial(randomness[i*2*64:(i+1)*2*64]) for i in range(self.module_dimension)]
        
        # parse
        B, C = rerandomized_commitment

        # C_ = (c*G + d) * a
        C_ = map_inner_product(B, a)

        # If all is good,  C - C_ == c*b - d*a
        # which is small according to our chunk-wise norm.
        # In this case, every chunk of 16 bits will be smaller (in
        # absolute value) than 2^14. Attempting to extract the
        # embedded message will result in the all-zero "message".
        k = extract_msg([C[i] - C_[i] for i in range(64)])
        return k == [0]*256

def test_commitment():
    sec_lvls = [128]#, 192, 256]
    for sec_lvl in sec_lvls:
        scheme = Rcs(sec_lvl)
        msg = os.urandom(16)
        com = scheme.Commit(msg)
        vfy = scheme.VerifyRaw(com, msg)
        assert(vfy), "Fail raw"
        rcom = scheme.Rerandomize(com)
        rvfy = scheme.VerifyRerandomized(rcom, msg)
        assert(rvfy), "Fail rerandomized"
    return True

