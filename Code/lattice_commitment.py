from lattice_operations import *
from CompactFIPS202 import SHAKE256 as XOF
import os

# Re-randomizable Commitment Scheme
class NaiveRcs:
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
        # sample randomness
        randomness = os.urandom(32)

        # seed PRNG with message||randomness and expand into enough bytes
        uniform_bytes = XOF(message + randomness, 2*64*self.module_dimension*2)
        uniform_chunks = [uniform_bytes[2*64*i:2*64*(i+1)] for i in range(2*self.module_dimension)]

        # sample secrets pseudorandomly
        a = [sample_short_polynomial(uniform_chunks[i]) for i in range(self.module_dimension)]
        b = [sample_short_polynomial(uniform_chunks[i]) for i in range(self.module_dimension, 2*self.module_dimension)]

        # G*a + b
        A = map_right(self.G, a, b)

        return A, randomness

    def VerifyRaw(self, message, randomness, commitment):
        # seed PRNG with message and expand into enough bytes
        uniform_bytes = XOF(message+randomness, 2*64*self.module_dimension*2)
        uniform_chunks = [uniform_bytes[i*2*64:(i+1)*2*64] for i in range(2*self.module_dimension)]

        # sample secrets pseudorandomly
        a = [sample_short_polynomial(uniform_chunks[i]) for i in range(self.module_dimension)]
        b = [sample_short_polynomial(uniform_chunks[i]) for i in range(self.module_dimension, 2*self.module_dimension)]
        
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

    def VerifyRerandomized(self, message, randomness, rerandomized_commitment):
        # seed PRNG with message and expand into enough bytes
        uniform_bytes = XOF(message + randomness, 2*64*self.module_dimension*2)
        uniform_chunks = [uniform_bytes[i*2*64:(i+1)*2*64] for i in range(2*self.module_dimension)]

        # sample secret (just the one) pseudorandomly
        a = [sample_short_polynomial(uniform_chunks[i]) for i in range(self.module_dimension)]
        #a = [sample_short_polynomial(randomness[i*2*64:(i+1)*2*64]) for i in range(self.module_dimension)]
        
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

class OptimizedRcs:
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
            self.Gf = [[gg for gg in g] for g in self.G]
            for g in self.Gf:
                for gg in g:
                    ntt_4_64.NTT_64(gg)
        else:
            M = self.module_dimension
            randomness = os.urandom(9*64*M*M)
            randomness_chunks = [randomness[i*9*64:(i+1)*9*64] for i in range(M*M)]
            self.Gf = [[sample_random_polynomial(randomness_chunks[j*M+i]) for i in range(M)] for j in range(M)]
            self.G = [[gg for gg in g] for g in self.Gf]
            for g in self.G:
                for gg in g:
                    ntt_4_64.iNTT_64(gg)

    def Commit(self, message):
        # sample randomness
        randomness = os.urandom(32)

        # seed PRNG with message||randomness and expand into enough bytes
        uniform_bytes = XOF(message + randomness, 2*64*self.module_dimension*2)
        uniform_chunks = [uniform_bytes[2*64*i:2*64*(i+1)] for i in range(2*self.module_dimension)]

        # sample secrets pseudorandomly
        af = [sample_short_polynomial(uniform_chunks[i]) for i in range(self.module_dimension)]
        bf = [sample_short_polynomial(uniform_chunks[i]) for i in range(self.module_dimension, 2*self.module_dimension)]

        # map to frequency domain
        for i in range(self.module_dimension):
            ntt_4_64.NTT_64(af[i])
            ntt_4_64.NTT_64(bf[i])

        # G*a + b
        Af = map_right_had(self.Gf, af, bf)

        return Af, randomness

    def VerifyRaw(self, message, randomness, commitment):
        # seed PRNG with message and expand into enough bytes
        uniform_bytes = XOF(message + randomness, 2*64*self.module_dimension*2)
        uniform_chunks = [uniform_bytes[i*2*64:(i+1)*2*64] for i in range(2*self.module_dimension)]

        # sample secrets pseudorandomly
        af = [sample_short_polynomial(uniform_chunks[i]) for i in range(self.module_dimension)]
        bf = [sample_short_polynomial(uniform_chunks[i]) for i in range(self.module_dimension, 2*self.module_dimension)]
 
        # map to frequency domain
        for i in range(self.module_dimension):
            ntt_4_64.NTT_64(af[i])
            ntt_4_64.NTT_64(bf[i])
       
        # G*a + b
        Af = map_right_had(self.Gf, af, bf)

        return Af == commitment

    def Rerandomize(self, commitment):
        # get random bytes
        randomness = os.urandom(2*64*self.module_dimension*2)
        randomness_chunks = [randomness[i*2*64:(i+1)*2*64] for i in range(2*self.module_dimension)]

        # sample secrets
        cf = [sample_short_polynomial(randomness_chunks[i]) for i in range(self.module_dimension)]
        df = [sample_short_polynomial(randomness_chunks[i]) for i in range(self.module_dimension, 2*self.module_dimension)]

        # map to frequency domain
        for i in range(self.module_dimension):
            ntt_4_64.NTT_64(cf[i])
            ntt_4_64.NTT_64(df[i])
    
        # B = c*G + d
        Bf = map_left_had(self.Gf, cf, df)
        # C = c*(G*a + b)
        Cf = map_inner_product_had(cf, commitment)

        return Bf, Cf

    def VerifyRerandomized(self, message, randomness, rerandomized_commitment):
        # seed PRNG with message and expand into enough bytes
        uniform_bytes = XOF(message + randomness, 2*64*self.module_dimension*2)
        uniform_chunks = [uniform_bytes[i*2*64:(i+1)*2*64] for i in range(2*self.module_dimension)]

        # sample secret (just the one) pseudorandomly
        af = [sample_short_polynomial(uniform_chunks[i]) for i in range(self.module_dimension)]

        # map to frequency domain
        for i in range(self.module_dimension):
            ntt_4_64.NTT_64(af[i])
        
        # parse
        Bf, Cf = rerandomized_commitment

        # C_ = (c*G + d) * a
        C_ = map_inner_product_had(Bf, af)

        # Diff = C - C_
        Difff = [Cf[i] - C_[i] for i in range(64)]

        # map to time domain
        ntt_4_64.iNTT_64(Difff)

        # If all is good,  C - C_ == c*b - d*a
        # which is small according to our chunk-wise norm.
        # In this case, every chunk of 16 bits will be smaller (in
        # absolute value) than 2^14. Attempting to extract the
        # embedded message will result in the all-zero "message".
        k = extract_msg(Difff)
        return k == [0]*256

def test_commitment():
    sec_lvls = [128, 192, 256]
    for sec_lvl in sec_lvls:
        #scheme = NaiveRcs(sec_lvl)
        scheme = OptimizedRcs(sec_lvl)
        msg = os.urandom(16)
        com, decom = scheme.Commit(msg)
        vfy = scheme.VerifyRaw(msg, decom, com)
        assert(vfy), "Fail raw"
        rcom = scheme.Rerandomize(com)
        rvfy = scheme.VerifyRerandomized(msg, decom, rcom)
        assert(rvfy), "Fail rerandomized"
    return True

