import os
import ntt_4_64

def ntt(array, root):
    if len(array) == 1:
        return array

    evens = ntt(array[::2], root*root % ((1<<64)-(1<<32)+1))
    odds = ntt(array[1::2], root*root % ((1<<64)-(1<<32)+1))
    half = len(array) // 2
    
    ret = [0] * len(array)
    acc = 1
    for i in range(len(array)):
        ret[i] = (evens[i % half] + acc * odds[i % half]) % ((1 << 64) - (1 << 32) + 1)
        acc = (acc * root) % ((1 << 64) - (1 << 32) + 1)
    return ret

def ring_mul(lhs, rhs):
    lhs_copy = [l for l in lhs]
    rhs_copy = [r for r in rhs]
    ntt_4_64.FastMul_64(lhs_copy, rhs_copy)
    return lhs_copy

def ring_mul_(lhs, rhs):
    # copy lhs and rhs
    # otherwise python will mutate in-place
    lhs = lhs[:]
    rhs = rhs[:]

    p128 = 17870292113338400769
    p128inv = 18442240469787213841
    p64 = 8
    p64inv = 16140901060737761281
    ninv = 18158513693329981441
    
    # scale
    acc = 1
    for i in range(64):
        lhs[i] = (lhs[i] * acc) % ((1<<64) - (1<<32) + 1)
        rhs[i] = (rhs[i] * acc) % ((1<<64) - (1<<32) + 1)
        acc = (acc * p128) % ((1 << 64) - (1<<32) + 1)

    # ntt
    lhs_ntt = ntt(lhs, p64)
    rhs_ntt = ntt(rhs, p64)

    # coordinate-wise product, and anticipate scale by n^-1 for intt
    prod = [(ninv*l*r) % ((1<<64) - (1<<32) + 1) for l,r in zip(lhs_ntt, rhs_ntt)]

    # inverse ntt
    product = ntt(prod, p64inv)

    # unscale
    acc = 1
    for i in range(64):
        product[i] = (product[i] * acc) % ((1<<64) - (1<<32) + 1)
        acc = (acc*p128inv) % ((1<<64) - (1<<32) + 1)

    return product

def ring_add(lhs, rhs):
    return [(l+r) % ((1<<64) - (1<<32) + 1) for (l, r) in zip(lhs, rhs)]

def ring_had(lhs, rhs):
    return [(l*r) % ((1<<64) - (1<<32) + 1) for (l, r) in zip(lhs, rhs)]

def sample_random_field_element(randomness):
    integer = sum((1 << 8*i) * int(randomness[i]) for i in range(9))
    return integer % ((1 << 64) - (1 << 32) + 1) # not perfectly uniform, but that's okay

def sample_random_polynomial(randomness):
    return [sample_random_field_element(randomness[i*9:(i+1)*9]) for i in range(64)]

def num_set_bits(integer):
    num_set_bits = 0
    # not constant time
    while integer:
        integer &= integer - 1
        num_set_bits += 1
    return num_set_bits

def sample_short_field_element(randomness):
    return sum((1 << (16*i)) * (num_set_bits(randomness[0]) - num_set_bits(randomness[1])) for i in range(4))

def sample_short_polynomial(randomness):
    return [sample_short_field_element(randomness[2*i:2*(i+1)]) for i in range(64)]

def map_right(G, a, b):
    module_dimension = len(G)
    v = [[0] * 64] * module_dimension
    for i in range(module_dimension):
        for j in range(module_dimension):
            prod = ring_mul(G[i][j][:], a[j][:])
            v[i] = ring_add(v[i], prod)
        v[i] = ring_add(v[i], b[i])
    return v

def map_left(G, a, b):
    module_dimension = len(G[0])
    v = [[0] * 64] * module_dimension
    for i in range(module_dimension):
        for j in range(module_dimension):
            prod = ring_mul(a[j], G[j][i])
            v[i] = ring_add(v[i], prod)
        v[i] = ring_add(v[i], b[i])
    return v

def map_inner_product(vector, a):
    s = [0] * 64
    module_dimension = len(a)
    for i in range(module_dimension):
        prod = ring_mul(vector[i][:], a[i][:])
        s = ring_add(s, prod)
    return s

def map_right_had(G, a, b):
    module_dimension = len(G)
    v = [[0] * 64] * module_dimension
    for i in range(module_dimension):
        for j in range(module_dimension):
            prod = ring_had(G[i][j][:], a[j][:])
            v[i] = ring_add(v[i], prod)
        v[i] = ring_add(v[i], b[i])
    return v

def map_left_had(G, a, b):
    module_dimension = len(G[0])
    v = [[0] * 64] * module_dimension
    for i in range(module_dimension):
        for j in range(module_dimension):
            prod = ring_had(a[j], G[j][i])
            v[i] = ring_add(v[i], prod)
        v[i] = ring_add(v[i], b[i])
    return v

def map_inner_product_had(vector, a):
    s = [0] * 64
    module_dimension = len(a)
    for i in range(module_dimension):
        prod = ring_had(vector[i][:], a[i][:])
        s = ring_add(s, prod)
    return s

def embed_msg(msg):
    assert(all(m == 0 or m == 1 for m in msg))
    embedding = [0]*(len(msg)//4)
    for i in range(len(msg)//4):
        for j in range(4):
            embedding[i] += msg[i*4+j] << (15 + 16 * j)

    return embedding

def extract_msg(embedding):
    msg = []
    for e in embedding:
        e = e % ((1<<64) - (1<<32) + 1)
        for i in range(4):
            chunk = e & 0xffff
            e >>= 16

            if chunk < (1<<14) or ((1<<16) - chunk) < (1<<14):
                msg += [0]
            else:
                msg += [1]
    return msg

def test_ntt_linearity():
    left = [3,2,1,0] + list(range(60))
    right = list(range(61)) + [1,2,3]
    
    p64 = 8

    left_ntt = ntt(left, p64)
    right_ntt = ntt(right, p64)

    left_ntt_plus_right_ntt = [(2*l + 3*r) % ((1<<64) - (1<<32) + 1) for (l, r) in zip(left_ntt, right_ntt)]
    left_plus_right = [(2*l + 3*r) % ((1<<64) - (1<<32) + 1) for (l, r) in zip(left, right)]
    left_plus_right_ntt = ntt(left_plus_right, p64)

    assert(left_ntt_plus_right_ntt == left_plus_right_ntt), "Fail"
    return True

def test_intt():
    array = [5,3,1] + list(range(61))

    p64 = 8
    p64inv = 16140901060737761281
    ninv = 18158513693329981441

    array_ntt = ntt(array, p64)
    array_ntt_ninv = [(a*ninv) % ((1<<64) - (1<<32) + 1) for a in array_ntt]
    array_ntt_intt = ntt(array_ntt_ninv, p64inv)

    assert(array == array_ntt_intt), "Fail"
    return True

def test_embed_extract():
    msg = [int(os.urandom(1)[0])%2 for i in range(256)]

    embedded = embed_msg(msg)

    perturbation = [0]*64
    for i in range(64):
        for j in range(4):
            noise_term = (int(os.urandom(1)[0])*256 + int(os.urandom(1)[0])) % (1<<14)
            if (os.urandom(1)[0])%2 == 1:
                noise_term *= -1
            perturbation[i] += noise_term * (1 << (16*j))

    noisy = [e + 0 for (e,p) in zip(embedded, perturbation)]

    decoded = extract_msg(noisy)

    assert(msg == decoded), "Fail"
    return True

def test_ring_mul():
    lhs = [1,2,3,4] + [0]*60
    rhs = list(reversed(list(range(10)))) + [0]*54
    expected = [9, 26, 50, 80, 70, 60, 50, 40, 30, 20, 11, 4] + [0] * 52

    observed = ring_mul(lhs[:], rhs[:])
    assert(expected == observed), "Fail"

    rhs_ = rhs[:-1] + [1]
    expected_ = expected[:]
    expected_[0] -= 2
    expected_[1] -= 3
    expected_[2] -= 4
    expected_[-1] = 1

    observed_ = ring_mul(lhs, rhs_)
    print(expected_)
    print(observed_)
    assert(expected_ == observed_), "Fail"

    return True

def test_ring_commutativity():
    a = sample_random_polynomial(os.urandom(9*64))
    b = sample_random_polynomial(os.urandom(9*64))
    a = list(range(5)) + [0]*59
    b = list(range(7)) + [0]*57
    ab = ring_mul(a,b)
    ba = ring_mul(b,a)
    assert(ab == ba), "Fail"
    return True

def test_ring_associativity():
    a = sample_random_polynomial(os.urandom(9*64))
    b = sample_random_polynomial(os.urandom(9*64))
    c = sample_random_polynomial(os.urandom(9*64))

    a = list(range(10)) + [0]*54
    b = list(reversed(range(5))) + [0]*59
    c = [1]*7 + [0]*57

    a_bc = ring_mul(a, ring_mul(b, c))
    ab_c = ring_mul(ring_mul(a, b), c)

    assert(all(a < ((1<<64) - (1<<32) + 1) and a >= 0 for a in a_bc))
    assert(all(a < ((1<<64) - (1<<32) + 1) and a >= 0 for a in ab_c))
    assert(a_bc == ab_c), "Fail"
    return True

def test_map_associativity():
    dimension = 4
    G = [[sample_random_polynomial(os.urandom(9*64)) for j in range(dimension)] for i in range(dimension)]
    a = [sample_random_polynomial(os.urandom(9*64)) for i in range(dimension)]
    b = [sample_random_polynomial(os.urandom(9*64)) for i in range(dimension)]

    zeros = [[0]*64]*dimension

    aG = map_left(G, a, zeros)
    aG_b = map_inner_product(aG, b) 
    Gb = map_right(G, b, zeros)
    a_Gb = map_inner_product(a, Gb)

    assert(a_Gb == aG_b), "Fail"

