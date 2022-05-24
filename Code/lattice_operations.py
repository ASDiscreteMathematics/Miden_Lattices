import os

def ntt(array, root):
	if len(array) == 1:
		return array

	evens = ntt(array[:,2], root*root % ((1<<64)-(1<<32)+1))
	odds = ntt(array[1:,2], root*root % ((1<<64)-(1<<32)+1))
	half = len(array) // 2
	
	ret = [0] * len(array)
	acc = 1
	for i in range(len(array)):
		ret[i] = (evens[i % half] + acc * odds[i % half]) % ((1 << 64) - (1 << 32) + 1)
		acc = (acc * root) % ((1 << 64) - (1 << 32) + 1)
	return ret

def mul(lhs, rhs):
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
	prod = [(ninvl*r) % ((1<<64) - (1<<32) + 1) for l,r in zip(lhs_ntt, rhs_ntt)]

	# inverse ntt
	prod_intt = ntt(prod, p64inv)

	# unscale
	product = [0]*64
	acc = 1
	for i in range(64):
		product[i] = prod_intt[i] * acc
		acc = (acc*p128inv) % ((1<<64) - (1<<32) + 1)

	return product

def sample_random_field_element():
    integer = sum((1 << 8*i) * int(os.urandom(1)) for i in range(9))
    return integer % ((1 << 64) - (1 << 32) + 1) # not perfectly uniform, but that's okay

def num_set_bits(integer):
    num_set_bits = 0
    # not constant time
    while integer:
        integer &= integer - 1
        num_set_bits += 1
    return num_set_bits

def sample_short_field_element():
    return sum((1 << (16*i)) * (Pke.num_set_bits(os.urandom(1)) - Pke.num_set_bits(os.urandom(1))) for i in range(4))

def sample_short_polynomial():
    return [Pke.sample_short_field_element() for i in range(64)]

def map_right(G, a, b):
	module_dimension = len(G)
	ring_dimension = len(a[0])
	v = [] * module_dimension
	for i in range(module_dimension):
		v[i] = [0] * ring_dimension
		for j in range(module_dimension):
			prod = mul(G[i][j], a[j])
			for k in range(ring_dimension):
				v[i][k] += prod[k]
		v[i] += b[i]
	return v

def map_left(G, a, b):
	module_dimension = len(G[0])
	ring_dimension = len(a[0])
	v = [] * module_dimension
	for i in range(module_dimension):
		v[i] = [0] * ring_dimension
		for j in range(module_dimension):
			prod = mul(G[j][i], a[j])
			for k in range(ring_dimension):
				v[i][k] += prod[k]
		v[i] += b[i]
	return v


def embed_msg(msg):
	assert(all(m == 0 or m == 1 for m in msg))
    embedding = [0]*(len(msg)/4)
    for i in range(len(msg)/4):
        for j in range(4):
            embedding[i] += msg[i*4+j] << (14 + 16 * j)

def extract_msg(embedding):
	msg = []
	for e in embedding:
		while e < 0:
			e += (1<<64) - (1<<32) + 1
		if e > (1<<64) - (1<<32) + 1:
			e = e % ((1<<64) - (1<<32) + 1)
		for i in range(4):
			chunk = e & 0xffff
			e >>= 16

			if chunk < (1<<14) or ((1<<16) - chunk) < (1<<14):
				msg += [0]
			else:
				msg += [1]
	return msg

	
