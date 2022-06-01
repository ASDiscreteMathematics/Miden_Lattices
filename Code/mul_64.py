import random

def bit_reverse_index(width, index):
    reversed_index = 0
    for i in range(width):
        if index & (1 << i) != 0:
            reversed_index |= 1 << (width-1 - i)
    return reversed_index

assert(set([bit_reverse_index(7, i) for i in range(128)]) == set(list(range(128))))
assert(all(bit_reverse_index(7, bit_reverse_index(7, i)) == i for i in range(128)))

# 128-th primitive root of unity mod q
# needed to multipy in Z[x]/(x^64+1)

psi_128 = 17870292113338400769

psi_128_array = [0]*128
acc = 1
for i in range(len(psi_128_array)):
    psi_128_array[i] = acc
    acc = (acc * psi_128) % ((1<<64) - (1<<32) + 1)
psi_128_inv_array = [0] * 128
acc = 1
for i in range(len(psi_128_inv_array)):
    psi_128_inv_array[i] = acc
    acc = (acc * psi_128_array[-1]) % ((1<<64) - (1<<32) + 1)

# bit reversed arrays of powers of psi_128 and its inverse 
psi_128_rev = [psi_128_array[bit_reverse_index(7, i)] for i in range(128)]
psi_128_inv_rev = [psi_128_inv_array[bit_reverse_index(7,i)] for i in range(128)]

assert(all((f*i) % ((1<<64) - (1<<32) + 1) == 1 for (f, i) in zip(psi_128_rev, psi_128_inv_rev)))

#psi_128_rev = [ 1, 18446462594437873665, 1099511627520, 16777216, 68719476736, 18442240469788262401, 18446744069414580225,
#999001684679808555, 14439691868389311614, 13787254465881465880, 8143771702088974879 ]

#psi_128_inv_rev = [ 1, 281474976710656, 18446744069397807105, 18446742969902956801, 17293822564807737345, 4096, 4503599626321920,
#259142034962052999, 14040629553323055984, 5263632251618544322, 8675931806573960433, 6529358023436602414,
#18209529147560584987 ]

def NTT_sb_64(a):
	t = 64
	m = 1
	while m < 64:
		t >>= 1
		for i in range(m):
			j1 = 2 * i * t
			j2 = j1 + t - 1
			S = psi_128_rev[m + i]
			for j in range(j1,j2+1):
				U = a[j]
				V = (a[j + t] * S) % ((1<<64) - (1<<32) + 1)
				a[j] = (U + V) % ((1<<64) - (1<<32) + 1)
				a[j + t] = (U - V) % ((1<<64) - (1<<32) + 1)
		m = 2*m
		
		
# second alternative version		
		
def NTT_64(a):
	k = 0
	m = 32
	N = 64
	while m > 0:
		start = 0
		while start < N:
			k = k+1;
			zeta = psi_128_rev[k]
			for j in range(start,start+m):
				t = (zeta * a[j + m]) % ((1<<64) - (1<<32) + 1)
				a[j + m] = (a[j] - t) % ((1<<64) - (1<<32) + 1)
				a[j] = (a[j] + t) % ((1<<64) - (1<<32) + 1)
			start = start + 2*m
		m >>= 1

def INTT_sb_64(a):
	t = 1
	m = 64
	while m > 1:
		j1 = 0
		h = m >> 1
		for i in range(h):
			j2 = j1 + t - 1
			S = psi_128_inv_rev[h + i]
			for j in range(j1,j2+1):
				U = a[j]
				V = a[j + t]
				a[j] = (U + V) % ((1<<64) - (1<<32) + 1)
				a[j + t] = ((U - V )*S) % ((1<<64) - (1<<32) + 1)
			j1 = j1 + 2*t
		t <<= 1
		m >>= 1

	#scaling by inverse of N mod q for N = 64
	invN = 18158513693329981441
	for j in range(64): 
		a[j] = (a[j] * invN) % ((1<<64) - (1<<32) + 1)


# a, b given as arrays of coefficients, note in place a := a*b, b also gets changed
def FastMul_64(a,b):
	NTT_sb_64(a)
	NTT_sb_64(b)

	for j in range(64):
		a[j] = (a[j]*b[j])  % ((1<<64) - (1<<32) + 1)

	INTT_sb_64(a)

# schoolbook multiply
def schoolbook_mult(a, b):

	c = [0 for _ in range(64)]
	
	for j in range(64):
		for k in range(64):
			i = (j+k) % 64
			prod = (a[j]*b[k]) % ((1<<64) - (1<<32) + 1)
			if i != (j+k):
				prod = -prod
			c[i] = (c[i] + prod) % ((1<<64) - (1<<32) + 1)
	return c
		


# sanity checking

N = 64
a = [random.randint(0,((1<<64)-(1<<32))) for _ in range(N)]
b = [random.randint(0,((1<<64)-(1<<32))) for _ in range(N)]

a = list(range(5)) + [0]*(N-5)
b = list(range(5)) + [0]*(N-5)

c = a.copy()
d = a.copy()

NTT_64(a)
NTT_sb_64(c)

print("Sanity check :", a == c)

INTT_sb_64(a)

print("Sanity check:", a == d)

c = schoolbook_mult(a, b)
FastMul_64(a,b) # sets a

print("Sanity check :", a == c)

print(a[:15])
print(c[:15])
