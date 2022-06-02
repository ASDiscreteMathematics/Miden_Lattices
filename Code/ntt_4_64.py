import random

q = 2**64 - 2**32 + 1

# 128-th primitive root of unity mod q
# needed to multipy in Z[x]/(x^256+1)

psi_128 = 4294901759

# inverse of N=64 modulo q

N = 64
inv_N = 18158513693329981441  

# bit reversed arrays of powers of psi_128 

psi_128_rev = [1, 18446462594437873665, 1099511627520, 16777216, 18446744000695107585, 4503599626321920, 4096, 17293822564807737345, 
18446744052234715141, 18446744069414322177, 4398046511104, 18158513693329981441, 18446744069414584257,
18014398509481984, 18446673700670423041, 18446744068340842497, 18446181119461294081, 18446744060824649729,
144115188075855872, 512, 18446744069412487169, 137438953440, 16140901060737761281, 18446708885042495489,
18444492269600899073, 18446744069414584313, 134217728, 18446735273321564161, 36028797010575360, 549755813888,
9223372032559808513, 18446744069414551553, 4294901759, 4295032831, 18374685375881805825, 72056494509522944,
4503599627370512, 18442240469787213841, 17591917604864, 17592454475776, 1125917086449664, 18445618186687873025,
18158513693262871553, 288230376151710720, 18446743794540871745, 18446743794532483137, 4611756386097823744,
13835128420805115905, 562949953421310, 562949953421314, 18446741870357774849, 2198989700608, 18437737007600893953,
18437736732722987009, 2305843009213685760, 2305843009213702144, 34360262648, 18446744035055370249, 576451956076183552,
576469548262227968, 18410715272395620481, 18410715272395620225, 140739635806208, 18446603334073745409]

		
# normal version forward NTT, output in bit reversed order		
		
def NTT_64(a):
	k = 0
	len = 32
	N = 64
	
	while len > 0:
		start = 0
		while start < N:
			k = k+1
			zeta = psi_128_rev[k]
			for j in range(start,start+len):
				t = (zeta * a[j + len])%q
				a[j + len] = (a[j] - t)%q
				a[j] = (a[j] + t)%q
			start = start + 2*len
		len >>= 1

# normal version backward NTT, input is assumed bit reversed, output normal ordering

def iNTT_64(a):

	k = 64
	len = 1
	N = 64
	
	while (len < N):
		start = 0
		while start < N:
			k = k-1
			zeta = -psi_128_rev[k]
			for j in range(start, start + len):
				t = a[j]
				a[j] = (t + a[j + len])%q
				a[j + len] = t - a[j + len]
				a[j + len] = (zeta * a[j + len])%q
			start = start + 2*len
		len <<= 1
 
	for j in range(N):
		a[j] = (inv_N*a[j])%q


# 4 element version
# a is now array of quadruples
# the original a has to be ordered as [[a[i], a[i+32], a[i+1], a[i+1+32]]

def NTT_4_64(a):
	k = 0
	len = 32
	N = 64

	while len > 2:
		start = 0
		while start < (N>>2):
			k = k+1
			zeta = psi_128_rev[k]
			for j in range(start,start+(len>>2)):
				A1 = a[j]  # read 4 tuple
				A2 = a[j + (len >> 2)] # second 4 tuple	
				t1 = (zeta * A1[1])%q
				t2 = (zeta * A2[1])%q
				t3 = (zeta * A1[3])%q
				t4 = (zeta * A2[3])%q
				a[j] = [(A1[0] + t1)%q, (A2[0] + t2)%q, (A1[2] + t3)%q, (A2[2] + t4)%q]
				a[j+(len>>2)] = [(A1[0] - t1)%q, (A2[0] - t2)%q, (A1[2] - t3)%q, (A2[2] - t4)%q]
			start += (len>>1)
		len >>= 1	

  # doing len = 2 
  
	for j in range(0,N>>2):
		k = k+1
		zeta = psi_128_rev[k]
		A1 = a[j]
		t1 = (zeta * A1[1])%q
		t2 = (zeta * A1[3])%q
		a[j] = [(A1[0] + t1)%q, (A1[2] + t2)%q, (A1[0] - t1)%q, (A1[2] - t2)%q]	
   
  # doing len = 1

	for j in range(0,N>>2):
		k = k+1
		zeta = psi_128_rev[k]
		A1 = a[j]
		t1 = (zeta * A1[1])%q
		k = k+1
		zeta = psi_128_rev[k]  
		t2 = (zeta * A1[3])%q
		a[j] = [(A1[0] + t1)%q, (A1[0] - t1)%q, (A1[2] + t2)%q, (A1[2] - t2)%q]


# 4 element version
# a is now array of quadruples where input is assumed in bit reversed order
# chopped up into 4 element arrays, i.e. [[a[i], a[i+1], a[i+2], a[i+3]]
# for i = 4*k and a the original array in bit reversed order

# the output is equal to map_4_lin_offset(original)

def iNTT_4_64(a):
	
	k = 64
	len = 1
	N = 64
   
	# len = 1 separate

	for j in range(N >> 2):
		k = k-1
		zeta1 = -psi_128_rev[k]
		k = k-1
		zeta2 = -psi_128_rev[k]
		
		A1 = a[j]
		u1 = (A1[0] + A1[1])%q
		u2 = (zeta1*(A1[0] - A1[1]))%q
			
		v1 = (A1[2] + A1[3])%q
		v2 = (zeta2*(A1[2] - A1[3]))%q
		
		a[j] = [u1, v1, u2, v2]
  
	len = 2  
   
	while (len < (N >> 1)):
		start = 0
		while (start < (N>>2)):  #for start := 0 to (N div 4)-1 by len do
			k = k-1
			zeta1 = -psi_128_rev[k]
			k = k-1
			zeta2 = -psi_128_rev[k]
			for j in range(start, start + (len>>1)):
				A1 = a[j]
				A2 = a[j + (len>>1)];
				
				# first pairs
		
				u1 = (A1[0] + A1[1])%q
				u2 = (zeta1*(A1[0] - A1[1]))%q
        
				v1 = (A2[0] + A2[1])%q
				v2 = (zeta2*(A2[0] - A2[1]))%q
		
				# second pairs
				
				u3 = (A1[2] + A1[3])%q
				u4 = (zeta1*(A1[2] - A1[3]))%q
        
				v3 = (A2[2] + A2[3])%q
				v4 = (zeta2*(A2[2] - A2[3]))%q
				
				a[j] = [u1, v1, u3, v3]
				a[j+(len>>1)] = [u2, v2, u4, v4]
			start = start + len
		len = 2*len
  
	# len = N/2 and scaling by N together
 
	k = k-1
	zeta = -psi_128_rev[k]
	zeta_invN = (zeta*inv_N)%q
 
	for j in range(N>>2):
		A = a[j]
		a[j] = [((A[0] + A[1])*inv_N)%q, (zeta_invN*(A[0] - A[1]))%q, ((A[2] + A[3])*inv_N)%q, (zeta_invN*(A[2] - A[3]))%q];
  
  
# maps a linear list to 4 element tuple list interleaving offset	
	
def map_lin_4_offset(a):
	res = []
	for i in range(N>>2):
		res.append([a[2*i], a[2*i+32], a[2*i+1], a[2*i+1+32]])
	return res	
	
# maps list of 4 tuples to linear list	
	
def map_4_lin(a):
	res = []
	for i in range(N>>2):
		res.extend(a[i])
	return res	
	
# maps list of 4 tuples to linear list with offset, inverse operation of map_lin_4_offset
	
def map_4_lin_offset(a):
	res = []
	for i in range(N>>2):
		res.append(a[i][0])
		res.append(a[i][2])

	for i in range(N>>2):
		res.append(a[i][1])
		res.append(a[i][3])
		
	return res		
		
		
# a, b given as arrays of coefficients, note in place a := a*b, b also gets changed
def FastMul_64(a,b):
	NTT_64(a)
	NTT_64(b)

	for j in range(N):
		a[j] = (a[j]*b[j])%q;

	iNTT_64(a);


# schoolbook multiply
def schoolbook_mult(a, b):

	c = [0 for _ in range(N)]
	
	for j in range(N):
		for k in range(N):
			i = (j+k)%N
			prod = (a[j]*b[k])%q;
			if i != (j+k):
				prod = -prod
			c[i] = (c[i] + prod)%q;
	return c
		


# sanity checking

a = [random.randint(0,q-1) for _ in range(N)]
b = [random.randint(0,q-1) for _ in range(N)]

c = schoolbook_mult(a, b)
FastMul_64(a,b)

print("Sanity check fastmult:", a == c)

		
# sanity checking packing 4 elements

a = [random.randint(0,q-1) for _ in range(N)]
a4 = map_lin_4_offset(a)

a_orig = a.copy() 
a4_orig = a4.copy()

NTT_64(a)
NTT_4_64(a4)

b = map_4_lin(a4)
print("Sanity check 4 way forward:", a == b)

iNTT_64(a);
print("Sanity check inverse NTT:", a == a_orig)

iNTT_4_64(a4)
print("Sanity check 4 way inverse:", a4_orig == a4)

c = map_4_lin_offset(a4)
print("Sanity check inverse maps:", c == a)