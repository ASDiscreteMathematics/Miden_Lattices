Qmiden = 2**64 - 2**32 + 1

# 512-th primitive root of unity mod Qmiden
# needed to multipy in Z[x]/(x^256+1)

psi_512 = 237214921853999334

# inverse of N=256 modulo Qmiden

inv_N = 18374686475393433601  

# bit reversed arrays of powers of psi_512 

psi_512_rev = [ 1, 18446462594437873665, 1099511627520, 16777216, 68719476736, 18442240469788262401, 18446744069414580225,
1152921504606846976, 262144, 18446744052234715141, 288230376084602880, 4398046511104, 18014398509481984, 64,
18446744068340842497, 70368744161280, 512, 18302628881338728449, 562949953290240, 8589934592, 35184372088832,
16140901060737761281, 18446744069412487169, 137438953440, 134217728, 18446735273321564161, 18446744069414584313,
2251799813685248, 9223372036854775808, 32768, 18446743519658770433, 36028797010575360, 17591917604864, 17592454475776,
18442240469787213841, 18442240469787213809, 18446744065119551490, 4294901759, 18374687574905061377,
18374685375881805825, 4611615648609468416, 4611756386097823744, 18446743794540871745, 18446743794532483137,
18445618152328134657, 1125882726711296, 288230376151712768, 18158513693262873601, 9007061813690368, 9007336691597312,
16140901060200898561, 16140901060200882177, 18446741870357774849, 2198989700608, 562949953421314, 18446181119461163011,
18410715272395620225, 36028797018963840, 18446603334073745409, 18446603329778778113, 17870274521152356353,
576451956076183552, 34360262648, 18446744035055370249, 18035314424752866021, 3328437340319972906, 16105685926854668541,
16933017626115159474, 16329239638270742865, 15113979899245772281, 6562114217670983589, 17311265416183374564,
4299803665592489687, 17330401598553671485, 10382722127243543029, 12053668962110821384, 8340939052496745868,
10561990880479197442, 4644772024090268603, 16192975500896648969, 10708950766175242252, 7059463857684370340,
416595521271101505, 18182056015521604139, 4195631349813649467, 9171943329124577373, 2495058814089251146,
8930739766887302688, 6336932523019185545, 281721071064741919, 3291437157293746400, 10265989416269385394,
9362914843564906265, 2843318466875884251, 16940035449150731648, 8215369291935911999, 5209436881246729393,
12110422903908887252, 9778634991702905054, 4497639551463306333, 12481021517947587610, 13039192753378044028,
5029422726070465669, 17449332314429639298, 10158338780952714962, 8494120110792728509, 14041890976876060974,
5575382163818481237, 18142929134658341675, 1362567150328163374, 7298973816981743824, 17090085178304640863,
10900537202625306992, 2430519478049941168, 7593472940535036657, 15395185741804386692, 7709569171718681254,
16792080670893602455, 10967010099451201909, 5834015391316509212, 17534372342291866343, 14004640413449681173,
13664737158269917819, 13797081185216407910, 10467450029535024137, 15104850399680027611, 10832292272906805046,
6366922389463153702, 237214921853999334, 11917386045977981907, 9770812262840623888, 13183111817796039999,
4406114516091528337, 8187602034452531322, 6045115764991696949, 14918307414590806615, 494216498237666005,
4459017075746761332, 10949047808060940701, 5290167988639048753, 12050543972821699434, 15181754998655957376,
4831070854124318830, 16589430531118646239, 10773575572760153082, 14276112633913910454, 3588235763047079665,
16691665375249202323, 5427855770283221382, 4641337882585395997, 14493012083513256281, 1221351532855077986,
13231174195295398387, 14067222244347930501, 16549024694582589649, 15341376048663650670, 8665994900238946994,
6979306088310177371, 1644572010096941946, 8286160002038086708, 743439328957095187, 16774647971309520093,
10006174791165856646, 2415297291837877958, 5602886161730919912, 3373377623857539246, 17032024114559111334,
5163568085532294797, 16755100833091933884, 1573035775395650770, 5464760906092500108, 8096577031901772269,
14780429931651188987, 10686628914424923326, 11441669947107069577, 13205888274518958430, 11706055037741049324,
10883769032692578351, 13413385849078745835, 700360770216364989, 9432384046970425189, 11622144959503752099,
13533145890581203496, 5862457866249378161, 875634265288439343, 12184322017746068437, 12499229437757822825,
13376768784840513824, 4415056545429189734, 11317759638843783896, 10517142914396393667, 9906467147968854674,
3863141824208378587, 4764679877301742210, 16508747943021518732, 408281368649950045, 12113519742882795430,
5932183857725394514, 3877499600194655066, 526448012694119458, 10094442559346256270, 3200815326405523330,
7721858563281021845, 502012629086366038, 8655137032736432017, 7433403846946633395, 10763480545232365762,
5095456396745892551, 4126998567329314197, 4545880015766881148, 3870163035137971766, 6125875985213995509,
4016101032690928304, 12012107771410162524, 11478179872302871445, 11286965527584982002, 3266250949199600360,
15503969011144524712, 5988353545162139946, 17222793189829815283, 4211584101552955664, 5873491337271928114,
13772306473425142486, 7882761346440596851, 17799792287555502819, 12418052014939319938, 8917938738259842505,
14424608849493817968, 16723337261179401086, 15122929597976639421, 12030096568512274289, 11779090253969091270,
4837070530626986986, 2454730591976115881, 9809941408468046069, 382428689435778886, 16901410098125234092,
13935318169262536835, 17608981172539450419, 17345757166192390690, 802080937612788754, 12362671770314801832,
9638848843637035273, 6702103175001071216, 3059429515486231088, 13754189079328553053, 16643667963227857075,
17255643403020241594, 4716406379463037818, 2443466371579597244, 5175614254872652016, 11336048296972946422,
1999001684679808555, 14439691868389311614, 13787254465881465880, 8143771702088974879 ]
		
# normal version forward NTT, output in bit reversed order		
		
def NTT_256(a):
	k = 0
	len = 128
	N = 256
	
	while len > 0:
		start = 0
		while start < N:
			k = k+1
			zeta = psi_512_rev[k]
			for j in range(start,start+len):
				t = (zeta * a[j + len])%Qmiden
				a[j + len] = (a[j] - t)%Qmiden
				a[j] = (a[j] + t)%Qmiden
			start = start + 2*len
		len >>= 1

# normal version backward NTT, input is assumed bit reversed, output normal ordering

def iNTT_256(a):

	k = 256
	len = 1
	N = 256

	while (len < N):
		start = 0
		while start < N:
			k = k-1
			zeta = -psi_512_rev[k]
			for j in range(start, start + len):
				t = a[j]
				a[j] = (t + a[j + len])%Qmiden
				a[j + len] = t - a[j + len]
				a[j + len] = (zeta * a[j + len])%Qmiden
			start = start + 2*len
		len <<= 1
 
	for j in range(N):
		a[j] = (inv_N*a[j])%Qmiden


# 4 element version
# a is now array of Quadruples
# the original a has to be ordered as [[a[i], a[i+128], a[i+1], a[i+1+128]]

def NTT_4_256(a):
	k = 0
	len = 128
	N = 256

	while len > 2:
		start = 0
		while start < (N>>2):
			k = k+1
			zeta = psi_512_rev[k]
			for j in range(start,start+(len>>2)):
				A1 = a[j]  # read 4 tuple
				A2 = a[j + (len >> 2)] # second 4 tuple	
				t1 = (zeta * A1[1])%Qmiden
				t2 = (zeta * A2[1])%Qmiden
				t3 = (zeta * A1[3])%Qmiden
				t4 = (zeta * A2[3])%Qmiden
				a[j] = [(A1[0] + t1)%Qmiden, (A2[0] + t2)%Qmiden, (A1[2] + t3)%Qmiden, (A2[2] + t4)%Qmiden]
				a[j+(len>>2)] = [(A1[0] - t1)%Qmiden, (A2[0] - t2)%Qmiden, (A1[2] - t3)%Qmiden, (A2[2] - t4)%Qmiden]
			start += (len>>1)
		len >>= 1	

  # doing len = 2 
  
	for j in range(0,N>>2):
		k = k+1
		zeta = psi_512_rev[k]
		A1 = a[j]
		t1 = (zeta * A1[1])%Qmiden
		t2 = (zeta * A1[3])%Qmiden
		a[j] = [(A1[0] + t1)%Qmiden, (A1[2] + t2)%Qmiden, (A1[0] - t1)%Qmiden, (A1[2] - t2)%Qmiden]	
   
  # doing len = 1

	for j in range(0,N>>2):
		k = k+1
		zeta = psi_512_rev[k]
		A1 = a[j]
		t1 = (zeta * A1[1])%Qmiden
		k = k+1
		zeta = psi_512_rev[k]  
		t2 = (zeta * A1[3])%Qmiden
		a[j] = [(A1[0] + t1)%Qmiden, (A1[0] - t1)%Qmiden, (A1[2] + t2)%Qmiden, (A1[2] - t2)%Qmiden]


# 4 element version
# a is now array of Quadruples where input is assumed in bit reversed order
# chopped up into 4 element arrays, i.e. [[a[i], a[i+1], a[i+2], a[i+3]]
# for i = 4*k and a the original array in bit reversed order

# the output is eQual to map_4_lin_offset(original)

def iNTT_4_256(a):
	
	k = 256
	len = 1
	N = 256
   
	# len = 1 separate

	for j in range(N >> 2):
		k = k-1
		zeta1 = -psi_512_rev[k]
		k = k-1
		zeta2 = -psi_512_rev[k]
		
		A1 = a[j]
		u1 = (A1[0] + A1[1])%Qmiden
		u2 = (zeta1*(A1[0] - A1[1]))%Qmiden
			
		v1 = (A1[2] + A1[3])%Qmiden
		v2 = (zeta2*(A1[2] - A1[3]))%Qmiden
		
		a[j] = [u1, v1, u2, v2]
  
	len = 2  
   
	while (len < (N >> 1)):
		start = 0
		while (start < (N>>2)):  #for start := 0 to (N div 4)-1 by len do
			k = k-1
			zeta1 = -psi_512_rev[k]
			k = k-1
			zeta2 = -psi_512_rev[k]
			for j in range(start, start + (len>>1)):
				A1 = a[j]
				A2 = a[j + (len>>1)];
				
				# first pairs
		
				u1 = (A1[0] + A1[1])%Qmiden
				u2 = (zeta1*(A1[0] - A1[1]))%Qmiden
        
				v1 = (A2[0] + A2[1])%Qmiden
				v2 = (zeta2*(A2[0] - A2[1]))%Qmiden
		
				# second pairs
				
				u3 = (A1[2] + A1[3])%Qmiden
				u4 = (zeta1*(A1[2] - A1[3]))%Qmiden
        
				v3 = (A2[2] + A2[3])%Qmiden
				v4 = (zeta2*(A2[2] - A2[3]))%Qmiden
				
				a[j] = [u1, v1, u3, v3]
				a[j+(len>>1)] = [u2, v2, u4, v4]
			start = start + len
		len = 2*len
  
	# len = N/2 and scaling by N together
 
	k = k-1
	zeta = -psi_512_rev[k]
	zeta_invN = (zeta*inv_N)%Qmiden
 
	for j in range(N>>2):
		A = a[j]
		a[j] = [((A[0] + A[1])*inv_N)%Qmiden, (zeta_invN*(A[0] - A[1]))%Qmiden, ((A[2] + A[3])*inv_N)%Qmiden, (zeta_invN*(A[2] - A[3]))%Qmiden];
  
  
# maps a linear list to 4 element tuple list interleaving offset	
	
def map_lin_4_offset(a):
	res = []
	for i in range(N>>2):
		res.append([a[2*i], a[2*i+128], a[2*i+1], a[2*i+1+128]])
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
