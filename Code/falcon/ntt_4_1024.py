import random 

four_elem = False

q = 2**64 - 2**32 + 1
N = 1024
hN = N >> 1

# 2048-th primitive root of unity mod q
# needed to multipy in Z[x]/(x^1024+1)

psi_2048 = 19582517437079335

# inverse of N modulo q

inv_N = 18428729670909296641  

# bit reversed arrays of powers of psi_2048

psi_2048_rev = [ 1, 18446462594437873665, 18446742969902956801, 18446744069397807105, 18442240469788262401, 18446744000695107585,
17293822564807737345, 18446744069414580225, 18158513693329981441, 18446739671368073217, 18446744052234715141,
18446744069414322177, 18446673700670423041, 18446744068340842497, 18428729670905102337, 18446744069414584257,
16140901060737761281, 18446708885042495489, 18446743931975630881, 18446744069412487169, 18446181119461294081,
18446744060824649729, 18302628881338728449, 18446744069414583809, 18410715272404008961, 18446743519658770433,
9223372032559808513, 18446744069414551553, 18446735273321564161, 18446744069280366593, 18444492269600899073,
18446744069414584313, 18446743794540871745, 18446743794532483137, 13834987683316760577, 4611615648609468416,
18158513693262873601, 18158513693262871553, 18445618152328134657, 1125882726711296, 18446744065119682562,
18446744065119551490, 18374685375881805825, 72056494509522944, 18442240469787213841, 18442240469787213809,
18446726476960108545, 17591917604864, 18446744035055370249, 18446744035054321673, 17870274521152356353,
576451956076183552, 18410715272395620481, 18410715272395620225, 18446603329778778113, 140735340838912,
16140901060200898561, 16140901060200882177, 18437736732722987009, 9007061813690368, 18446181119461163011,
18446181119461163007, 18446741870357774849, 2198989700608, 4195631349813649467, 9171943329124577373,
15951685255325333175, 9516004302527281633, 18182056015521604139, 18030148548143482816, 10708950766175242252,
7059463857684370340, 1506708620263852673, 10231374777478672322, 2843318466875884251, 9083829225849678056,
6336932523019185545, 281721071064741919, 15155306912120837921, 8180754653145198927, 12053668962110821384,
8064021942171041292, 4299803665592489687, 17330401598553671485, 13801972045324315718, 2253768568517935352,
10561990880479197442, 10105805016917838453, 16329239638270742865, 15113979899245772281, 11884629851743600732,
1135478653231209757, 16933017626115159474, 2341058142559915780, 18035314424752866021, 3328437340319972906,
10832292272906805046, 6366922389463153702, 3341893669734556710, 10467450029535024137, 912371727122717978,
4442103655964903148, 13664737158269917819, 13797081185216407910, 3051558327610197629, 7593472940535036657,
7546206866789277329, 16016224591364643153, 10967010099451201909, 5834015391316509212, 1654663398520981866,
7709569171718681254, 5965722551466996711, 5407551316036540293, 5029422726070465669, 17449332314429639298,
13949104517951277988, 9778634991702905054, 13237307188167854928, 6336321165505697069, 7298973816981743824,
17090085178304640863, 17084176919086420947, 18142929134658341675, 8288405288461869359, 9952623958621855812,
14041890976876060974, 5575382163818481237, 1691643236322650437, 16873708294018933551, 5464760906092500108,
8096577031901772269, 7760115154989660995, 14780429931651188987, 13205888274518958430, 7005074122307514744,
17032024114559111334, 5163568085532294797, 15073366445557045075, 5602886161730919912, 17703304740457489134,
1672096098105064228, 10006174791165856646, 2415297291837877958, 7128984430570800425, 4415056545429189734,
9906467147968854674, 7929601155018190654, 12499229437757822825, 13376768784840513824, 6262422051668515884,
875634265288439343, 6740689031673534997, 7562975036722005970, 13413385849078745835, 700360770216364989,
6824599109910832222, 9432384046970425189, 5862457866249378161, 4913598178833380825, 3953731985901328040,
17225392536559506335, 4641337882585395997, 13018888299131362939, 10773575572760153082, 14276112633913910454,
14858508306367504656, 1755078694165381998, 6979306088310177371, 9780749169175637327, 10160584067376497613,
1644572010096941946, 1897719374831994672, 3105368020750933651, 14067222244347930501, 5215569874119185934,
494216498237666005, 4459017075746761332, 7497696261353643620, 13156576080775535568, 15181754998655957376,
6396200096592884887, 1857313538295938082, 4831070854124318830, 12401628304422887372, 3528436654823777706,
8187602034452531322, 14040629553323055984, 237214921853999334, 11917386045977981907, 8675931806573960433,
5263632251618544322, 16447742384734775766, 4007052201025272707, 13787254465881465880, 8143771702088974879,
7110695772441637899, 5175614254872652016, 13730337689951546503, 16003277697834987077, 12362671770314801832,
17644663131801795567, 11744640894413513105, 9638848843637035273, 15387314553928353233, 4692554990086031268,
16643667963227857075, 17255643403020241594, 6667653815445493051, 12030096568512274289, 1723406808235183235,
3323814471437944900, 12418052014939319938, 646951781859081502, 4022135219920766353, 8917938738259842505,
1545333971289350229, 4511425900152047486, 17608981172539450419, 17345757166192390690, 18064315379978805435,
9809941408468046069, 13609673538787597335, 15992013477438468440, 15245928743009060991, 10094442559346256270,
502012629086366038, 10724885506133562476, 8655137032736432017, 7433403846946633395, 7683263524182218559,
13351287672668691770, 526448012694119458, 14569244469219929255, 12113519742882795430, 5932183857725394514,
13682064192112842111, 3863141824208378587, 408281368649950045, 1937996126393065589, 4211584101552955664,
5873491337271928114, 4674437595989441835, 10563982722973987470, 17222793189829815283, 12458390524252444375,
3266250949199600360, 15503969011144524712, 13900864053647703173, 4126998567329314197, 6125875985213995509,
14576581034276612555, 4016101032690928304, 12012107771410162524, 6968564197111712876, 7159778541829602319,
3929858786378642316, 7112675246154377750, 8362186383759854260, 4689912295553665450, 14358435046894549184,
6820186327231405039, 6108922642702621141, 7305983592013185897, 15049383599936516047, 12216811346274483113,
3589423675261482283, 6414348153479289383, 224350547607727331, 5006481804801239664, 12489358087930152296,
6743454643571072270, 9714604383004622450, 5500770423122943299, 10268645332677273943, 14421297089005146422,
1794804380861818648, 3158366299580748670, 7681144356368296763, 17054149009739409518, 4187015958668887546,
17668002479022071670, 11977893002791800486, 3107636527861734213, 11557258861835081117, 625810225600154958,
1561169760991269037, 5454617847800030114, 10141440556758839363, 12486396704535672088, 11984005542840547177,
10019076149651226019, 18147478832086943132, 12577098593386243920, 9906341360885134636, 2727123837437860044,
17740512949860132546, 11724314991892485077, 6816548736552749790, 8515228971291783927, 10659847895797062167,
14031575217582598302, 5919394105455887829, 15030590866359316324, 12796895112978970121, 1560799588066959011,
17638901753592829678, 12781599562090518453, 11491806888718160052, 1572137324173280490, 10461664704817933990,
9564262514387024666, 16052622170793454809, 8383068400017029755, 5463754609422739804, 3370246630088296031,
3638323995651455811, 6365632919551470868, 7657453289212455099, 11102195893002206701, 3885345703086309931,
13294663425500451833, 17527399748276289503, 12316227567088954246, 15948194847534211277, 13869829016883058002,
152897937997792376, 5164132063260791647, 6113546424387384073, 2225341748615664720, 9785468031858712064,
16909802868642731951, 14948939724807468932, 13475313378280530262, 2308232038958038546, 9592291974280344910,
12014883256269903942, 17802733988925317760, 4496767977211359228, 6151214463239765361, 8911053381972245530,
15568786679171320491, 19112242249724047, 2951359516584421996, 16905094363786184290, 278167718576958090,
1223183503982339008, 4419568367257164534, 11091989500308225777, 6296100189638712363, 14123587056930693059,
5810722514138689194, 15137851097357772055, 2849220811487664857, 4295002282146690441, 11575701465310827636,
12582249520745188423, 5940943517668292086, 12066191983457963400, 13315466594398149922, 12053974342864933269,
11285503742479007084, 15919780095311700439, 3639634848410716242, 16625729085584007730, 2975131003309595864,
16329435310479291959, 8854965448075557493, 4198074395846544547, 16497053662173719388, 16677776346006097586,
10670334717871145615, 3878624198769971593, 5354303957062182591, 1508273997932245425, 15499491376360706981,
8424275818888585779, 10634060002517168046, 4295815520590785595, 14290012408112277771, 15913274187758939207,
371891375413699483, 4347022422486734535, 8024399707039913807, 9233735841019365682, 11888040935601451494,
4573069520345433394, 3192518480993723602, 7847524879092096009, 16260897004766174524, 12458073038457296305,
12890081553990608899, 4179502387700367909, 7679740417818447560, 4106679476439837717, 13308480401157259412,
15975288260888972401, 1406998020037882997, 4518113032494938455, 17783460465441878945, 14989275032188358951,
6097691134303827517, 14406691742104117415, 14234122862185153691, 17121841670624273282, 11255984160303063976,
17698160190544923319, 13140475237632941313, 7439966824493015109, 959967552227305945, 7430863960585448835,
10886932084851949587, 18137812093348882831, 7093403778535204495, 2870607137738690347, 18363833618917996149,
345137759837927448, 1941397000334789249, 5280444194155204294, 663576273645521453, 9115369905823964012,
6337038648111976301, 1415419453610921553, 6715029048772165709, 11534607820881582817, 18188848021460212523,
16799868753440642108, 5486745524883165993, 5907035176470557038, 5575393374484204350, 13568943604939006010,
14804671509201811970, 43142219979740931, 16383575685779609937, 5271741541623046617, 7000476060236159302,
10362793272935287662, 7709658857044466158, 16317828492439126475, 7756907657126989834, 17582727038347959133,
13802821046066641766, 11323355628887372424, 16826744251348157030, 5350065414412465710, 5308610189164171624,
15531176002678313992, 15685641990711164737, 9546597805441465650, 17198755642670596954, 2540820207615693247,
16121944306381782101, 3235915244053982668, 14586854504982200837, 16434255353882186006, 4232816070675458120,
4184390855894463221, 11221484848131637518, 327930691828598087, 12645811551425139186, 15038540732087693240,
17233511790631916809, 12362461035457730117, 16207038811842065314, 15028382777741121447, 15984902507394762860,
2623445534628784696, 8932772064328191883, 9627861440214039994, 8740885839153244225, 6665967936588919331,
529102008834432265, 7440577883017277023, 6014371623370100770, 2346834345155397801, 15415784495989080639,
1879817591510961655, 18295090034566750882, 8462836655462685385, 15860937903541196405, 5828091010015769947,
4665691128499368837, 6686338362028015651, 1913039459307282826, 12273731855508974132, 14390139890846703192,
18290750489319984117, 16672792867292992783, 10755587837161802966, 17078493612874372559, 8463154943360171265,
15594331550120231815, 3650541573257562281, 11754060979178594938, 3456327113326256432, 14383800816696994133,
12257726419636086444, 7500740417092890225, 12365007338637617157, 14074187984474348594, 10757588516645913927,
1798767486355837899, 9203872837195467135, 4389942117088447138, 5956134496998871451, 4440654710286119610,
17198795428657782689, 4255134452441852017, 16597218757394956566, 15304315674458262608, 432040889165782054,
8715504128117593387, 19582517437079335, 5667459462650274795, 6548472451862452516, 14833754696748909971,
977285899799902802, 13931017301309072163, 7909919691978484502, 12023728924795965445, 7206065378950026365,
6141391951880571024, 8172770412429231461, 13202144340943831319, 13273867392076105962, 8576353051786437773,
12229268290740484781, 17193462953441506881, 2308290823356457957, 12237647476215399550, 10041931091190098725,
13383434380477728947, 13957218789535926091, 13270592206047749221, 5600425978850956643, 8420495141629964801,
7818287198399222416, 767673993985071378, 7939125327584123053, 3956111051294801955, 15494291476070451486,
7989573157503773842, 8446187562373029718, 18290083929917949641, 11763737490839323877, 9957242804506486702,
10645132529729130693, 1411969236303298865, 12806605743856568005, 9595098171547258567, 17562133852338181982,
17073490902438230581, 7966028552578642196, 5343728688343265695, 17171922384768587588, 4345517660586282245,
17205695403739328796, 16579054845752790076, 10074241128319282014, 3437307432281229033, 8387996212385384605,
5856341367916956918, 8248170592246610457, 16317397215275673639, 8518354744012540121, 3505230280120230361,
6806952748895918828, 9051715388835247943, 10219125603779622435, 2973809094719731252, 11369862332803365609,
7460718733603754401, 11374083960174708260, 11295753890426390920, 5870966158393556332, 16570564489772914910,
16678858123743380913, 375666990542046292, 2538429194941679965, 16626071583888734006, 11889071400710163917,
13456847399023759765, 7652235498811113389, 10137158258732105336, 4584062799452533527, 12688006344449746194,
10127726119132064550, 3142086129344290269, 14885515920950843192, 12604169412247670087, 5595943325276378367,
2464236106469512186, 18225758326205683895, 9270330408525047947, 7234832675398179116, 6689944965339737831,
8403662951119239610, 8599634950908439091, 7874058463381858294, 1267144782341513167, 2878850858608389731,
15421058845117156515, 5877651782245154149, 7310289792198505404, 1860689490118855399, 3881364185207781801,
3005335924336370336, 14143087565369627264, 1089775420116882570, 9402811212061249573, 18098689330632295042,
535888234053312425, 13213827362080690254, 18282406491446742922, 15559851734640773081, 391943999578398962,
15581468049507881811, 7929139079472734785, 18153075338036508171, 6637671903602949247, 14617984856762654786,
15850102909997410879, 11484107350653274400, 4041349390177852804, 13971279979575548562, 8092880427538125317,
16097394218389975121, 16207887089994425334, 6263414437613732362, 16120358863491781106, 18085882527567857916,
13884051052008238111, 13476898549572600427, 17132043445671853129, 13798349460638678722, 3135551996627191696,
15662306159156270089, 4287105872426499400, 1435513418831659300, 9728540708479523761, 8996233907908426289,
4232190452164704782, 14647341494806332478, 17025280296452151564, 6699048990509617996, 6855511880884072764,
13559410886336109653, 8020563897425393614, 4464021796080112361, 14477646780045217513, 805325463025554705,
15253999561031414629, 15092652196876062542, 1260038877477225157, 12605771966736925554, 14536840115127415776,
17265430299226314567, 5140709823874234178, 6442603704204437640, 11351532071763811106, 10060753158520994410,
10080311019817801256, 8612455386822482827, 5614256504531820282, 16698903785247775326, 17950606908243413470,
16241566743615955619, 8824278971159395949, 6498267541963153898, 7075033885715122265, 15410779547903053935,
1817105014390926972, 8628887334295837389, 17141527482521913171, 7918792107371928848, 12644047188125007864,
5996397764581325090, 10006679713294021897, 3388640217848422690, 39115722593613692, 14834575544913119340,
13238203290721534494, 13958789178738624629, 2503406245991276288, 8740604997609669605, 16009025055173377493,
8699858785941968005, 1153532687503936734, 7996139942817448794, 13671905978699354347, 989849013421491106,
1580505898515625983, 14584607772633603877, 17391736024899514018, 14258638079291991077, 9228261500031493872,
11077693977821432078, 6266461428693837892, 8662377673372797199, 312925780748909536, 8010104650731677821,
8918657157927141307, 8005011374273215121, 4755877603291638172, 13771851712143767403, 6226945910957093071,
5821067165119911282, 15334529499874140019, 9897770647513668340, 17516019197232489140, 12170147349297687194,
638751328927570930, 6268023080778906846, 14220584458004245700, 4126179439639438934, 3986596912535370878,
3613417179382635628, 3732452312971992203, 11142912843547685803, 4045205754698933616, 13250696507402086126,
3084211247546459674, 14562691447700927151, 13446031230868382703, 10460593365646500703, 11412874434361353303,
15356326470723149140, 13914901968176884607, 5395188902451009436, 11000945091957822873, 5127458447308575947,
5110010631420567440, 9675049182130121614, 11995771582505614226, 12922079148827575926, 505650719337366702,
2233320433665088896, 677213148894866745, 14223071440482507279, 16560517619117915738, 2882849005092529066,
3202224669556314531, 6853447461186539370, 1904562137429705200, 34895631776017014, 2028194088043702453,
14345523919403080297, 11211512379013626874, 6386112789565834127, 8408715666915300935, 6448153390442303038,
4641444800750985224, 279165054208136112, 16225552704349619624, 4083726938737136450, 15905122754450677708,
14195414177697504374, 11929493127078654517, 14691738984709255662, 238070267178713150, 4616047971325648207,
7171053287035931927, 17934091550663146318, 15236497099437641600, 3104107107372552306, 3356932467041235657,
5417705191158933960, 580180600093873153, 17123959646811123237, 10036420970426885618, 12903336218102904697,
7825583292067718279, 2143390476713559247, 11594715736241105207, 9834508212829624140, 13228982645303394011,
8049782023765701561, 4192044342847360408, 2219227261000078046, 16551406175760922259, 14156778904344190484,
2775240818140193189, 15141643747224812638, 10871226768963172092, 9058023981881859525, 15089610673364298943,
17753818088000624368, 3284040920185287825, 2573766818266017946, 3755182475706961191, 10452685561310995178,
13182837874047039452, 17147123813708473976, 524005542855920051, 4889089424978655836, 13598140815354230483,
10992969397750315971, 7264434128297993269, 6504391485756747660, 10582275380827688672, 10492810266216088765,
5182229877542655744, 11401781171663148294, 12063173382964878836, 12986439443942495773, 10377558731886900618,
12657898587945958274, 2362831314881886490, 1028641288663435027, 80972341836603996, 16895514643714203733,
3647251597124060792, 10290976279272702297, 15726589663754288343, 18068062982682034159, 10989673530509938917,
8229130309307480216, 647778734692831968, 6036908663811539617, 10731268707577902015, 8540833956523281092,
15132252893546800818, 15417295375554183025, 14130411966421174052, 11657795204467044579, 9233493577436867660,
9029468356494744587, 455906449640507599, 17427273095646849068, 4271666716646109083, 4564350881512077310,
8291238217344211485, 1320446646352403721, 17851261006466688773, 5702391734953591701, 14220392861548041004,
3418554491047856910, 10358133758933250987, 14318335397559550085, 14793292948442149158, 15873302663502414709,
17284518142217611933, 12483006042496573711, 5986801160153745741, 14464933718152766765, 6214683737759992527,
17229316179578437891, 7725134980519083461, 16305956891531811746, 9148936651838805217, 7630327992899668083,
11000921142400797286, 5039005328734628194, 12823981763250771574, 8707320950725412881, 6460847635908914725,
8901691858968270959, 9078093793807670612, 3866218763988894754, 7665879171049687338, 8725645740799564966,
3082678475896822106, 13682879565831419937, 7883170898595354553, 11705557770412325915, 369010131771762096,
14016298947912065284, 4132703590028571263, 4492946629632996312, 11889593103854070891, 17279230114556215215,
15569952869032286760, 10847423255292999153, 4617451800662579863, 17512827236137298821, 352803869678799306,
11599417334472130768, 6238612790024380338, 5169904363978189823, 7160809539609098601, 12992409764685655940,
46126266471470262, 10975409403196300321, 2822430957430394448, 561618328704124539, 13015414181365874062,
4465746772996349942, 1946244108629035845, 17496828967649386175, 2883024483759645523, 9106632430547631473,
13879158535770788154, 1449927166809016346, 14614884650813985783, 2952081054174096768, 17036002253188898606,
4531681534472198908, 11577486448795369550, 10837072743969543298, 12051547528202478533, 5783728287759472282,
18098289833339012493, 13252235870351971830, 14103712702359684079, 1223729028314539628, 14592417029992571650,
18039612578869902675, 17193897619119733648, 11043125045711151195, 14982534958960668951, 3089369946320278360,
5120286904412619024, 9789832226516317024, 6058871823453067274, 15189692145057131153, 8423972467055778937,
14558024088030872276, 9179815255197845682, 6268215501147642559, 4068807096471783550, 9376338163246609614,
15659110180810009697, 13784166615742853035, 2149237202389966706, 12909605674098009100, 4178659878546906659,
386171243290034795, 640035863051577378, 13531439356002704432, 10333841085396651538, 13250664030524184275,
10984823263288084821, 4398816809351718018, 2227301125824267897, 10582190547628267274, 7690326403748531133,
4823114757291188337, 13420063566851055161, 13177408549284069980, 12565544035126800166, 17939014829891499155,
2052414212683224346, 15729787035875247197, 984852478312379439, 1691429919500338054, 15126788187735519683,
13185548047199638235, 8290631933941479723, 14384910153229902993, 16419313701465794768, 15157831870514471650,
7878819826499035512, 16743790405399159823, 17818409006594143176, 10870548103367800908, 6182379021744496101,
13771591897120552595, 14091609828646341284, 8883752405514875020, 2428949568465870470, 9038655354939666725,
14822175004957379533, 11785598603923442959, 15525090770401718210, 17129396081869776328, 4844988293365633673,
3394087862029566757, 283013482195826647, 7923449144205210053, 14931345664767205936, 14307438406331844917,
18112862860532905408, 16408547874516976536, 15928373626736996427, 7835252777128489741, 11821867505128027872,
8047360945397927461, 8770300901650141562, 3779042834167253410, 15775694398361153017, 2141174510233722041,
16746524597408465490, 7341790008784164965, 2341219693951301371, 7907960169056120377, 1866418208095900742,
8705958826821949735, 2264107857566613176, 2051068484314622067, 13520261746726239754, 7896935623171530338,
1477733438141003484, 12900185010441369428, 2215986326227644212, 15188579624681553725, 1445879912504428202,
10001140629665727188, 919275336056348351, 881910076455206420, 10763956104934902577, 12883701938510673118,
3493389299362541501, 1102012684889457917, 6363660147494145051, 12838404370061053708, 302594053210483323,
12695916392667139321, 4491642455408651053, 10835895161012463339, 9500370325485747687, 8816101479115663336,
14015793041123991766, 10473514613415508059, 2420752425683866584, 9333610794264192963, 17486395573854624103,
6222148759667480220, 7354202688450786808, 7055280611641651360, 12324672561820883332, 10828172580964923874,
11567039300035425616, 17727890609821153696, 7478984332956550502 ]

# normal version forward NTT, output in bit reversed order		
		
def NTT_1024(a):
	k = 0
	len = 512
	N = 1024
	
	while len > 0:
		start = 0
		while start < N:
			k = k+1
			zeta = psi_2048_rev[k]
			for j in range(start,start+len):
				t = (zeta * a[j + len])%q
				a[j + len] = (a[j] - t)%q
				a[j] = (a[j] + t)%q
			start = start + 2*len
		len >>= 1

# normal version backward NTT, input is assumed bit reversed, output normal ordering

def iNTT_1024(a):

	k = 1024
	len = 1
	N = 1024

	while (len < N):
		start = 0
		while start < N:
			k = k-1
			zeta = -psi_2048_rev[k]
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
# the original a has to be ordered as [[a[i], a[i+128], a[i+1], a[i+1+128]]

def NTT_4_1024(a):
	k = 0
	len = 512
	N = 1024
	
	while len > 2:
		start = 0
		while start < (N>>2):
			k = k+1
			zeta = psi_2048_rev[k]
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
		zeta = psi_2048_rev[k]
		A1 = a[j]
		t1 = (zeta * A1[1])%q
		t2 = (zeta * A1[3])%q
		a[j] = [(A1[0] + t1)%q, (A1[2] + t2)%q, (A1[0] - t1)%q, (A1[2] - t2)%q]	
		
  # doing len = 1

	for j in range(0,N>>2):
		k = k+1
		zeta = psi_2048_rev[k]
		A1 = a[j]
		t1 = (zeta * A1[1])%q
		k = k+1
		zeta = psi_2048_rev[k]  
		t2 = (zeta * A1[3])%q
		a[j] = [(A1[0] + t1)%q, (A1[0] - t1)%q, (A1[2] + t2)%q, (A1[2] - t2)%q]
		
# 4 element version
# a is now array of quadruples where input is assumed in bit reversed order
# chopped up into 4 element arrays, i.e. [[a[i], a[i+1], a[i+2], a[i+3]]
# for i = 4*k and a the original array in bit reversed order

# the output is equal to map_4_lin_offset(original)

def iNTT_4_1024(a):
	
	k = 1024
	len = 1
	N = 1024
   
	# len = 1 separate

	for j in range(N >> 2):
		k = k-1
		zeta1 = -psi_2048_rev[k]
		k = k-1
		zeta2 = -psi_2048_rev[k]
		
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
			zeta1 = -psi_2048_rev[k]
			k = k-1
			zeta2 = -psi_2048_rev[k]
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
	zeta = -psi_2048_rev[k]
	zeta_invN = (zeta*inv_N)%q
 
	for j in range(N>>2):
		A = a[j]
		a[j] = [((A[0] + A[1])*inv_N)%q, (zeta_invN*(A[0] - A[1]))%q, ((A[2] + A[3])*inv_N)%q, (zeta_invN*(A[2] - A[3]))%q];
  
  
# maps a linear list to 4 element tuple list interleaving offset	
	
def map_lin_4_offset(a):
	res = []
	for i in range(N>>2):
		res.append([a[2*i], a[2*i+hN], a[2*i+1], a[2*i+1+hN]])
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
		

# makes copy to avoid overwriting, this can be avoided if inputs can be overwritten		
		
def fast_mul_1024(a, b):
	if (four_elem):
		a4 = map_lin_4_offset(a)
		b4 = map_lin_4_offset(b)
		NTT_4_1024(a4)
		NTT_4_1024(b4)
		for j in range(256):
			for i in range(4):
				a4[j][i] = (a4[j][i]*b4[j][i])%q;
		iNTT_4_1024(a4)
		return map_4_lin_offset(a4)
	else:
		ac = a.copy()
		bc = b.copy()
		NTT_1024(ac)
		NTT_1024(bc)
		for j in range(1024):
			ac[j] = (ac[j]*bc[j])%q;
		iNTT_1024(ac)
		return ac

"""
N = 1024
a = [random.randint(0,q-1) for _ in range(N)]
b = [random.randint(0,q-1) for _ in range(N)]

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
	
c = schoolbook_mult(a, b)
d = fast_mul_1024(a, b)

print("Sanity check :", d == c)	
"""