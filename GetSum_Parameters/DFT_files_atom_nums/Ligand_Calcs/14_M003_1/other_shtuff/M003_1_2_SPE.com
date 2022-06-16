%nprocs=28
%mem=128GB
%chk=M003_1_2_SPE.chk
# PBE1PBE gen pseudo=read int=(grid=ultrafine) empiricaldispersion=GD3BJ nmr=giao 

M003_1_2 title 

0 1
Fe    -0.0062403935    -1.4820415156    1.682749451
P    1.7931195319    -0.0596371444    -0.791460811
P    -1.7898273307    0.29775939    -0.5460468326
C    1.5048323527    -1.407853469    0.3905751188
C    0.5637822107    -2.4419150158    -0.0055125036
C    0.4763529415    -3.4026550887    1.0465034313
C    1.3386888691    -2.9786266443    2.1011165119
C    2.0045156965    -1.7559193306    1.7204303548
C    -1.5210086732    -0.2686540076    1.1776498452
C    -2.0287852142    -1.3661087504    2.0035274901
C    -1.3179322759    -1.3039511662    3.2560634033
C    -0.4155313129    -0.1981487936    3.2400033097
C    -0.5304463622    0.4432476102    1.9683594369
C    3.3451606431    -0.5977483326    -1.6157360258
C    4.2466620311    0.3414144378    -2.146972269
C    5.4318576128    -0.0952894572    -2.7505841069
C    5.7220270248    -1.4628521165    -2.8502858778
C    4.8089645444    -2.3941536104    -2.3428032796
C    3.6223679041    -1.9673785581    -1.7302767532
C    2.2738940083    1.5200308167    0.0276784927
C    3.5087395945    1.7009163662    0.6792287483
C    3.8230706346    2.9549364433    1.2248740224
C    2.9241920769    4.028572921    1.1326007314
C    1.7075585936    3.8474534339    0.4689271818
C    1.3905170925    2.604495787    -0.0962028946
C    -2.4730695389    1.9832919373    -0.2110052664
C    -3.1399857208    2.2143144271    1.0032445242
C    -3.6804926225    3.47639371    1.2860006085
C    -3.5410065364    4.5239084128    0.3709837997
C    -2.8659367576    4.2937129878    -0.8369287314
C    -2.3374506996    3.033151217    -1.1389401204
C    -3.2171108066    -0.6620997871    -1.1968587798
C    -4.5395732559    -0.3225973502    -0.8855291509
C    -5.5973735435    -1.0972534835    -1.3879688372
C    -5.3444178043    -2.2005898015    -2.2089218688
C    -4.0180547994    -2.5198864285    -2.5400514545
C    -2.9581312084    -1.7548471649    -2.0422345741
H    -0.0020653438    -2.4403752516    -0.9385718192
H    0.0429262525    1.3075342992    1.6342389403
H    4.0254791917    1.4096653446    -2.1049527528
H    6.6431974201    -1.7976011081    -3.3331261144
H    2.9203315134    -2.7111668425    -1.345560662
H    4.2087157583    0.8584121908    0.7528063122
H    3.182434479    5.002788011    1.5537313186
H    0.4608934041    2.4793174802    -0.662963765
H    -3.2419568688    1.413116408    1.7384214582
H    -3.9469153465    5.5132057289    0.5979417391
H    -1.8347101556    2.8583565592    -2.0948904468
H    -4.7599540982    0.5495102807    -0.2635992297
H    -6.1716455127    -2.7916264895    -2.6084147247
H    -1.9318483107    -1.9765899541    -2.3474206163
Cl    1.502266104    0.8138007113    -3.8111319889
Cl    -1.8287354227    0.8512782044    -3.6307423939
Pd    -0.0649077395    0.3749846667    -2.0997863714
C    3.0703087905    -1.0276439211    2.5112274501
H    0.2946522342    0.0580683096    4.0262579009
H    -1.4450622548    -2.0118088474    4.0731975049
C    -3.2004005206    -2.3147256151    1.7802221313
H    2.8484374767    0.054843367    2.4408420322
N    4.3953866839    -1.1800988395    1.842506211
C    5.4751606674    -0.5728641909    2.6263226868
C    4.775210836    -2.5340493536    1.440077446
H    6.3602496197    -0.4482911199    1.9793680491
H    5.1699616368    0.417984586    2.9936758551
H    5.766272521    -1.1853130624    3.5075746061
H    5.6168362955    -2.4701899636    0.7308108648
H    5.1011547838    -3.1737506256    2.2907894506
H    3.9417796961    -3.0332562496    0.9247820403
C    3.1078583477    -1.8187699458    6.7969666022
C    3.3729196234    -2.8639140527    5.9004031294
C    3.3515396738    -2.6347723625    4.5174991177
C    3.0602330998    -1.3541513328    4.0038648099
C    2.8290595376    -0.3089739969    4.9185767371
C    2.8410639084    -0.5349987724    6.301419682
H    3.1219435814    -2.0018825107    7.8760037463
H    3.6007954065    -3.8659422386    6.2786365097
H    3.5674656991    -3.4624329164    3.8361116742
H    2.6552649602    0.7044405695    4.5348982365
H    2.6518617313    0.2939568558    6.9912018004
N    -4.1803378113    -2.134246125    2.8829986115
H    -3.6645093618    -2.0647651504    0.804189142
C    -5.3633569    -2.9747696079    2.6525628343
C    -4.6137272747    -0.7403573861    3.0185213487
H    -5.9119895524    -2.6989913756    1.7234833224
H    -5.0709064709    -4.033534756    2.5880188739
H    -6.0530932967    -2.8554948707    3.5048470235
H    -5.1256265064    -0.3477718345    2.1091880672
H    -5.3285581135    -0.6721176759    3.8554908346
H    -3.7577840228    -0.0909169727    3.2589623353
C    -1.8620620164    -6.4420650142    1.4867926898
C    -1.8951077468    -5.8036813681    2.738197886
C    -2.3452502918    -4.4813403679    2.8471328706
C    -2.7569756217    -3.772806095    1.7031962951
C    -2.7523838678    -4.4295330753    0.462918174
C    -2.3045174194    -5.7533489637    0.3483062255
H    -1.5086803622    -7.4749738059    1.405108126
H    -2.4090663376    -3.9987766212    3.8269155122
H    -3.1075591174    -3.9025082083    -0.4251765218
C    -4.3484424937    3.6892959974    2.6223595748
F    -3.4280144405    3.8760300767    3.6087661298
F    -5.1673656679    4.7707121249    2.6168483932
F    -5.0933872794    2.6022259151    2.9783474216
C    -2.7242874796    5.4426449549    -1.8100380799
F    -3.9332751145    5.7821159574    -2.3368136744
F    -2.2432143326    6.5502760608    -1.1766736791
F    -1.8908271342    5.1505672182    -2.8339051451
C    -7.0058424267    -0.7313197902    -0.9880841333
F    -7.2384731566    0.5960162911    -1.1569207401
F    -7.9368056859    -1.4170697978    -1.6954711946
F    -7.2197345892    -1.0005472256    0.3376639478
C    -3.7318438589    -3.7294631642    -3.3980217023
F    -3.6899155287    -4.8688553116    -2.6338028561
F    -4.6969843287    -3.9207470478    -4.3315920043
F    -2.5386888277    -3.6322789887    -4.0310603989
C    6.434300028    0.92286785    -3.2399562484
F    7.3001947804    1.2672557304    -2.242677496
F    7.1797159401    0.4392626833    -4.2682512247
F    5.8304292082    2.0639279835    -3.6576609227
C    5.1374060711    -3.8658058966    -2.3653498871
F    5.9535081089    -4.1895297826    -3.3998759978
F    5.7804478898    -4.2405703146    -1.2158567153
F    4.0177967767    -4.6315585435    -2.455903036
C    5.1221542288    3.1454338978    1.9711346404
F    6.1236199818    2.3991911289    1.4407309037
F    5.5205523985    4.4417788628    1.9746919297
F    4.9855684338    2.7665567661    3.2816881585
C    0.6952668852    4.9634968221    0.3620114301
F    0.2494462667    5.0955244714    -0.9142940399
F    -0.3918905804    4.6953138561    1.1439272543
F    1.2006905859    6.1567828817    0.7585423897
H    -1.5761136938    -6.3433978928    3.6362485768
H    -2.3119993157    -6.239296437    -0.6325616664
H    -0.1855029598    -4.2661876347    1.0681737358
H    1.4240155358    -3.45621732    3.0756712234

C N F P Cl H 0
def2TZVP
****
Pd Fe 0
SDD
****

Pd Fe 0
SDD

--Link1--
%nprocs=28
%mem=128GB
%chk=M003_1_2_SPE.chk
# PBE1PBE gen pseudo=read int=(grid=ultrafine) empiricaldispersion=GD3BJ pop=nbo guess=read geom=check

M003_1_2 title 

0 1

C N F P Cl H 0
def2TZVP
****
Pd Fe 0
SDD
****

Pd Fe 0
SDD

