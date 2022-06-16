%nprocs=28
%mem=128GB
%chk=W001_1_6_F_2_SPE.chk
# PBE1PBE gen pseudo=read int=(grid=ultrafine) empiricaldispersion=GD3BJ nmr=giao 

W001_1_6_F_2 title 

0 1
C    0.1093222749    -3.917480125    4.1432312141
C    0.0873616624    -3.3438947335    5.4601449646
C    1.4111067422    -3.4448296612    6.0104178482
C    2.2523575338    -4.0783613846    5.0325637995
C    1.4439767853    -4.3739340833    3.8782572335
H    -0.770032047    -2.8677433371    5.9371404972
H    1.7333620193    -3.059183207    6.9781870168
H    3.3202313428    -4.272488852    5.1399286283
H    1.7810215596    -4.8471204923    2.9545216889
H    -0.7333741537    -3.9706366951    3.4526644923
Fe    1.4380680851    -2.3567169838    4.2790266178
C    1.7076407305    -0.4139894044    4.8645698906
C    2.8751600987    -0.9164626437    4.2119801009
C    2.5197129982    -1.3242325191    2.8732650312
C    1.0829827458    -1.1081412516    2.7214460868
C    0.6113033501    -0.5185059652    3.9520432816
H    1.6515678915    -0.0564058624    5.8928432417
H    -0.4207246064    -0.2483632985    4.1690334787
C    0.1813055531    -1.4882625515    1.5627583627
P    0.6667711529    -0.4790050323    0.0514746236
H    0.3439068296    -2.5522450327    1.3169170668
C    -1.3085986107    -1.2768383487    1.887214157
H    -1.9518606498    -1.6846943762    1.0935291115
H    -1.5496402733    -1.7991331369    2.8277574769
H    -1.5396189245    -0.2045578872    1.9967226902
C    -2.7947776073    -0.6098485045    -3.0648763738
C    -2.1762088238    0.6071269769    -2.7575099299
C    -1.1605911642    0.6729400494    -1.7909076329
C    -0.7551360542    -0.4957573091    -1.1290937046
C    -1.3548101257    -1.7255293387    -1.4560415084
C    -2.3790659242    -1.777298085    -2.4093173092
H    -3.5889103606    -0.6514327037    -3.8141825001
C    -2.6420090743    1.8816735868    -3.4241365172
H    -0.7004849923    1.6312441174    -1.5386108808
H    -1.0261998573    -2.6541680688    -0.97887355
C    -2.975772004    -3.1103403785    -2.789357975
C    3.3065520258    -3.1135428297    -2.763767148
C    2.8331020526    -3.6651283387    -1.5667814088
C    2.0599224136    -2.9004256426    -0.6839994706
C    1.755413333    -1.5646595066    -0.9851088812
C    2.2131159077    -1.0164661971    -2.1967245326
C    2.9888616422    -1.7852738228    -3.0746884849
H    3.9128198796    -3.711513002    -3.4491572707
C    3.1445267564    -5.1094022645    -1.2554703853
H    1.7261649115    -3.3603525247    0.246720876
H    1.9535841826    0.0133975104    -2.4640054988
C    3.5471606359    -1.153125637    -4.326854518
F    4.7648230894    -0.5763679868    -4.0750936239
F    3.7389625181    -2.0662709191    -5.3120606502
F    2.7370758015    -0.1762662408    -4.8052955306
F    4.4733925767    -5.3714999368    -1.4089774098
F    2.8046740446    -5.4471891337    0.0202823746
F    2.4762848433    -5.9492347515    -2.0895518095
F    -1.6105574915    2.7427482764    -3.619880881
F    -3.2100053028    1.6296209241    -4.6367385765
F    -3.5710176163    2.5210959165    -2.6675531021
F    -2.2797150881    -3.6916221317    -3.8058561275
F    -2.9553856579    -3.9822660269    -1.7424780557
F    -4.2638033362    -2.9859654809    -3.2023371754
P    3.6696447635    0.9677317447    0.9996964852
C    5.2930076556    1.9861077258    5.2354639239
C    4.0261172856    2.4493982344    4.8556210058
C    3.5208819631    2.1688775203    3.5784154815
C    4.2903041372    1.4172983224    2.6741188117
C    5.5649669739    0.9483286035    3.0567233412
C    6.0644244285    1.2373744687    4.3323524524
H    5.6832609644    2.2087426111    6.2340218639
H    3.4204232719    3.0343606049    5.55501748
H    2.5370498699    2.5363216851    3.2773135799
H    6.1628551641    0.351668188    2.359880331
H    7.0569154965    0.8761476042    4.6210363485
Cl    2.2507760496    3.7788330574    1.0107960021
Cl    -0.8018907663    2.2373967223    0.9829123243
H    3.8821832952    -0.9671170456    4.6256427359
C    5.6101126979    -2.5998456874    0.1346806202
C    5.0852288361    -3.4852499997    1.0805592258
C    4.0707407306    -3.0523768358    1.9422755783
C    3.5544897054    -1.7446020619    1.8882814251
C    4.1356616028    -0.8288913667    0.9594112136
C    5.1466016225    -1.2783128294    0.0884619685
H    6.3890495304    -2.9261199614    -0.5615380215
H    5.4450713195    -4.5168134127    1.1347435295
H    3.6492557845    -3.732563679    2.6841140016
H    5.5942549535    -0.5826022028    -0.6253466581
C    6.3901037795    3.1864293137    -2.0394250716
C    5.5301128044    2.166927945    -2.4734582149
C    4.7264423275    1.4925164161    -1.548259665
C    4.7940983113    1.8093501963    -0.1774726932
C    5.6534225274    2.8369213696    0.2529954444
C    6.4452999358    3.5204235983    -0.67942431
H    7.0114119917    3.7244615646    -2.7624120614
H    5.474827328    1.8921242259    -3.5309729557
H    4.054428391    0.7045848144    -1.8949942722
H    5.6902479868    3.1119426863    1.3098423055
H    7.1040071758    4.3252011711    -0.3376948048
Pd    1.4767138803    1.5573509575    0.6736678196

H P C Cl F 0
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
%chk=W001_1_6_F_2_SPE.chk
# PBE1PBE gen pseudo=read int=(grid=ultrafine) empiricaldispersion=GD3BJ pop=nbo guess=read geom=check

W001_1_6_F_2 title 

0 1

H P C Cl F 0
def2TZVP
****
Pd Fe 0
SDD
****

Pd Fe 0
SDD

