%nprocs=28
%mem=128GB
%chk=XantPhos_2_SPE.chk
# PBE1PBE gen pseudo=read int=(grid=ultrafine) empiricaldispersion=GD3BJ nmr=giao 

XantPhos_2 title 

0 1
P    -1.62940157    1.6450186821    -3.4632554537
C    -2.1684513049    2.236278409    -1.7957049443
C    -3.5198612774    2.3176229737    -1.3994828427
C    -3.8542170425    2.5139184907    -0.0545444126
C    -2.8564044288    2.5785720177    0.9325229471
C    -1.502389973    2.4740203421    0.5853907393
C    -1.2017472699    2.3299753139    -0.780639186
O    0.1207977689    2.1652870494    -1.1466650232
C    -0.3113563005    2.4828297544    1.5586242233
C    -0.749766344    2.2200970972    3.006778556
C    0.3811030994    3.8741507109    1.4781927762
C    -2.4313378737    -0.0294079784    -3.3913744239
C    -3.2219027597    -0.5344890082    -4.4406194103
C    -3.8204423928    -1.7978895095    -4.3311586896
C    -3.6357054184    -2.5712873988    -3.1788008979
C    -2.8470084467    -2.0752227757    -2.1299487277
C    -2.2499982441    -0.8156632655    -2.2344139971
C    -2.6087651937    2.4930662599    -4.7625256467
C    -2.3072328335    2.1439646176    -6.0954785459
C    -4.0326993354    3.6524412238    -6.8943969098
C    -3.0235739974    2.7131439503    -7.1536370307
C    -4.3105308312    4.0265797936    -5.5734406504
C    -3.6039123424    3.4510925103    -4.5084519963
H    -4.3125552467    2.1754369069    -2.1391947478
H    -3.1459217199    2.6887064817    1.9815641784
H    -4.9070765898    2.5810630197    0.2357821516
H    -1.2436208111    1.2396124527    3.1154588287
H    -1.4469421219    3.0036273246    3.3462312866
H    0.7090550154    4.1016734642    0.4515084655
H    -0.3207021156    4.6620892698    1.8006865143
H    -1.4807714865    1.4547266936    -6.3025212394
H    -4.5847118668    4.1080282166    -7.7225943269
H    -5.0727906854    4.7845977188    -5.3657572315
H    -3.7992847796    3.7876703356    -3.4880321879
H    -3.3914624114    0.066807874    -5.3377024983
H    -4.4424068716    -2.1713351896    -5.1515805646
H    -4.1051580781    -3.5568180717    -3.0947601585
H    -2.6910938994    -2.6733868909    -1.2267636934
H    -1.6471584411    -0.4375197702    -1.4023417673
P    1.6402149661    0.1395295491    -2.8010053897
C    1.6859548713    0.3760554687    -0.962582909
C    2.4872002761    -0.4235901523    -0.1188079957
C    2.3683960583    -0.3276791461    1.2730331233
C    1.4541419935    0.5664844122    1.8555395691
C    0.6792726737    1.4177836912    1.0553255065
C    0.8422921064    1.3124228969    -0.3381063801
C    0.7886395387    -1.4715939856    -3.061023298
C    0.7996089916    -2.519180805    -2.1221281443
C    0.2321044367    -3.7574522223    -2.447302042
C    -0.3452248983    -3.9603005197    -3.7092207492
C    -0.362401635    -2.9183249427    -4.6468195846
C    0.196555627    -1.6762306708    -4.3224814079
C    3.3887473225    -0.3132789082    -3.1558743117
C    4.4311739335    0.472261285    -2.6213043873
C    6.0712213449    -0.9754529108    -3.6793215914
C    5.7641432192    0.1331570172    -2.8760685002
C    5.0381152858    -1.7487282397    -4.2238204461
C    3.6997073518    -1.4231735859    -3.9622246958
H    3.2028723226    -1.1218072573    -0.5636573857
H    1.362391136    0.6061650264    2.9448905241
H    2.989805122    -0.9605327687    1.9139546159
H    0.1184049602    2.2528900308    3.6852420452
H    1.2644700177    3.8939819021    2.1391391448
H    4.1962849767    1.3495166549    -2.0136377901
H    2.900434504    -2.0395810911    -4.3832619763
H    1.24358806    -2.367042539    -1.1340698223
H    0.2400819971    -4.5675918237    -1.7104206365
H    -0.7909749859    -4.9285986623    -3.958690873
H    -0.8284454928    -3.065156584    -5.6260030617
H    0.1753438619    -0.8479026768    -5.0403715485
H    5.2693831623    -2.6146785949    -4.8526110318
H    7.1157864497    -1.2329139344    -3.8828125244
H    6.5663237805    0.7486949494    -2.4566328879
H    -2.7750110988    2.4403494665    -8.1842514098
Pd    0.6715873016    1.9958538996    -3.7990545389
Cl    2.8695394594    2.7085244148    -4.2902739439
Cl    -0.0948225765    4.1957422601    -4.3120367274

C O P Cl H 0
def2TZVP
****
Pd 0
SDD
****

Pd 0
SDD

--Link1--
%nprocs=28
%mem=128GB
%chk=XantPhos_2_SPE.chk
# PBE1PBE gen pseudo=read int=(grid=ultrafine) empiricaldispersion=GD3BJ pop=nbo guess=read geom=check

XantPhos_2 title 

0 1

C O P Cl H 0
def2TZVP
****
Pd 0
SDD
****

Pd 0
SDD

