%NProcShared=32
%mem=18GB
%chk=R-C4-Tunephos_1_SPE_NoPd.chk
# PBE1PBE/def2TZVP int=(grid=ultrafine) empiricaldispersion=GD3BJ nmr=giao 

 title 

0 1
P    -0.7155233465    1.6868614779    -0.0169775876
C    -0.0515713246    2.1557307033    -1.6543209942
C    -0.8618341142    1.9260538647    -2.7831894463
C    -0.3972420749    2.2886485355    -4.0536463859
C    0.8682866576    2.8761235261    -4.2020226534
C    1.6719655806    3.1088359536    -3.0763864003
C    1.2148866102    2.7530024858    -1.8012045522
H    1.8564723141    2.9010457122    -0.9275110057
H    1.2300536905    3.1515750623    -5.1982597527
H    -1.8531037223    1.4787174451    -2.6469189093
C    -1.3495241961    3.2386981778    0.7340331516
C    -0.8440376004    4.4930113713    0.3499283968
C    -1.2847079359    5.6535111937    1.001618869
C    -2.2210803227    5.5661323893    2.0404106424
C    -2.7277633133    4.3148083079    2.4232899493
C    -2.3026399027    3.1535379791    1.7682846757
H    -2.7240020872    2.1796197845    2.0360347988
H    -2.565599035    6.4745376026    2.5453264851
H    -0.1185331503    4.5671596083    -0.4656717082
C    0.7473303573    1.298876957    1.0513920379
C    0.8576551555    1.9255387764    2.3096196591
C    1.9120779994    1.6006601334    3.1678798177
C    2.8708400284    0.6514758584    2.7984111218
C    2.763381309    0.0177141487    1.5525613879
C    1.7026517396    0.3308815806    0.6644038041
C    1.7022299115    -0.3499386261    -0.6624388956
C    2.7721913037    -0.052314348    -1.5447105115
C    2.87873504    -0.6900041937    -2.7886410269
C    1.9093608919    -1.6267077567    -3.1625374074
C    0.8457381132    -1.9362209647    -2.3099353755
C    0.7371564896    -1.3065572748    -1.0531644389
P    -0.7293927697    -1.6800671883    0.0150792462
C    -0.0633086215    -2.154805415    1.6498418597
C    1.1881844405    -2.7850922961    1.7884540575
C    1.6498533226    -3.1414765234    3.0617583872
C    0.8662755704    -2.8751129553    4.1941347115
C    -0.3839487865    -2.2538484963    4.0543437523
C    -0.853365985    -1.8917283908    2.7855543213
H    -1.8327883997    -1.4173747988    2.6562024616
H    1.2319325282    -3.1506971316    5.1889104748
H    1.8156916232    -2.9580387046    0.9090555948
C    -1.381172949    -3.223146923    -0.7382922702
C    -2.327058347    -3.1219220237    -1.7777752121
C    -2.7675956479    -4.2755755722    -2.4358397173
C    -2.2838659079    -5.5352668317    -2.0505179478
C    -1.3553148812    -5.6384614777    -1.0062245965
C    -0.8989711792    -4.4856268648    -0.3515733759
H    -0.1801656929    -4.5721512957    0.4686945672
H    -2.6405466912    -6.4376899503    -2.5577035542
H    -2.730784369    -2.1408715328    -2.0469259055
H    0.1024276911    -2.6744772277    -2.6172598301
H    1.9843296171    -2.121235748    -4.1362971965
H    3.6961101171    -0.4491976445    -3.4728926031
O    3.6243391799    0.9379357703    -1.1147104036
C    4.999686816    0.9053173391    -1.522468204
H    5.3699990213    1.918535467    -1.2907739226
H    5.093343582    0.7610561786    -2.6146659581
C    5.8017420497    -0.1474909724    -0.7541780138
H    6.835392928    -0.1229499821    -1.1452336012
H    5.3973797056    -1.146848982    -1.0014458077
C    5.7994625527    0.0679926045    0.7805448273
H    6.830507629    0.0273816483    1.1770865184
H    5.4091296245    1.0733438204    1.0260661932
C    4.9766198785    -0.9724939804    1.5434771259
H    5.3327910134    -1.9912705228    1.3139203843
O    3.6040466861    -0.9832624654    1.1257156419
H    5.0649826764    -0.8300775148    2.6363124575
H    3.6802917558    0.3981517142    3.4875362705
H    1.9877397815    2.0925930886    4.1428998966
H    0.1221681501    2.6726567254    2.6143197622
H    -0.8972762833    6.629057106    0.6897888171
H    -3.47366032    4.2423114127    3.2211183122
H    -1.026453758    2.1088166172    -4.931269032
H    2.6635107067    3.559948917    -3.1896618291
H    -3.5076977985    -4.1904839719    -3.2378052158
H    -0.9863086386    -6.6204410078    -0.6922588025
H    2.6296650971    -3.6190547359    3.1683734696
H    -0.9972976089    -2.0472727568    4.9372840498

--Link1--
%NProcShared=32
%mem=18GB
%chk=R-C4-Tunephos_1_SPE_NoPd.chk
# PBE1PBE/def2TZVP int=(grid=ultrafine) empiricaldispersion=GD3BJ pop=nbo guess=read geom=check

 title 

0 1

