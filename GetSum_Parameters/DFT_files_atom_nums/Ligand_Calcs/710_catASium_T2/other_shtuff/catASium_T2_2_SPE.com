%NProcShared=40
%mem=18GB
%chk=catASium_T2_2_SPE.chk
# PBE1PBE gen pseudo=read int=(grid=ultrafine) empiricaldispersion=GD3BJ nmr=giao 

 title 

0 1
C    17.5007432261    4.0962999374    6.5203725225
C    17.6872586834    2.8555878802    7.1612572867
C    18.966087239    2.4469552751    7.5664382988
C    20.0547645726    3.3030215744    7.3109397858
C    19.8930935269    4.5355601258    6.6593381238
C    18.600435089    4.9257487654    6.2604616661
C    15.9022567286    3.7527378139    4.2033278363
C    16.1074206279    4.5830457094    3.091813223
C    16.1866057381    4.0428728857    1.795024935
C    16.053565369    2.6531371178    1.6465814655
C    15.8679077898    1.7978855441    2.751310311
C    15.806671674    2.3610615164    4.0347297122
C    14.6548432327    3.665130843    6.9057944031
C    14.7209617965    3.7265944986    8.4296863011
C    13.8649031707    2.4735223404    8.8100271721
C    12.6626012519    2.9234536299    7.884547506
C    13.4110235871    3.2108557768    6.5585703311
C    13.781255357    4.9159307819    8.8161466136
C    12.3786552999    4.3565070116    8.4764150042
C    14.5103806571    1.1345236871    8.417013384
C    13.5173627975    2.3818885141    10.3054975075
C    11.4112702652    2.0572069834    7.8461688794
C    12.8930017532    3.1224171922    5.1735976459
C    12.5566744873    1.899144878    4.605753436
C    12.5579658178    3.7812751681    2.9097038415
C    12.8714574332    4.2174199418    4.193185107
C    12.5916882559    0.5177505681    5.1917870408
C    12.4845165055    4.4535832589    1.5663441976
C    11.5828954008    6.2391096365    5.7173018823
C    10.4621751336    5.3973282334    5.5507021692
C    9.3119936289    5.5945646648    6.3265781006
C    9.2700396169    6.6290285453    7.2724045748
C    10.3784624328    7.4730929932    7.4312578939
C    11.5328340579    7.2847298672    6.6599575709
C    12.7682509392    7.0020387261    3.1952210349
C    13.8575656773    7.5268250703    2.4749186379
C    13.6384133866    8.2564468212    1.3009032813
C    12.3290913575    8.4940596799    0.8549399659
C    11.2407133133    8.0143639457    1.596879554
C    11.456814497    7.2706775737    2.7655833264
P    15.8468493282    4.5364952979    5.8651201768
P    13.0991376011    5.9882994954    4.6926480147
S    12.2353462612    2.0652050168    2.9064432969
Pd    15.0476893252    6.6742461784    5.76239112
Cl    16.7989264072    7.478883253    7.1372267062
Cl    14.2506789164    8.9171326439    5.4965073201
H    16.8319197958    2.1984315856    7.3462128779
C    19.1708336054    1.1190840302    8.2591282089
H    21.058849429    2.9991069582    7.6334106873
C    21.0689822025    5.4482396407    6.4019834544
H    18.446369319    5.8974442329    5.783439361
H    16.1904289939    5.6651792044    3.2458880725
C    16.4013653221    4.9464485753    0.6025794579
H    16.0952733341    2.2185126952    0.6398217111
C    15.6877617715    0.3116146273    2.551712993
H    15.6590636242    1.7139905217    4.9065156906
H    10.6044492281    6.9024149947    3.3448931844
H    12.1573699446    9.0751372243    -0.0568756669
H    10.2169600708    8.2246944814    1.2702760428
H    14.8708823098    7.3915084306    2.8630283772
H    14.4913067203    8.6652577722    0.7494004213
H    12.3982512527    7.9448663761    6.7778248629
H    8.3746980097    6.7749492185    7.8856654561
H    10.3526870294    8.283300885    8.1670531322
H    10.4965472504    4.5768107761    4.8278067993
H    8.4501002765    4.9320414855    6.1942869556
H    13.0913512578    -0.1890763158    4.5064258367
H    13.1504873681    0.5154396989    6.1377619867
H    11.5795484087    0.1227372825    5.3930637099
H    11.5691185464    5.0554742043    1.442935423
H    13.3395815279    5.1281927947    1.413905087
H    12.5053826383    3.69660381    0.7650941038
H    10.7223199634    2.4020344871    7.0553553664
H    11.6311740602    0.9929655638    7.6651678311
H    10.8743092014    2.1266004852    8.8077368801
H    13.8058905103    0.2977609144    8.5754050268
H    14.8373016451    1.1141701944    7.3658363615
H    15.3941432988    0.9364247186    9.049959477
H    12.8143867249    1.5488072343    10.4851815327
H    14.4309333354    2.171419942    10.8892435617
H    13.0669437189    3.2964232578    10.7166373634
H    15.7319788769    3.7866599862    8.8586111049
H    14.0418846601    5.8201968037    8.2392312722
H    13.8869912182    5.1665385401    9.8850249632
H    11.7465062819    4.2533020528    9.3753889422
H    11.8214736327    4.9749964921    7.7653107341
H    19.7025428217    1.2432565912    9.2192195454
H    18.2087400432    0.6219634373    8.4680005524
H    19.7769484963    0.4307852928    7.6421008606
H    21.2119203414    5.6230927535    5.3204628921
H    20.9084029545    6.4365262084    6.8678912881
H    22.004927476    5.026832914    6.8038886737
H    15.9074257664    -0.249002371    3.4756473721
H    14.6455613669    0.0803349361    2.2607372841
H    16.3427529733    -0.0743991366    1.7521812102
H    16.2354319186    4.4095624344    -0.3457987389
H    15.7180011891    5.8143065685    0.6287411199
H    17.4307890735    5.3484876461    0.5844811321

S Cl C P H 0
def2TZVP
****
Pd 0
SDD
****

Pd 0
SDD

--Link1--
%NProcShared=40
%mem=18GB
%chk=catASium_T2_2_SPE.chk
# PBE1PBE gen pseudo=read int=(grid=ultrafine) empiricaldispersion=GD3BJ pop=nbo guess=read geom=check

 title 

0 1

S Cl C P H 0
def2TZVP
****
Pd 0
SDD
****

Pd 0
SDD

