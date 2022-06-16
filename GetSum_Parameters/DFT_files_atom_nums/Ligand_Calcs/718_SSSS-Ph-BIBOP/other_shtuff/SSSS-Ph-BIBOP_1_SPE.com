%NProcShared=32
%mem=18GB
%chk=SSSS-Ph-BIBOP_1_SPE.chk
# PBE1PBE gen pseudo=read int=(grid=ultrafine) empiricaldispersion=GD3BJ nmr=giao 

 title 

0 1
P    11.2195918776    -5.9872549894    7.2637737982
P    13.2395018649    -3.8238971229    6.2665034446
O    10.5044998872    -6.4005059808    4.6932483076
O    11.7957201357    -3.6277098865    3.993456248
C    11.5778008542    -7.5956904738    6.4683139447
C    12.1331117751    -8.7879849309    6.9909623827
C    12.0553884571    -9.9512495844    6.1951831625
H    12.4972722103    -10.8757463624    6.5782214017
C    11.4744182086    -9.9249843329    4.9197643527
H    11.4423337512    -10.8420214856    4.3228909223
C    10.9633790652    -8.7382862369    4.3780496849
H    10.5272888323    -8.6946643844    3.3767396659
C    11.026462572    -7.5838802031    5.1652554008
C    10.9164978602    -5.2986794319    5.5330385975
H    10.0972797453    -4.5618801703    5.5326115246
C    12.1802233248    -4.6494545332    4.9406680888
H    12.7525386406    -5.4110745347    4.3871516126
C    11.7574936582    -2.3957099063    4.606911846
C    11.1767414904    -1.3060433732    3.9498900494
H    10.7265602116    -1.4387618121    2.9627619242
C    11.1974096168    -0.0659089013    4.6014342225
H    10.7451288983    0.8024208485    4.1118525832
C    11.7547055372    0.0747771171    5.8799259877
H    11.7125564592    1.0394335416    6.3940607798
C    12.332322582    -1.0225958933    6.5541824003
C    12.3556466287    -2.2683400442    5.8828655845
C    9.5028209029    -6.1657062243    8.0571986335
C    9.1665821978    -4.7931525237    8.6688888896
H    8.1629364259    -4.842453331    9.1291965582
H    9.1417355564    -3.995539284    7.9033422642
H    9.8976251805    -4.5204527129    9.447512041
C    8.4477234405    -6.5707349571    7.0115438731
H    8.6898500313    -7.5321559333    6.5297186406
H    8.3036239521    -5.8107339682    6.2260486195
H    7.4825451117    -6.6932796409    7.5361549663
C    9.5972680279    -7.25792172    9.1399348492
H    8.621393339    -7.3221571696    9.6555413153
H    10.37436746    -7.0262705357    9.8845501895
H    9.8099103679    -8.2438183735    8.6931378062
C    14.9712625111    -3.6438980501    5.506715592
C    15.6756692857    -2.471443025    6.2161731884
H    15.7163057418    -2.6209511554    7.306124458
H    16.7092688473    -2.4014169144    5.8300671363
H    15.1688997273    -1.5150897723    6.0030870195
C    14.8872933864    -3.3532601127    3.9970307335
H    15.9184696621    -3.2240832275    3.6204768389
H    14.335588173    -2.4229317486    3.7843763715
H    14.4257337407    -4.1740887365    3.4238182387
C    15.6960429171    -4.9774231103    5.7665900871
H    15.1707738394    -5.830288463    5.2979972923
H    15.7941627687    -5.1668609142    6.8480372112
H    16.7076410383    -4.927303445    5.3241253167
Pd    12.9075833483    -4.8374160075    8.2863288432
Cl    15.0200391549    -4.1672115483    9.194538448
Cl    12.1908102781    -5.3673597012    10.5093606161
C    14.1284561721    -8.9353241965    10.7960040975
C    13.1407149705    -9.8825184137    10.4834701809
C    12.488670726    -9.8376802855    9.2460076084
C    12.8129764152    -8.8360893989    8.3073660165
C    13.8192230153    -7.902875298    8.6261983003
C    14.4672813777    -7.9449398909    9.8654491637
H    14.6238195481    -8.9623554408    11.7717779412
H    12.8664419611    -10.6516862035    11.2129623072
H    11.6997931793    -10.5604562001    9.0098576209
H    14.101237254    -7.1365058703    7.8948912029
H    15.2098955849    -7.1773015644    10.1018062759
C    13.8170620543    -0.5255754807    10.5579243528
C    12.9326302175    -1.5709552798    10.2637171461
C    12.4476628603    -1.7268589342    8.9606486077
C    12.8530703874    -0.8538371646    7.9317880537
C    13.7334228407    0.2037897788    8.241318551
C    14.2142970441    0.3624596139    9.5459528016
H    14.2090993862    -0.4093815084    11.5733336494
H    12.6368726842    -2.2938313098    11.0296064416
H    11.741697663    -2.535987519    8.7402216036
H    14.0618023944    0.8804326272    7.4445616488
H    14.9136717049    1.1744733124    9.7707024736

O H P C Cl 0
def2TZVP
****
Pd 0
SDD
****

Pd 0
SDD

--Link1--
%NProcShared=32
%mem=18GB
%chk=SSSS-Ph-BIBOP_1_SPE.chk
# PBE1PBE gen pseudo=read int=(grid=ultrafine) empiricaldispersion=GD3BJ pop=nbo guess=read geom=check

 title 

0 1

O H P C Cl 0
def2TZVP
****
Pd 0
SDD
****

Pd 0
SDD

