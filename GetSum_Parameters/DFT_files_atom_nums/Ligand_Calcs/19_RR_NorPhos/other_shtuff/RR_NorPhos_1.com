%nprocs=14
%mem=128GB
%chk=RR_NorPhos_1.chk
# PBEPBE gen 6D pseudo=read scf=maxcycles=500 denfit int=(grid=ultrafine) empiricaldispersion=GD3BJ opt freq=noraman

 RR_NorPhos

0 1
P    -1.61350    0.01400    -0.04840
C    -2.36130    1.70650    0.16920
C    -1.80580    2.81820    -0.49750
C    -2.36130    4.10150    -0.32990
C    -3.48100    4.28080    0.50430
C    -4.04560    3.17460    1.16710
C    -3.48730    1.89310    0.99730
H    -3.92630    1.04590    1.50510
H    -4.90930    3.30730    1.80240
H    -3.91030    5.26430    0.63060
H    -1.93110    4.94810    -0.84550
H    -0.95030    2.69170    -1.14470
C    -3.01560    -1.22080    -0.06240
C    -2.98720    -2.39870    0.71640
C    -4.06230    -3.30810    0.67880
C    -5.17630    -3.04960    -0.14140
C    -5.21170    -1.88100    -0.92440
C    -4.13530    -0.97430    -0.88360
H    -4.16830    -0.07980    -1.49000
H    -6.06250    -1.67920    -1.55910
H    -6.00080    -3.74720    -0.17270
H    -4.03170    -4.20630    1.27840
H    -2.14350    -2.62150    1.35070
C    -0.60960    -0.30390    1.44340
C    -1.04980    -0.23240    2.92520
C    -1.20420    1.19900    3.40900
C    0.01040    1.76080    3.42960
H    0.26110    2.77680    3.70020
C    0.99390    0.74330    2.88270
H    2.04130    0.93180    3.12080
C    0.35840    -0.51650    3.50150
H    0.81670    -1.45180    3.17280
H    0.37940    -0.50650    4.59380
C    0.61370    0.63100    1.38110
P    1.65290    -0.04650    0.03700
C    2.36590    -1.64230    0.67960
C    1.89710    -2.87750    0.18660
C    2.42330    -4.08770    0.67860
C    3.42630    -4.06940    1.66630
C    3.90430    -2.83980    2.15750
C    3.37600    -1.63170    1.66330
H    3.75040    -0.69020    2.03930
H    4.67900    -2.82200    2.91040
H    3.83320    -4.99720    2.04220
H    2.06060    -5.03050    0.29490
H    1.13190    -2.90340    -0.57550
C    3.08820    1.12620    -0.19330
C    3.08990    2.42830    0.35390
C    4.19140    3.28470    0.16030
C    5.30270    2.84710    -0.58410
C    5.30960    1.55190    -1.13390
C    4.20700    0.69860    -0.93760
H    4.21930    -0.29510    -1.36340
H    6.15920    1.21100    -1.70790
H    6.14770    3.50340    -0.73470
H    4.18390    4.27910    0.58260
H    2.25040    2.78710    0.92880
Pd    0.06040    -0.21590    -1.76730
Cl    1.76470    -0.45420    -3.33410
Cl    -1.57080    -0.34620    -3.42230
H    0.31080    1.62020    1.03430
H    -2.14620    1.65770    3.67440
H    -1.82670    -0.92860    3.24250
H    -0.28630    -1.33440    1.29160

H Cl C P 0
6-31G(d)
****
Pd 0
LANL2DZ
****

Pd 0
LANL2DZ

