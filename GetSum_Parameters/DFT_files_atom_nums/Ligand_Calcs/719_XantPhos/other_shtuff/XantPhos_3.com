%nprocs=14
%mem=128GB
%chk=XantPhos_3.chk
# PBEPBE gen 6D pseudo=read scf=maxcycles=500 denfit int=(grid=ultrafine) empiricaldispersion=GD3BJ opt freq=noraman

 XantPhos

0 1
P    -1.65770    1.60310    -3.54460
C    -2.08790    2.30900    -1.85870
C    -3.39370    2.68080    -1.47310
C    -3.64870    3.12320    -0.15930
C    -2.61310    3.15170    0.79730
C    -1.30390    2.75880    0.43350
C    -1.06090    2.40610    -0.91340
O    0.19750    2.07150    -1.30350
C    -0.13150    2.62340    1.42850
C    -0.64200    2.37940    2.86950
C    0.69180    3.92300    1.41750
C    -2.36270    -0.13480    -3.51170
C    -2.29320    -0.93490    -4.67240
C    -2.85390    -2.22640    -4.69320
C    -3.47380    -2.74280    -3.54090
C    -3.52760    -1.96580    -2.36940
C    -2.97440    -0.67080    -2.35610
C    -2.89880    2.58940    -4.57070
C    -2.90920    3.99920    -4.47320
C    -4.75710    4.14440    -6.05440
C    -3.82610    4.77390    -5.20920
C    -4.76760    2.74170    -6.15650
C    -3.84680    1.97230    -5.41870
H    -4.20660    2.61540    -2.18240
H    -2.85190    3.45530    1.80420
H    -4.65140    3.40980    0.12380
H    -1.23380    1.46530    2.93050
H    -1.25690    3.20250    3.23060
H    1.07590    4.14120    0.42030
H    0.08970    4.77680    1.72900
H    -2.19940    4.49520    -3.82670
H    -5.46170    4.73570    -6.62120
H    -5.48470    2.25310    -6.80020
H    -3.89130    0.89950    -5.51090
H    -1.82030    -0.54940    -5.56380
H    -2.79780    -2.82540    -5.59030
H    -3.89640    -3.73690    -3.55240
H    -3.99440    -2.36160    -1.47930
H    -3.03100    -0.09400    -1.44540
P    1.53870    -0.06710    -2.86700
C    1.53130    0.14420    -0.99960
C    2.21180    -0.72700    -0.12280
C    2.18390    -0.50410    1.26860
C    1.45840    0.58010    1.80380
C    0.75740    1.45650    0.94270
C    0.84130    1.23230    -0.45140
C    0.82460    -1.75580    -3.24690
C    0.24740    -2.56420    -2.24400
C    -0.24080    -3.84910    -2.55130
C    -0.16350    -4.33650    -3.86940
C    0.39170    -3.53140    -4.88130
C    0.87510    -2.24570    -4.57030
C    3.40380    -0.34130    -2.96670
C    4.26920    0.67330    -2.49930
C    6.22300    -0.67980    -3.03680
C    5.66740    0.51120    -2.53610
C    5.37600    -1.70520    -3.49450
C    3.97800    -1.53710    -3.45640
H    2.76580    -1.56500    -0.52200
H    1.44540    0.71040    2.87420
H    2.71550    -1.17450    1.92870
H    0.17540    2.29910    3.58470
H    1.54810    3.85580    2.08910
H    3.85430    1.59200    -2.11030
H    3.35790    -2.34890    -3.80100
H    0.17420    -2.20530    -1.22850
H    -0.67940    -4.45810    -1.77450
H    -0.53870    -5.32180    -4.10470
H    0.44750    -3.89990    -5.89530
H    1.31210    -1.63910    -5.35050
H    5.79810    -2.62460    -3.87380
H    7.29540    -0.80720    -3.06710
H    6.31190    1.30190    -2.18070
H    -3.81220    5.85090    -5.12510
Pd    0.67600    1.77270    -4.20540
Cl    2.83290    2.09350    -5.04330
Cl    0.08990    3.56830    -5.57920

H Cl P C O 0
6-31G(d)
****
Pd 0
LANL2DZ
****

Pd 0
LANL2DZ

