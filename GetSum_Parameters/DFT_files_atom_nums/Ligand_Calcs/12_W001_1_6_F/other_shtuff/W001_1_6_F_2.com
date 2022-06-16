%nprocs=14
%mem=128GB
%chk=W001_1_6_F_2.chk
# PBEPBE gen 6D pseudo=read scf=maxcycles=500 denfit int=(grid=ultrafine) empiricaldispersion=GD3BJ opt freq=noraman

 W001_1_6_F

0 1
C    -0.52710    -3.70550    4.97250
C    -0.24320    -2.89710    6.11910
C    1.15070    -3.02500    6.41790
C    1.72720    -3.91550    5.45790
C    0.69070    -4.33450    4.56160
H    -0.94800    -2.26480    6.64080
H    1.68180    -2.50800    7.20490
H    2.76990    -4.19620    5.40650
H    0.80910    -4.99320    3.71300
H    -1.48820    -3.79720    4.48700
Fe    0.88620    -2.27250    4.51590
C    0.83830    -0.22790    4.57050
C    2.16050    -0.69130    4.26870
C    2.12020    -1.45690    3.05480
C    0.75140    -1.45910    2.63530
C    -0.04850    -0.71860    3.56030
H    0.55180    0.34880    5.43920
H    -1.11840    -0.56450    3.52870
C    0.23960    -2.06930    1.34420
P    0.71950    -0.75100    0.13810
H    0.75330    -3.01060    1.15670
C    -1.27520    -2.32790    1.39490
H    -1.61920    -2.88950    0.52890
H    -1.54890    -2.92420    2.26400
H    -1.84540    -1.39850    1.43450
C    -2.69010    -0.28360    -3.09800
C    -2.15220    0.86360    -2.48150
C    -1.14330    0.71780    -1.50480
C    -0.66670    -0.55570    -1.12580
C    -1.18940    -1.68530    -1.79250
C    -2.20590    -1.56300    -2.76160
H    -3.46570    -0.17560    -3.84200
C    -2.66240    2.24610    -2.86110
H    -0.72840    1.60880    -1.05960
H    -0.80570    -2.67040    -1.57050
C    -2.74850    -2.81220    -3.43810
C    3.59260    -2.45530    -3.19210
C    3.42410    -3.18690    -2.00120
C    2.59420    -2.67190    -0.98090
C    1.92580    -1.43700    -1.13500
C    2.12140    -0.71830    -2.33580
C    2.94460    -1.21790    -3.36710
H    4.22080    -2.85020    -3.97740
C    4.13530    -4.52220    -1.83820
H    2.48230    -3.24720    -0.07740
H    1.62100    0.22800    -2.48440
C    3.11740    -0.44440    -4.66560
F    4.35340    -0.60330    -5.13960
F    2.24860    -0.90750    -5.56210
F    2.89650    0.85740    -4.48410
F    5.43830    -4.30830    -1.66380
F    3.66480    -5.20850    -0.79590
F    3.96480    -5.25680    -2.93650
F    -1.70390    3.16100    -2.71430
F    -3.07040    2.26410    -4.13010
F    -3.68880    2.56300    -2.07420
F    -1.83390    -3.28080    -4.28540
F    -3.00040    -3.74460    -2.51900
F    -3.87270    -2.56080    -4.10860
P    3.72760    0.65760    1.31640
C    6.55570    2.41050    4.67090
C    5.17020    2.32430    4.89760
C    4.32410    1.80540    3.89860
C    4.84100    1.37130    2.65910
C    6.23610    1.45850    2.44990
C    7.08890    1.97520    3.44430
H    7.20700    2.81020    5.43460
H    4.75370    2.65920    5.83640
H    3.26430    1.74760    4.09240
H    6.66320    1.12840    1.51440
H    8.15220    2.03830    3.26390
Cl    1.89480    3.27680    2.35520
Cl    -0.87100    1.83490    1.55030
H    3.04530    -0.50760    4.86280
C    5.67260    -3.05560    1.25690
C    4.78820    -3.90240    1.94210
C    3.61150    -3.37580    2.50140
C    3.29610    -2.00570    2.36870
C    4.15060    -1.16140    1.61120
C    5.35480    -1.69570    1.10120
H    6.59220    -3.44880    0.84630
H    5.02440    -4.95120    2.05480
H    2.97100    -4.05040    3.04440
H    6.05620    -1.05770    0.58570
C    5.82730    2.43220    -2.51120
C    5.62700    1.04290    -2.42200
C    5.00110    0.49300    -1.28650
C    4.57950    1.31760    -0.22150
C    4.78190    2.71040    -0.32740
C    5.40260    3.26740    -1.46190
H    6.30300    2.85690    -3.38380
H    5.94970    0.40150    -3.22950
H    4.84220    -0.57240    -1.24360
H    4.46340    3.36080    0.47600
H    5.55310    4.33560    -1.52540
Pd    1.38070    1.23410    1.35200

F Cl H P C 0
6-31G(d)
****
Fe Pd 0
LANL2DZ
****

Fe Pd 0
LANL2DZ
