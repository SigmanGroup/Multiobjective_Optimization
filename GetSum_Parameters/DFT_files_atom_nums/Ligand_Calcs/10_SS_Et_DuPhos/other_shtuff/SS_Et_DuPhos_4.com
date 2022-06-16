%nprocs=14
%mem=128GB
%chk=SS_Et_DuPhos_4.chk
# PBEPBE gen 6D pseudo=read scf=maxcycles=500 denfit int=(grid=ultrafine) empiricaldispersion=GD3BJ opt freq=noraman

 SS_Et_DuPhos

0 1
Cl    -1.22910    -3.83410    1.27990
P    -1.64560    -0.64770    0.11840
P    1.49600    -0.11030    0.27360
C    -1.01470    1.06010    -0.33270
C    0.38260    1.29320    -0.28180
C    -3.02390    -0.63550    1.33280
H    -2.82680    -1.43930    2.04340
C    -2.61050    -1.11240    -1.37040
H    -3.00040    -0.20230    -1.82650
C    2.41510    0.64950    1.66690
H    2.49780    1.72070    1.48210
C    2.92020    -0.41810    -0.84520
H    3.01550    -1.50010    -0.94580
C    -4.25120    -1.04360    0.49650
H    -4.96290    -1.61300    1.09530
H    -4.77490    -0.14710    0.16200
C    -3.80630    -1.83930    -0.74770
H    -4.62350    -1.94400    -1.46200
H    -3.51520    -2.84720    -0.44900
C    3.81510    0.06020    1.47000
H    4.56280    0.61330    2.03920
H    3.84000    -0.96910    1.83020
C    4.13290    0.08350    -0.03920
H    5.02400    -0.50370    -0.26390
H    4.35820    1.11120    -0.32720
C    -3.20890    0.66580    2.13630
H    -3.47080    1.48720    1.47060
H    -4.06300    0.54380    2.80330
C    -1.86420    -1.87510    -2.48520
H    -1.28040    -1.15590    -3.05880
H    -2.59850    -2.26990    -3.18830
C    1.81120    0.49470    3.07850
H    1.00360    1.21820    3.18540
H    2.56010    0.79470    3.81220
C    2.81110    0.17680    -2.26210
H    2.77690    1.26430    -2.21390
H    3.72170    -0.06500    -2.81140
C    -1.98850    1.06650    2.97500
H    -2.20720    1.94660    3.57980
H    -1.13480    1.31260    2.34670
H    -1.69030    0.26470    3.65070
C    -0.93080    -3.01110    -2.05220
H    -1.46470    -3.76880    -1.48020
H    -0.50840    -3.50770    -2.92590
H    -0.09260    -2.64580    -1.46340
C    1.26820    -0.88710    3.45990
H    2.03850    -1.65300    3.38080
H    0.41780    -1.17290    2.84540
H    0.92270    -0.88570    4.49380
C    1.60720    -0.33010    -3.06650
H    1.60200    -1.41820    -3.13080
H    0.66640    -0.01460    -2.61980
H    1.62840    0.06280    -4.08300
C    -1.87990    2.09730    -0.73740
H    -2.94470    1.92600    -0.78560
C    -1.36610    3.36160    -1.08050
H    -2.03550    4.15160    -1.38780
C    0.02000    3.59520    -1.02300
H    0.41640    4.56520    -1.28500
C    0.88990    2.56370    -0.62350
H    1.95210    2.75180    -0.58080
Cl    2.07620    -3.35030    1.19080
Pd    0.18610    -2.06130    0.74020

H Cl C P 0
6-31G(d)
****
Pd 0
LANL2DZ
****

Pd 0
LANL2DZ

