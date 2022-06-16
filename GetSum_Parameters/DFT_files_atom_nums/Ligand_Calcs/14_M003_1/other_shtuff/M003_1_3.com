%nprocs=14
%mem=128GB
%chk=M003_1_3.chk
# PBEPBE gen 6D pseudo=read scf=maxcycles=500 denfit int=(grid=ultrafine) empiricaldispersion=GD3BJ opt freq=noraman

 M003_1

0 1
Fe    0.04400    -1.66260    1.72780
P    1.80110    -0.10330    -0.73930
P    -1.78020    0.16630    -0.57010
C    1.55140    -1.53780    0.42380
C    0.67120    -2.60870    0.04210
C    0.62210    -3.57000    1.09450
C    1.45010    -3.11570    2.15940
C    2.05830    -1.89190    1.72900
C    -1.46200    -0.42850    1.19630
C    -1.98480    -1.51580    1.99750
C    -1.33140    -1.54470    3.26980
C    -0.44800    -0.42940    3.31620
C    -0.53080    0.24930    2.06300
C    3.42410    -0.60780    -1.54720
C    4.19840    0.34930    -2.23890
C    5.43240    0.00690    -2.82970
C    5.89230    -1.32160    -2.74760
C    5.13070    -2.29920    -2.07820
C    3.90470    -1.93520    -1.48210
C    2.27710    1.50630    0.10440
C    3.58040    1.74610    0.59110
C    3.93900    2.99150    1.15300
C    2.97320    4.01220    1.24280
C    1.67240    3.80410    0.74520
C    1.34170    2.56000    0.16910
C    -2.39470    1.90580    -0.21470
C    -2.90880    2.27630    1.04900
C    -3.39140    3.57850    1.29180
C    -3.35190    4.53490    0.25750
C    -2.83880    4.19120    -1.00890
C    -2.37020    2.87900    -1.23600
C    -3.37890    -0.62060    -1.14960
C    -4.64450    -0.16770    -0.71520
C    -5.83390    -0.78130    -1.16120
C    -5.75920    -1.86620    -2.05650
C    -4.50940    -2.32680    -2.51370
C    -3.32990    -1.70330    -2.05250
H    0.10700    -2.67420    -0.87640
H    0.04730    1.12220    1.81750
H    3.84850    1.36810    -2.32070
H    6.83410    -1.59530    -3.20000
H    3.34500    -2.69590    -0.95820
H    4.32890    0.97150    0.51210
H    3.24350    4.96100    1.68370
H    0.34560    2.43230    -0.22380
H    -2.95770    1.56180    1.85890
H    -3.71640    5.53670    0.43200
H    -2.00000    2.63040    -2.21900
H    -4.72650    0.65050    -0.01680
H    -6.66440    -2.35120    -2.39250
H    -2.37640    -2.06920    -2.40490
Cl    1.47200    0.19440    -4.04680
Cl    -1.71630    0.46580    -3.87690
Pd    -0.04400    0.16400    -2.27740
C    3.07010    -1.09760    2.52950
H    0.19580    -0.17630    4.14440
H    -1.49370    -2.26400    4.05920
C    -3.12300    -2.46050    1.66360
H    2.79560    -0.05670    2.36540
N    4.40580    -1.22950    1.94490
C    5.39450    -0.40230    2.63420
C    4.93070    -2.59200    1.85810
H    6.32540    -0.35960    2.06700
H    5.03730    0.61760    2.76430
H    5.62620    -0.79480    3.62610
H    5.84430    -2.60080    1.26300
H    5.18290    -3.00170    2.83500
H    4.22830    -3.26950    1.37140
C    2.90310    -1.50440    6.87780
C    3.20730    -2.62940    6.08990
C    3.26230    -2.51390    4.68750
C    2.99940    -1.27770    4.05230
C    2.71820    -0.15090    4.85830
C    2.66260    -0.26320    6.26100
H    2.86270    -1.59160    7.95420
H    3.40340    -3.58180    6.56130
H    3.50350    -3.38820    4.10280
H    2.54370    0.81470    4.40620
H    2.44060    0.60470    6.86570
N    -4.17210    -2.37570    2.68460
H    -3.58030    -2.16470    0.72780
C    -5.32970    -3.20130    2.33550
C    -4.67080    -1.01390    2.88600
H    -5.77060    -2.89520    1.38580
H    -5.06750    -4.25710    2.25940
H    -6.10140    -3.12800    3.10330
H    -5.10020    -0.59960    1.97520
H    -5.44400    -0.99900    3.65550
H    -3.88490    -0.34030    3.22920
C    -1.83060    -6.58910    1.15540
C    -1.90210    -6.02960    2.44530
C    -2.30230    -4.68920    2.60970
C    -2.62800    -3.89330    1.48690
C    -2.56170    -4.46760    0.19870
C    -2.16410    -5.80860    0.03230
H    -1.52890    -7.61910    1.02850
H    -2.39320    -4.27830    3.60500
H    -2.83510    -3.89660    -0.67150
C    -3.94470    3.92630    2.66540
F    -4.54830    5.11490    2.66680
F    -4.82950    3.00000    3.03440
F    -2.94870    3.94770    3.54930
C    -2.80660    5.21770    -2.13150
F    -1.76130    5.00010    -2.93030
F    -3.92480    5.11070    -2.84630
F    -2.71930    6.45540    -1.64500
C    -7.18150    -0.29430    -0.65100
F    -8.14610    -0.55130    -1.53380
F    -7.46010    -0.92170    0.49010
F    -7.14930    1.01910    -0.42560
C    -4.43310    -3.51030    -3.46630
F    -4.43990    -4.63410    -2.75200
F    -5.47590    -3.51540    -4.29570
F    -3.31270    -3.46760    -4.18700
C    6.25030    1.06320    -3.55730
F    6.07680    2.25810    -2.99090
F    7.54910    0.76480    -3.52310
F    5.85450    1.12610    -4.82700
C    5.64540    -3.72840    -1.98720
F    6.51770    -3.81690    -0.98440
F    4.65150    -4.59030    -1.77110
F    6.25890    -4.07030    -3.12030
C    5.34810    3.24390    1.66870
F    6.21890    2.38160    1.14450
F    5.73550    4.47990    1.35670
F    5.35180    3.11300    2.99410
C    0.61660    4.89520    0.81710
F    1.11420    6.04480    1.27000
F    0.10660    5.09300    -0.39710
F    -0.35970    4.49490    1.62900
H    -1.66230    -6.63190    3.30990
H    -2.12260    -6.24260    -0.95700
H    0.02510    -4.46750    1.11160
H    1.55390    -3.61200    3.11120

H Cl P N C F 0
6-31G(d)
****
Fe Pd 0
LANL2DZ
****

Fe Pd 0
LANL2DZ
