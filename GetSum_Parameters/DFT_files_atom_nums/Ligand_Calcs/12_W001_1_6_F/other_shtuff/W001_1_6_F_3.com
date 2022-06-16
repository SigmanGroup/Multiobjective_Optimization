%nprocs=14
%mem=128GB
%chk=W001_1_6_F_3.chk
# PBEPBE gen 6D pseudo=read scf=maxcycles=500 denfit int=(grid=ultrafine) empiricaldispersion=GD3BJ opt freq=noraman

 W001_1_6_F

0 1
C    -0.67300    -3.88720    4.84200
C    -0.49520    -3.05800    5.99500
C    0.87550    -3.14540    6.39760
C    1.54350    -4.03240    5.49580
C    0.58720    -4.48910    4.53110
H    -1.25250    -2.43790    6.45420
H    1.33300    -2.60500    7.21450
H    2.59380    -4.28680    5.52620
H    0.78410    -5.15500    3.70290
H    -1.59300    -4.00960    4.28850
Fe    0.73370    -2.42330    4.47090
C    0.62920    -0.38100    4.49410
C    1.98220    -0.80820    4.29460
C    2.05370    -1.59570    3.09470
C    0.72060    -1.62960    2.57420
C    -0.16740    -0.90880    3.43040
H    0.26500    0.19970    5.33050
H    -1.23550    -0.78350    3.31740
C    0.32170    -2.22810    1.23980
P    0.76110    -0.80110    0.14050
H    0.92370    -3.11250    1.03700
C    -1.16580    -2.61300    1.20330
H    -1.41570    -3.16450    0.29940
H    -1.42690    -3.26830    2.03280
H    -1.81670    -1.73880    1.25040
C    -2.67540    -0.05110    -3.01410
C    -2.12760    1.03670    -2.30580
C    -1.11140    0.80350    -1.35420
C    -0.63580    -0.49900    -1.09060
C    -1.16580    -1.56430    -1.85080
C    -2.19160    -1.35590    -2.79540
H    -3.45670    0.12200    -3.73960
C    -2.63150    2.45000    -2.55970
H    -0.68820    1.65190    -0.83970
H    -0.77700    -2.56410    -1.72670
C    -2.74380    -2.53850    -3.57600
C    3.77470    -2.30620    -3.17010
C    3.71490    -2.98670    -1.93930
C    2.83200    -2.53200    -0.93570
C    2.00670    -1.40400    -1.14070
C    2.10350    -0.72690    -2.37800
C    2.97250    -1.17200    -3.39640
H    4.44410    -2.65190    -3.94410
C    4.60760    -4.19550    -1.70170
H    2.81100    -3.06640    -0.00080
H    1.49660    0.14670    -2.56500
C    3.03160    -0.45700    -4.73820
F    2.27340    -1.11970    -5.60960
F    2.58140    0.79470    -4.64600
F    4.28370    -0.42710    -5.19550
F    5.82830    -3.77820    -1.37040
F    4.13070    -4.96630    -0.72330
F    4.69080    -4.92960    -2.81090
F    -1.65320    3.33770    -2.37980
F    -3.08610    2.57220    -3.80690
F    -3.62250    2.71880    -1.71180
F    -3.89940    -2.23850    -4.16850
F    -1.86030    -2.89790    -4.50540
F    -2.94340    -3.56740    -2.75180
P    3.76200    0.49150    1.38480
C    6.69860    2.85980    4.22450
C    6.36970    3.38300    2.96100
C    5.49510    2.67820    2.11130
C    4.93590    1.44260    2.50810
C    5.27530    0.93150    3.78060
C    6.15000    1.63130    4.63420
H    7.36830    3.40040    4.87760
H    6.78600    4.32760    2.64220
H    5.25670    3.10370    1.14940
H    4.86700    -0.01070    4.11070
H    6.39930    1.22530    5.60370
Cl    1.98970    2.88170    2.88740
Cl    -0.81080    1.72000    1.76670
H    2.80890    -0.58360    4.95330
C    5.83850    -3.12300    1.79300
C    4.88530    -3.99000    2.34760
C    3.62730    -3.49040    2.72510
C    3.29500    -2.13270    2.52290
C    4.21020    -1.28610    1.84460
C    5.50160    -1.77930    1.55810
H    6.82480    -3.48730    1.54160
H    5.13390    -5.02820    2.51690
H    2.94090    -4.17450    3.19420
H    6.25680    -1.11650    1.16410
C    5.17350    1.83730    -2.91670
C    5.68910    0.64090    -2.38820
C    5.28410    0.19740    -1.11420
C    4.37620    0.95340    -0.33670
C    3.85730    2.14480    -0.89120
C    4.24980    2.58730    -2.16850
H    5.47740    2.17280    -3.89830
H    6.38790    0.05510    -2.96840
H    5.67700    -0.74200    -0.77100
H    3.15830    2.74170    -0.32240
H    3.84680    3.50500    -2.57280
Pd    1.42210    1.05190    1.55530

F Cl H P C 0
6-31G(d)
****
Fe Pd 0
LANL2DZ
****

Fe Pd 0
LANL2DZ
