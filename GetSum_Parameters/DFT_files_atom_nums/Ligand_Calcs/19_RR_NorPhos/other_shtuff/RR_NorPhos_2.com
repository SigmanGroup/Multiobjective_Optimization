%nprocs=14
%mem=128GB
%chk=RR_NorPhos_2.chk
# PBEPBE gen 6D pseudo=read scf=maxcycles=500 denfit int=(grid=ultrafine) empiricaldispersion=GD3BJ opt freq=noraman

 RR_NorPhos

0 1
P    -1.66390    0.12300    -0.05010
C    -2.44680    1.80400    0.20000
C    -1.77160    2.94850    -0.27590
C    -2.31800    4.23410    -0.09870
C    -3.55520    4.38830    0.55300
C    -4.24490    3.25400    1.01890
C    -3.69410    1.97020    0.83950
H    -4.24900    1.11630    1.19110
H    -5.20170    3.36680    1.50780
H    -3.97930    5.37320    0.68580
H    -1.79080    5.10170    -0.46850
H    -0.82530    2.84680    -0.78630
C    -2.92470    -1.26340    -0.00380
C    -3.73550    -1.48890    1.13000
C    -4.66370    -2.54760    1.15640
C    -4.78890    -3.40100    0.04520
C    -3.98180    -3.19260    -1.08790
C    -3.05500    -2.13230    -1.10790
H    -2.44320    -1.99490    -1.98540
H    -4.07310    -3.84520    -1.94420
H    -5.50140    -4.21310    0.06150
H    -5.28050    -2.70460    2.02950
H    -3.65440    -0.84670    1.99250
C    -0.63940    -0.23130    1.42690
C    -1.04510    -0.15230    2.91960
C    -1.17190    1.28240    3.40050
C    0.05300    1.82220    3.40670
H    0.32450    2.83450    3.67080
C    1.01220    0.78710    2.84790
H    2.06660    0.95670    3.06870
C    0.36640    -0.46070    3.47300
H    0.80130    -1.40320    3.13310
H    0.40490    -0.45400    4.56480
C    0.60610    0.67450    1.35460
P    1.61450    -0.00640    -0.00430
C    2.24170    -1.66250    0.57240
C    1.69860    -2.85150    0.04330
C    2.16120    -4.10570    0.48650
C    3.17480    -4.17830    1.46080
C    3.72610    -2.99540    1.98840
C    3.26080    -1.74300    1.54360
H    3.69000    -0.83760    1.94860
H    4.50810    -3.04750    2.73220
H    3.53210    -5.14010    1.79970
H    1.74010    -5.01260    0.07700
H    0.92240    -2.80790    -0.70660
C    3.10260    1.10750    -0.18810
C    3.14240    2.40460    0.36930
C    4.28220    3.21710    0.21110
C    5.39370    2.74040    -0.50860
C    5.36240    1.45020    -1.06950
C    4.22150    0.64090    -0.90850
H    4.20400    -0.34870    -1.34340
H    6.21160    1.07990    -1.62530
H    6.26800    3.36300    -0.63240
H    4.30390    4.20800    0.64110
H    2.30310    2.79280    0.92490
Pd    0.00250    -0.00090    -1.79540
Cl    1.69640    -0.15460    -3.38480
Cl    -1.61860    0.09700    -3.46360
H    0.32200    1.67030    1.01370
H    -2.10370    1.75830    3.67220
H    -1.81350    -0.84270    3.26570
H    -0.35860    -1.27410    1.27030

H Cl C P 0
6-31G(d)
****
Pd 0
LANL2DZ
****

Pd 0
LANL2DZ
