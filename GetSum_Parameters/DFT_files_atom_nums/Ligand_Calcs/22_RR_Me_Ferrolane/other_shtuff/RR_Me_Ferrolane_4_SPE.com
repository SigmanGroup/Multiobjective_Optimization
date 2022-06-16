%nprocs=28
%mem=128GB
%chk=RR_Me_Ferrolane_4_SPE.chk
# PBE1PBE gen pseudo=read int=(grid=ultrafine) empiricaldispersion=GD3BJ nmr=giao 

RR_Me_Ferrolane_4 title 

0 1
Fe    -0.0370662432    2.5918621113    0.0403403441
P    -1.7938316979    -0.2680628413    -0.0241488298
P    1.7832073222    -0.2284114458    0.0875397825
C    -1.5201008736    1.4255227836    -0.6332657467
C    -0.6251460669    1.7631870706    -1.7169616912
C    -0.5911478322    3.1906854953    -1.8486885895
C    -1.4762706175    3.7502756661    -0.8677294592
C    -2.0573175739    2.6706789509    -0.1218040615
C    1.4717171503    1.4548903729    0.7069340263
C    1.9809378963    2.7148399956    0.2031252311
C    1.3759362244    3.7766076118    0.9555674839
C    0.503603099    3.1913973679    1.9330273282
C    0.569498005    1.7658460123    1.7926001783
C    -3.3340205925    -0.8519827398    -0.9631508792
C    -4.3533607439    -1.1157621684    0.170292584
C    -2.6535102785    -0.0168499067    1.6337756014
C    -4.1558664658    -0.0561257159    1.2689840535
C    3.3361702837    -0.7834336025    1.0230157336
C    4.361031775    -1.0174576524    -0.1119826834
C    2.6368730545    0.052040487    -1.5688029544
C    4.1397721555    0.0441553056    -1.2042177519
H    -0.0386709956    1.0441719569    -2.2900260648
H    0.0360537674    3.7528619681    -2.5413713847
H    -1.6375099969    4.8126725018    -0.6807712819
H    -2.747071936    2.7756940186    0.7157550496
H    2.6680949768    2.8403559701    -0.6337473508
H    -0.0007024633    1.0304109951    2.361234621
Cl    -1.642077455    -3.4367519133    0.0208706956
Cl    1.7023665096    -3.3993356385    0.0230489283
Pd    0.0109793128    -1.7053080792    0.0272895583
H    -3.0526125912    -1.8006114207    -1.4453596896
C    -3.8323769442    0.1525014674    -2.0087735029
H    -2.3589484092    0.9808311537    2.0035951462
C    -2.2434914258    -1.0780226812    2.6645588907
H    -4.1662543567    -2.1219591724    0.582313153
H    -5.384036104    -1.0995126567    -0.2282534973
H    -4.4780968141    0.934444811    0.8955710837
H    -4.7640452481    -0.2786889538    2.1644472636
H    1.5133661412    4.8434635309    0.7750965533
H    -0.1359879863    3.7351779837    2.6290586934
H    -3.06608191    0.3711612267    -2.7710065663
H    -4.1364161997    1.1130039396    -1.5580356212
H    -4.7110739187    -0.2726501057    -2.5256055449
H    -2.7846266859    -0.9054604829    3.6122432369
H    -2.4587096473    -2.097476171    2.3065290235
H    -1.1590809314    -1.0408570053    2.8681963569
H    3.0761003051    -1.7410307383    1.4994440263
C    3.8120511598    0.2255736381    2.074736583
H    2.3200341831    1.0451158429    -1.9325546604
C    2.2505767554    -1.011750379    -2.606030475
H    4.1964296182    -2.025064451    -0.5301175904
H    5.391128779    -0.9805863466    0.2866841029
H    4.4398111609    1.0393963191    -0.8247821506
H    4.7526780635    -0.1593055984    -2.1010009969
H    4.7000883613    -0.182955968    2.5890038611
H    4.0944912815    1.195353315    1.6298499152
H    3.0411379654    0.4224119569    2.8382592784
H    2.4885785192    -2.0282938718    -2.2541995949
H    2.7876134959    -0.8213731918    -3.5526343697
H    1.1655836837    -0.997610928    -2.8094653788

P Cl H C 0
def2TZVP
****
Pd Fe 0
SDD
****

Pd Fe 0
SDD

--Link1--
%nprocs=28
%mem=128GB
%chk=RR_Me_Ferrolane_4_SPE.chk
# PBE1PBE gen pseudo=read int=(grid=ultrafine) empiricaldispersion=GD3BJ pop=nbo guess=read geom=check

RR_Me_Ferrolane_4 title 

0 1

P Cl H C 0
def2TZVP
****
Pd Fe 0
SDD
****

Pd Fe 0
SDD

