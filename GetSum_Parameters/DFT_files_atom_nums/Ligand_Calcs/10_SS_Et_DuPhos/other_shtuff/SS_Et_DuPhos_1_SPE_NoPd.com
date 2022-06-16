%nprocs=28
%mem=128GB
%chk=SS_Et_DuPhos_1_SPE_NoPd.chk
# PBE1PBE/def2TZVP int=(grid=ultrafine) empiricaldispersion=GD3BJ nmr=giao 

SS_Et_DuPhos_1 title 

0 1
P    -1.5534086934    -0.8283122526    0.0536547218
P    1.5233168    -0.2145646161    0.3788789695
C    -0.9517698817    0.8734874658    -0.3697701895
C    0.4116949182    1.1695910796    -0.1258036721
C    -2.9734455073    -0.7337670365    1.2745165755
H    -3.0312162011    -1.8101134343    1.5417432179
C    -2.657082636    -1.3363184667    -1.3990195133
H    -2.5065609753    -0.5447575035    -2.1583528405
C    2.5588570146    0.3936988029    1.8340681567
H    2.3995441706    1.4917383027    1.8541969955
C    2.9643952792    -0.2803681893    -0.8405315115
H    3.1645712892    -1.3678515315    -0.8925661343
C    -4.1963308421    -0.4009192124    0.4029890385
H    -5.1268077007    -0.5943866086    0.9676480074
H    -4.2094581569    0.671617581    0.1351664094
C    -4.110054283    -1.2914782837    -0.8469238703
H    -4.8096910422    -0.9690417477    -1.6377322579
H    -4.3922379412    -2.321640432    -0.5630491623
C    4.0066181982    0.09216669    1.3959962926
H    4.7195501122    0.6756620624    2.0054520149
H    4.20723554    -0.9804529585    1.5690791374
C    4.1309499333    0.4075151706    -0.0984731427
H    5.0961470793    0.066283125    -0.513438623
H    4.0784957013    1.5025063393    -0.2639104342
C    -2.819237625    0.0499113841    2.5864889482
H    -3.7547015604    -0.1183805678    3.1546724141
H    -2.0223166576    -0.4227764401    3.1850524112
C    -2.3293367778    -2.7061137119    -2.0193129967
H    -3.1373423516    -2.9282804284    -2.7440649912
H    -2.3852177541    -3.4707945672    -1.2242068452
C    2.1905883648    -0.1999071525    3.1988598626
H    2.9265169246    0.1812853143    3.9331306508
H    2.3199029852    -1.2958953463    3.147251494
C    2.7189366456    0.2732870802    -2.2494964456
H    2.4910319041    1.35433419    -2.1896244396
H    3.6791979412    0.1960071804    -2.7954821344
C    -2.5606493647    1.5576406858    2.4829694326
H    -3.3045558311    2.0640372599    1.8434177749
H    -1.5628765586    1.7759592433    2.0698356836
H    -2.6151410609    2.018991458    3.4836263114
C    -0.9700235853    -2.795829481    -2.7096096759
H    -0.8574103132    -2.0334717503    -3.5015478496
H    -0.8332225933    -3.7891540192    -3.1687194027
H    -0.154968691    -2.6590783595    -1.9749207844
C    0.7701957447    0.1339877684    3.6530408929
H    0.042742405    -0.3466209263    2.9748067978
H    0.5781923523    1.2226846565    3.6531712709
H    0.5753943601    -0.2490495579    4.6687217557
C    1.6216638452    -0.4452324974    -3.0353883675
H    1.5652084639    -0.0678609959    -4.0706535583
H    1.8036187936    -1.532344797    -3.0759539921
H    0.633697259    -0.2917509763    -2.5699700751
C    -1.7895348785    1.8768042971    -0.8948495247
H    -2.8366632946    1.6553409707    -1.1198871728
C    -1.2977301165    3.1675746786    -1.1255515058
H    -1.9624260994    3.9388099194    -1.5273029466
C    0.0397538531    3.4725932417    -0.8355878445
H    0.4234832834    4.4834265848    -1.0057531474
C    0.8922869916    2.4745373287    -0.3467178641
H    1.9427845582    2.711142581    -0.1459104603

--Link1--
%nprocs=28
%mem=128GB
%chk=SS_Et_DuPhos_1_SPE_NoPd.chk
# PBE1PBE/def2TZVP int=(grid=ultrafine) empiricaldispersion=GD3BJ pop=nbo guess=read geom=check

SS_Et_DuPhos_1 title 

0 1

