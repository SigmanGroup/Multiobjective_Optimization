%NProcShared=52
%mem=18GB
%chk=S-35tBu-C3-Tunephos_2_SPE_NoPd.chk
# PBE1PBE/def2TZVP int=(grid=ultrafine) empiricaldispersion=GD3BJ nmr=giao 

 title 

0 1
P    -5.0561439301    -0.126841147    2.5328676737
P    -2.2717728185    0.0789979749    0.7475083024
O    -1.4776188988    3.4741115095    3.7821262411
O    -4.1383408813    4.3827000336    3.0816423434
C    -6.1204959204    -1.1058531543    3.6568506265
C    -5.7246358253    -2.4196201571    3.9729495735
C    -6.4743534027    -3.1899021109    4.8698489225
C    -7.6245770035    -2.6046450217    5.4472433237
C    -8.0418004213    -1.2997021625    5.1462370939
C    -7.2674733679    -0.5538122698    4.2342239265
C    -6.095900528    1.1895344333    1.821202771
C    -6.5610070495    2.2750928523    2.5892786283
C    -7.3849224697    3.2467545806    2.0069080916
C    -7.7007101604    3.1142921363    0.63490943
C    -7.2353350877    2.0507715434    -0.1539985146
C    -6.4312277408    1.0764907971    0.4682053257
C    -3.8989315329    0.7155401313    3.7079034398
C    -3.7373314007    0.1831206083    5.0025530674
C    -2.8233156077    0.7536866091    5.8950849694
C    -2.0632780219    1.8658798562    5.5186561294
C    -2.229443252    2.4070192057    4.2370383671
C    -3.1318627708    1.8361308641    3.3082193051
C    -3.2403802279    2.4762268073    1.9650746195
C    -3.7108107338    3.8079556546    1.9033201607
C    -3.8126994476    4.4928207232    0.6850731798
C    -3.4281061985    3.8503265336    -0.4953114616
C    -2.9510796006    2.5350129962    -0.4603183375
C    -2.8508714475    1.8395845324    0.7607982675
C    -1.4761191127    -0.0565188583    -0.8963174539
C    -2.2789449354    -0.4456915087    -1.9857643475
C    -1.7470377078    -0.4912373723    -3.2795285936
C    -0.3911744464    -0.1291893765    -3.4522359987
C    0.431863272    0.2535027842    -2.3827229683
C    -0.1348647096    0.2835909056    -1.0918165663
C    -0.9549102505    -0.006410663    2.0043227983
C    -0.0152821207    1.0266468587    2.1657190516
C    1.0010382604    0.9111085525    3.1271110817
C    1.015636156    -0.2450864149    3.9326238629
C    0.0631297846    -1.2752706312    3.8160598147
C    -0.917382927    -1.1473504848    2.8194872592
C    -3.4932451332    5.5989836741    3.5155935291
H    -4.8488482902    -2.8376324869    3.4719354517
H    -8.2184841121    -3.2045782607    6.142154082
H    -7.5766833984    0.4577677173    3.9591804294
H    -6.2453966312    2.3656206855    3.6311890177
H    -8.3298574837    3.8795137257    0.1715293241
H    -6.0532736617    0.2130564052    -0.0901759582
H    -4.333976795    -0.6779282877    5.3130362117
H    -2.7012886525    0.3257540112    6.8954477387
H    -1.3329007989    2.298910169    6.2080488343
H    -4.2123724237    5.5115124806    0.6715329869
H    -3.5079486913    4.3734926581    -1.4537114226
H    -2.651180238    2.0444137073    -1.3891003433
H    -3.3106431814    -0.7483206549    -1.793475148
H    0.0316266284    -0.1671482663    -4.4598611058
H    0.4828355506    0.5580570269    -0.2329450345
H    -0.1052992759    1.9306766169    1.5589630494
H    1.7971563435    -0.3399941904    4.6948223043
H    -1.6748928743    -1.9193925418    2.6582183287
H    -2.8173371561    5.9725100181    2.7241198099
C    -2.7387135012    5.3351593715    4.8224700695
C    -1.392161442    4.6322128078    4.6315945794
H    -0.9443572058    4.3773762659    5.6116376381
H    -3.39422764    4.7348867287    5.478830245
H    -2.5530173424    6.2942295178    5.3420609568
C    -7.2725084812    -5.5622841542    4.7469361185
C    -6.1069983976    -4.647353259    5.1941516724
C    -5.8823107414    -4.7953024367    6.7170543847
C    -4.8295187518    -5.1013334631    4.461354533
H    -7.4444271997    -5.4668468445    3.6615031403
H    -7.0344219075    -6.6176913781    4.9702035955
H    -8.2134623483    -5.3104997086    5.2654978908
H    -5.0445807994    -4.1597078264    7.0539791108
H    -6.7779962867    -4.5106294312    7.2954751464
H    -5.6396267247    -5.843758261    6.9661265394
H    -3.9506467242    -4.5036679537    4.7613067079
H    -4.615880946    -6.1546046333    4.7124166819
H    -4.9378135413    -5.0304379995    3.3653300864
C    -10.0514444018    -1.6127401134    6.6978531477
C    -9.2961885862    -0.6527683962    5.7598840798
C    -10.2564914624    -0.219835038    4.6259691091
C    -8.8715560043    0.5936770604    6.5735143652
H    -10.3946512211    -2.5162588138    6.164706769
H    -10.942087218    -1.1067282381    7.1083026098
H    -9.4249571037    -1.9298088728    7.5498486403
H    -9.7858938007    0.5097398611    3.94483755
H    -10.5689644441    -1.0887779528    4.0222678152
H    -11.1609088946    0.2504826763    5.0521021322
H    -8.3498969389    1.3326150537    5.9406815675
H    -9.7579923977    1.0878811051    7.0101808472
H    -8.1908140783    0.3145581474    7.3963554258
C    -9.4798533474    4.4193603893    2.7530060021
C    -7.9341237693    4.4410211803    2.8062811644
C    -7.4067919394    5.7539109909    2.1801495514
C    -7.4978512693    4.403381526    4.2837273393
H    -9.8734635681    3.4905467171    3.2019964771
H    -9.8944609028    5.276264574    3.3136855001
H    -9.8562056794    4.4790402895    1.7179268538
H    -6.3049100281    5.7788064171    2.2274517612
H    -7.7061120478    5.8528722424    1.1230245792
H    -7.8029569189    6.6278965212    2.727842286
H    -6.3991489725    4.4270962615    4.3785678733
H    -7.9084233777    5.2803657313    4.8133836192
H    -7.8731553202    3.5000576824    4.7971500272
C    -8.4458551058    3.0420762802    -2.1759498689
C    -7.5482719635    1.9060883744    -1.6521336307
C    -8.2652132431    0.5555867841    -1.8922819363
C    -6.2163138763    1.9291713856    -2.441197864
H    -7.967144702    4.0304832274    -2.0586383962
H    -8.6435987936    2.8921311227    -3.2510578057
H    -9.4203353049    3.0641208582    -1.6572350978
H    -7.6601461034    -0.3005438228    -1.5495371154
H    -9.2282520663    0.5226265079    -1.3541375913
H    -8.4670481901    0.4200627153    -2.9698446119
H    -5.5410103259    1.1180724893    -2.1204045133
H    -5.6853589501    2.8848952791    -2.2896317885
H    -6.4099566108    1.8032062585    -3.5216712173
C    -4.0160922091    -1.3476425864    -4.0828775897
C    -2.5808893996    -0.9554627187    -4.4851807317
C    -1.9008610631    -2.1947360577    -5.1151781635
C    -2.652901314    0.1862289407    -5.5260686379
H    -4.5745830118    -0.4899786408    -3.6680527427
H    -4.5645504547    -1.699512341    -4.9736179741
H    -4.0253275131    -2.1603892091    -3.3362576723
H    -0.8767864958    -1.9712649193    -5.4599728252
H    -1.8429434405    -3.0183029231    -4.3836212798
H    -2.4821876252    -2.5444610256    -5.9869939725
H    -1.6499912815    0.492715364    -5.8699821608
H    -3.2284687557    -0.1397144374    -6.4108268081
H    -3.1519947626    1.0742864616    -5.099723468
C    2.1017344175    2.1095845663    -2.1052115244
C    1.9114700546    0.6407004219    -2.5546952789
C    2.3851272834    0.5128814231    -4.0145310044
C    2.7883032609    -0.2811194264    -1.6734208779
H    1.8049806441    2.2528675202    -1.0516142311
H    3.1614092819    2.4069761445    -2.2018030468
H    1.4929948067    2.7929114758    -2.7223340039
H    2.2871310767    -0.5223582618    -4.3847205772
H    1.8173182908    1.1775487753    -4.6888464025
H    3.449511443    0.7954042585    -4.085621906
H    2.6654387258    -1.3377387801    -1.9659780841
H    3.854156935    -0.0117511411    -1.7836141322
H    2.5278964438    -0.1969542749    -0.6043909988
C    3.4624794082    1.3989331936    3.0349224151
C    2.0679280246    2.0012466381    3.3270253944
C    2.013195051    2.5050887281    4.7886404269
C    1.8502076332    3.2075920206    2.3934525384
H    3.5230492714    1.041482494    1.9919723539
H    4.2492321003    2.1598867715    3.1852981734
H    3.6851284986    0.54498136    3.6968666508
H    1.018616956    2.9284151149    5.0088130636
H    2.2047022553    1.6938088504    5.510840832
H    2.7724632537    3.2903591345    4.9537063137
H    0.8602931763    3.6669851082    2.5552743527
H    2.6223077178    3.9715925188    2.5890124343
H    1.9304360646    2.9211667501    1.3296249154
C    1.4073533598    -3.2809036729    4.5093197539
C    0.1053599197    -2.4858146752    4.7631300432
C    -1.0999317533    -3.423237931    4.5510870174
C    0.0765947981    -1.9885078252    6.2279366311
H    2.3022140147    -2.655189677    4.6706901291
H    1.4684777102    -4.1451068375    5.1947532181
H    1.439031861    -3.6588908077    3.4732915803
H    -2.0543355053    -2.8938608316    4.722082858
H    -1.1230224734    -3.8494065356    3.5335483119
H    -1.0468572395    -4.2629290397    5.2654195234
H    -0.8380709062    -1.4016245501    6.4202146007
H    0.0913243716    -2.8481723974    6.9212031496
H    0.9450286438    -1.3509587961    6.4659769349
H    -0.6807290036    5.2894490569    4.1034134076
H    -4.2882996136    6.3485458313    3.6787305709

--Link1--
%NProcShared=52
%mem=18GB
%chk=S-35tBu-C3-Tunephos_2_SPE_NoPd.chk
# PBE1PBE/def2TZVP int=(grid=ultrafine) empiricaldispersion=GD3BJ pop=nbo guess=read geom=check

 title 

0 1

