%nprocs=28
%mem=128GB
%chk=XantPhos_5_SPE.chk
# PBE1PBE gen pseudo=read int=(grid=ultrafine) empiricaldispersion=GD3BJ nmr=giao 

XantPhos_5 title 

0 1
P    -1.4620568065    1.4944080069    -3.6695148714
C    -2.0736515974    1.6587610823    -1.9221004185
C    -3.4353258101    1.5234355096    -1.5786023119
C    -3.8378306284    1.6208308889    -0.2404118302
C    -2.9060532225    1.8920370571    0.7780313334
C    -1.5530699447    2.1001978123    0.4733657466
C    -1.1884205097    1.9777237204    -0.8800347238
O    0.1265045192    2.1804832911    -1.2225212338
C    -0.4470807834    2.5218006772    1.4652336669
C    -0.8206046779    2.1800340757    2.9155707045
C    -0.2441747491    4.0584157575    1.3322095253
C    -1.8940083212    -0.1938898718    -4.2478384643
C    -2.7877125769    -1.0615183845    -3.5989969194
C    -3.1416462623    -2.2769002474    -4.2011468611
C    -2.6142572183    -2.6290957417    -5.4511669336
C    -1.7077665256    -1.7723132398    -6.0941097191
C    -1.3402467739    -0.5627284191    -5.492661105
C    -2.7298573125    2.5360849828    -4.507921643
C    -3.6690857395    1.9963885922    -5.4035093605
C    -4.7170935997    4.1784852283    -5.6253728177
C    -4.659352424    2.8205687376    -5.9603379591
C    -3.7773310044    4.7189080106    -4.7332679214
C    -2.7835417096    3.9064915853    -4.180462397
H    -4.1759616227    1.3481089255    -2.3653795559
H    -3.2524806888    1.9655754389    1.8132379711
H    -4.8946253991    1.4975947958    0.0152925816
H    -0.9772608056    1.0967300227    3.0525926861
H    -1.7404173187    2.7088886178    3.2148078968
H    0.0361572632    4.3363352118    0.3035607255
H    -1.1751714519    4.5885718066    1.5966293606
H    -3.6340143489    0.9362524764    -5.6683898959
H    -5.4898189747    4.8193879993    -6.0623765636
H    -3.8081332681    5.7830994216    -4.4787842526
H    -2.0351860929    4.336636891    -3.5090248714
H    -3.1934307423    -0.8011445229    -2.6177715398
H    -3.8301839201    -2.9539536884    -3.6849303819
H    -2.9009261305    -3.576139966    -5.9199685358
H    -1.2821866849    -2.0469071417    -7.0647726375
H    -0.6260788997    0.1101280162    -5.982790036
P    1.7191311091    0.3067298024    -2.6328876289
C    2.0772917828    0.8873799965    -0.9152262033
C    3.0803951315    0.4152277181    -0.0522587774
C    2.9764100138    0.6514557718    1.3273784673
C    1.8550623034    1.3011395258    1.8724191422
C    0.8543981364    1.8209975053    1.0343240234
C    1.0522299302    1.6559962223    -0.342435826
C    0.6562307548    -1.1699227674    -2.2592656596
C    -0.0561013446    -1.3288746252    -1.0563919434
C    -0.8323394534    -2.476943089    -0.83816605
C    -0.9000036823    -3.4801292894    -1.8120387691
C    -0.199767923    -3.324403375    -3.0180059755
C    0.5653778387    -2.1767392625    -3.2445458365
C    3.2160132595    -0.5552250244    -3.2668013858
C    3.7039561117    -1.6884026436    -2.5790818207
C    5.4304434873    -1.9756292446    -4.2616456974
C    4.8153233089    -2.3847187358    -3.0690758812
C    4.9288844716    -0.8706055072    -4.9623195374
C    3.8273568455    -0.1572152438    -4.4694030094
H    3.9213736634    -0.1573968621    -0.4530144244
H    1.7663936056    1.3965275889    2.9586647263
H    3.7616910184    0.2851886714    1.9957044538
H    -0.0268459778    2.5087624035    3.6061403245
H    0.5595969265    4.3918290969    2.0106283553
H    3.1988306989    -2.0418797098    -1.674545206
H    3.4518768136    0.7231130966    -4.9971694378
H    0.0087655768    -0.5747882665    -0.2679058921
H    -1.3735300781    -2.5877058506    0.1074102483
H    -1.4986223029    -4.3798775981    -1.6361829137
H    -0.2573515486    -4.0955718867    -3.7922303911
H    1.1107531203    -2.0671959108    -4.1867587969
H    5.3976654514    -0.5531353607    -5.8992020666
H    6.2949294573    -2.5254380531    -4.6481341688
H    5.1922908235    -3.2563210447    -2.5239729613
H    -5.3862172104    2.3928680889    -6.6586797638
Pd    0.7953296767    2.0130345081    -3.9194003744
Cl    3.0092115056    2.8952107078    -3.7117995908
Cl    0.2029494624    3.7615852306    -5.3935901889

C O P Cl H 0
def2TZVP
****
Pd 0
SDD
****

Pd 0
SDD

--Link1--
%nprocs=28
%mem=128GB
%chk=XantPhos_5_SPE.chk
# PBE1PBE gen pseudo=read int=(grid=ultrafine) empiricaldispersion=GD3BJ pop=nbo guess=read geom=check

XantPhos_5 title 

0 1

C O P Cl H 0
def2TZVP
****
Pd 0
SDD
****

Pd 0
SDD

