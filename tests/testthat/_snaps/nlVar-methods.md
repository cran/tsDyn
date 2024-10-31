# simple methods

    Code
      sapply(mods, class)
    Output
      [[1]]
      [1] "VAR"   "nlVar"
      
      [[2]]
      [1] "VAR"   "nlVar"
      
      [[3]]
      [1] "VAR"   "nlVar"
      
      [[4]]
      [1] "VAR"   "nlVar"
      
      [[5]]
      [1] "VAR"   "nlVar"
      
      [[6]]
      [1] "VAR"   "nlVar"
      
      [[7]]
      [1] "VAR"   "nlVar"
      
      [[8]]
      [1] "VAR"   "nlVar"
      
      [[9]]
      [1] "VECM"  "VAR"   "nlVar"
      
      [[10]]
      [1] "VECM"  "VAR"   "nlVar"
      
      [[11]]
      [1] "VECM"  "VAR"   "nlVar"
      
      [[12]]
      [1] "VECM"  "VAR"   "nlVar"
      
      [[13]]
      [1] "VECM"  "VAR"   "nlVar"
      
      [[14]]
      [1] "VECM"  "VAR"   "nlVar"
      
      [[15]]
      [1] "VECM"  "VAR"   "nlVar"
      
      [[16]]
      [1] "VECM"  "VAR"   "nlVar"
      
      [[17]]
      [1] "TVAR"  "nlVar"
      
      [[18]]
      [1] "TVAR"  "nlVar"
      
      [[19]]
      [1] "TVAR"  "nlVar"
      
      [[20]]
      [1] "TVAR"  "nlVar"
      
      [[21]]
      [1] "TVAR"  "nlVar"
      
      [[22]]
      [1] "TVAR"  "nlVar"
      
      [[23]]
      [1] "TVAR"  "nlVar"
      
      [[24]]
      [1] "TVAR"  "nlVar"
      
      [[25]]
      [1] "TVAR"  "nlVar"
      
      [[26]]
      [1] "TVAR"  "nlVar"
      
      [[27]]
      [1] "TVAR"  "nlVar"
      
      [[28]]
      [1] "TVAR"  "nlVar"
      
      [[29]]
      [1] "TVAR"  "nlVar"
      
      [[30]]
      [1] "TVAR"  "nlVar"
      
      [[31]]
      [1] "TVAR"  "nlVar"
      
      [[32]]
      [1] "TVAR"  "nlVar"
      
      [[33]]
      [1] "TVECM" "nlVar"
      
      [[34]]
      [1] "TVECM" "nlVar"
      
      [[35]]
      [1] "TVECM" "nlVar"
      
      [[36]]
      [1] "TVECM" "nlVar"
      
      [[37]]
      [1] "TVECM" "nlVar"
      
      [[38]]
      [1] "TVECM" "nlVar"
      
      [[39]]
      [1] "TVECM" "nlVar"
      
      [[40]]
      [1] "TVECM" "nlVar"
      
      [[41]]
      [1] "TVECM" "nlVar"
      
      [[42]]
      [1] "TVECM" "nlVar"
      
      [[43]]
      [1] "TVECM" "nlVar"
      
      [[44]]
      [1] "TVECM" "nlVar"
      
      [[45]]
      [1] "TVECM" "nlVar"
      
      [[46]]
      [1] "TVECM" "nlVar"
      
      [[47]]
      [1] "TVECM" "nlVar"
      
      [[48]]
      [1] "TVECM" "nlVar"
      

---

    Code
      sapply(mods, print)
    Output
                       Intercept         Trend  dolcan -1   cpiUSA -1
      Equation dolcan 0.01044232 -2.172806e-05  0.9871079 0.000145502
      Equation cpiUSA 0.37899558 -7.098272e-04 -0.1395487 1.002279482
                       Intercept  dolcan -1    cpiUSA -1
      Equation dolcan 0.01278017  0.9870607 6.443813e-05
      Equation cpiUSA 0.45537005 -0.1410917 9.996312e-01
                      dolcan -1    cpiUSA -1
      Equation dolcan 1.0004113 1.200969e-05
      Equation cpiUSA 0.3346056 9.977632e-01
                              Trend  dolcan -1    cpiUSA -1
      Equation dolcan -0.0000766771 0.99195961 0.0003319243
      Equation cpiUSA -0.0027041580 0.03653954 1.0090455285
                       Intercept         Trend dolcan -1    cpiUSA -1  dolcan -2
      Equation dolcan 0.01360877 -2.594234e-05 1.1503909 -0.003742394 -0.1660825
      Equation cpiUSA 0.28320638  3.472133e-04 0.5331984  1.457512609 -0.6049136
                         cpiUSA -2
      Equation dolcan  0.003911853
      Equation cpiUSA -0.459101557
                       Intercept dolcan -1    cpiUSA -1  dolcan -2    cpiUSA -2
      Equation dolcan 0.01636855 1.1506724 -0.003709029 -0.1664027  0.003781728
      Equation cpiUSA 0.24626928 0.5294298  1.457066046 -0.6006284 -0.457359952
                      dolcan -1   cpiUSA -1  dolcan -2    cpiUSA -2
      Equation dolcan 1.1609666 -0.00194731 -0.1602002  0.001956871
      Equation cpiUSA 0.6843077  1.48357158 -0.5073094 -0.484815406
                              Trend dolcan -1    cpiUSA -1  dolcan -2    cpiUSA -2
      Equation dolcan -9.611456e-05 1.1534930 -0.003171405 -0.1628884  0.003578896
      Equation cpiUSA -1.113111e-03 0.5977555  1.469395219 -0.5384425 -0.466030590
                              ECT  Intercept         Trend dolcan -1    cpiUSA -1
      Equation dolcan -0.01506974 0.01098652 -4.515916e-05 0.1655667 -0.003976831
      Equation cpiUSA -0.05033679 0.19306592 -3.133701e-04 0.5871833  0.456867918
                               ECT   Intercept dolcan -1   cpiUSA -1
      Equation dolcan -0.001106223 0.002123837 0.1595354 -0.00300748
      Equation cpiUSA  0.046559351 0.131565673 0.5453308  0.46359448
                                ECT dolcan -1  cpiUSA -1
      Equation dolcan -0.0006768995 0.1681196 0.00242903
      Equation cpiUSA  0.0731547016 1.0770919 0.80037086
                              ECT        Trend dolcan -1    cpiUSA -1
      Equation dolcan 0.001372394 7.603074e-06 0.1610732 -0.001227277
      Equation cpiUSA 0.238600556 6.138198e-04 0.5082187  0.505185783
                              ECT  Intercept         Trend dolcan -1    cpiUSA -1
      Equation dolcan -0.01524527 0.01156071 -4.607701e-05 0.1678529 -0.002542136
      Equation cpiUSA -0.03983784 0.18485483 -2.787031e-04 0.6918440  0.453852695
                        dolcan -2   cpiUSA -2
      Equation dolcan -0.02196011 -0.00308873
      Equation cpiUSA -0.63011337  0.01100919
                               ECT   Intercept dolcan -1    cpiUSA -1   dolcan -2
      Equation dolcan -0.001031649 0.002457958 0.1635127 -0.001892005 -0.03033589
      Equation cpiUSA  0.046135183 0.129795583 0.6655917  0.457785100 -0.68077540
                         cpiUSA -2
      Equation dolcan -0.002298449
      Equation cpiUSA  0.015789314
                                ECT dolcan -1  cpiUSA -1   dolcan -2   cpiUSA -2
      Equation dolcan -0.0007525866  0.172287 0.00152943 -0.02518247 0.001222138
      Equation cpiUSA  0.0608714297  1.128933 0.63845831 -0.40864241 0.201698403
                              ECT        Trend dolcan -1     cpiUSA -1   dolcan -2
      Equation dolcan 0.001657334 8.424941e-06 0.1655519 -0.0006685009 -0.03020354
      Equation cpiUSA 0.230433428 5.927787e-04 0.6550524  0.4838119692 -0.76192515
                         cpiUSA -2
      Equation dolcan -0.001017307
      Equation cpiUSA  0.044131076
      Model TVAR with  1  thresholds
      
      $Bdown
                        Intercept         Trend dolcan -1    cpiUSA -1
      Equation dolcan  0.08210166  0.0007054584 0.9822979 -0.002408619
      Equation cpiUSA -0.55588003 -0.0017977750 0.7181117  1.003209833
      
      $Bup
                        Intercept         Trend  dolcan -1    cpiUSA -1
      Equation dolcan -0.01354042 -0.0001653546  0.9914561 0.0007168335
      Equation cpiUSA  2.20206986  0.0073398658 -0.7468369 0.9702240258
      
      
      Threshold value[1] "1.1439"
      Model TVAR with  2  thresholds
      
      $Bdown
                        Intercept         Trend dolcan -1    cpiUSA -1
      Equation dolcan  0.08210166  0.0007054584 0.9822979 -0.002408619
      Equation cpiUSA -0.55588003 -0.0017977750 0.7181117  1.003209833
      
      $Bmiddle
                      Intercept         Trend  dolcan -1    cpiUSA -1
      Equation dolcan -0.021278 -0.0002016937  0.9952333 0.0008405793
      Equation cpiUSA  2.041256  0.0064145304 -0.6685873 0.9732236367
      
      $Bup
                       Intercept        Trend dolcan -1    cpiUSA -1
      Equation dolcan  0.3579121  0.001426282 0.8674328 -0.005605241
      Equation cpiUSA -1.2735570 -0.003504812 0.1599709  1.021491583
      
      
      Threshold value[1] "1.1439 1.3893"
      Model TVAR with  1  thresholds
      
      $Bdown
                         Intercept dolcan -1    cpiUSA -1
      Equation dolcan  0.007691528 0.9939684 1.357332e-05
      Equation cpiUSA -1.012509544 1.3415788 9.969201e-01
      
      $Bup
                       Intercept  dolcan -1   cpiUSA -1
      Equation dolcan 0.03045117  0.9689214 0.000145664
      Equation cpiUSA 0.46501765 -0.2412925 1.000711673
      
      
      Threshold value[1] "1.2252"
      Model TVAR with  2  thresholds
      
      $Bdown
                        Intercept dolcan -1    cpiUSA -1
      Equation dolcan -0.02978663 1.0310970 8.444228e-06
      Equation cpiUSA -0.27074633 0.5937531 9.970502e-01
      
      $Bmiddle
                       Intercept dolcan -1    cpiUSA -1
      Equation dolcan 0.06718762  0.943391 1.089975e-05
      Equation cpiUSA 1.91845343 -1.055018 9.957302e-01
      
      $Bup
                       Intercept  dolcan -1   cpiUSA -1
      Equation dolcan 0.03045117  0.9689214 0.000145664
      Equation cpiUSA 0.46501765 -0.2412925 1.000711673
      
      
      Threshold value[1] "1.1439 1.2252"
      Model TVAR with  1  thresholds
      
      $Bdown
                      dolcan -1     cpiUSA -1
      Equation dolcan 1.0018378 -5.982604e-06
      Equation cpiUSA 0.3056463  9.994945e-01
      
      $Bup
                      dolcan -1    cpiUSA -1
      Equation dolcan 0.9949031 9.159023e-05
      Equation cpiUSA 0.1554722 9.998859e-01
      
      
      Threshold value[1] "1.2252"
      Model TVAR with  2  thresholds
      
      $Bdown
                      dolcan -1    cpiUSA -1
      Equation dolcan 0.9994742 8.089218e-05
      Equation cpiUSA 0.3063167 9.977088e-01
      
      $Bmiddle
                      dolcan -1    cpiUSA -1
      Equation dolcan 0.9998115 1.872682e-05
      Equation cpiUSA 0.5559962 9.959537e-01
      
      $Bup
                      dolcan -1    cpiUSA -1
      Equation dolcan 0.9949031 9.159023e-05
      Equation cpiUSA 0.1554722 9.998859e-01
      
      
      Threshold value[1] "1.1439 1.2252"
      Model TVAR with  1  thresholds
      
      $Bdown
                             Trend dolcan -1    cpiUSA -1
      Equation dolcan 2.085515e-05 1.0038460 -0.000084888
      Equation cpiUSA 3.944688e-03 0.6854827  0.984569738
      
      $Bup
                              Trend    dolcan -1   cpiUSA -1
      Equation dolcan -0.0002075003  0.970279328 0.000997695
      Equation cpiUSA -0.0013705132 -0.007164849 1.005870624
      
      
      Threshold value[1] "1.2252"
      Model TVAR with  2  thresholds
      
      $Bdown
                             Trend dolcan -1    cpiUSA -1
      Equation dolcan 0.0002877063 1.0280167 -0.001015842
      Equation cpiUSA 0.0010306707 0.4085665  0.993779851
      
      $Bmiddle
                              Trend dolcan -1    cpiUSA -1
      Equation dolcan -0.0001727056 0.9846775 0.0006501501
      Equation cpiUSA  0.0016630729 0.7017294 0.9898733656
      
      $Bup
                              Trend    dolcan -1   cpiUSA -1
      Equation dolcan -0.0002075003  0.970279328 0.000997695
      Equation cpiUSA -0.0013705132 -0.007164849 1.005870624
      
      
      Threshold value[1] "1.1439 1.2252"
      Model TVAR with  1  thresholds
      
      $Bdown
                        Intercept         Trend   dolcan -1   cpiUSA -1  dolcan -2
      Equation dolcan  0.09231288  0.0006969576  1.11640771 0.005716144 -0.1476973
      Equation cpiUSA -0.38461556 -0.0006441612 -0.05789139 1.209808395  0.6664581
                         cpiUSA -2
      Equation dolcan -0.008053514
      Equation cpiUSA -0.210683537
      
      $Bup
                         Intercept         Trend dolcan -1    cpiUSA -1  dolcan -2
      Equation dolcan 0.0006464748 -0.0001118549  1.141001 -0.003541846 -0.1549285
      Equation cpiUSA 1.6398792639  0.0059730476  0.572203  1.334289684 -1.1084546
                        cpiUSA -2
      Equation dolcan  0.00404909
      Equation cpiUSA -0.35828351
      
      
      Threshold value[1] "1.1439"
      Model TVAR with  2  thresholds
      
      $Bdown
                        Intercept         Trend   dolcan -1   cpiUSA -1  dolcan -2
      Equation dolcan  0.09231288  0.0006969576  1.11640771 0.005716144 -0.1476973
      Equation cpiUSA -0.38461556 -0.0006441612 -0.05789139 1.209808395  0.6664581
                         cpiUSA -2
      Equation dolcan -0.008053514
      Equation cpiUSA -0.210683537
      
      $Bmiddle
                        Intercept        Trend dolcan -1    cpiUSA -1  dolcan -2
      Equation dolcan -0.01537662 -0.000157434 1.1806500 0.0009972436 -0.1865112
      Equation cpiUSA  1.50018865  0.005426318 0.5086428 1.3201624875 -0.9522052
                          cpiUSA -2
      Equation dolcan -0.0003233887
      Equation cpiUSA -0.3424836706
      
      $Bup
                      Intercept       Trend dolcan -1   cpiUSA -1   dolcan -2
      Equation dolcan 0.3044483 0.001169578 0.9527351 -0.02974125 -0.07704594
      Equation cpiUSA 1.2815143 0.009182135 2.0104789  1.30743327 -2.39121154
                        cpiUSA -2
      Equation dolcan  0.02531077
      Equation cpiUSA -0.33933921
      
      
      Threshold value[1] "1.1439 1.3893"
      Model TVAR with  1  thresholds
      
      $Bdown
                        Intercept dolcan -1   cpiUSA -1   dolcan -2   cpiUSA -2
      Equation dolcan  0.01395885 1.1438210 0.002560297 -0.15743623 -0.00252461
      Equation cpiUSA -0.64727526 0.8354151 1.394781489  0.01498696 -0.39677145
      
      $Bup
                       Intercept dolcan -1   cpiUSA -1  dolcan -2   cpiUSA -2
      Equation dolcan 0.04017077 1.1335259 -0.01023868 -0.1713070  0.01039812
      Equation cpiUSA 0.28300510 0.4492155  1.31096745 -0.5821929 -0.31059283
      
      
      Threshold value[1] "1.2252"
      Model TVAR with  2  thresholds
      
      $Bdown
                        Intercept dolcan -1   cpiUSA -1   dolcan -2   cpiUSA -2
      Equation dolcan  0.01395885 1.1438210 0.002560297 -0.15743623 -0.00252461
      Equation cpiUSA -0.64727526 0.8354151 1.394781489  0.01498696 -0.39677145
      
      $Bmiddle
                        Intercept dolcan -1    cpiUSA -1  dolcan -2    cpiUSA -2
      Equation dolcan  0.02333395 1.2135735 -0.004115204 -0.2373632  0.004240793
      Equation cpiUSA -0.05032674 0.5533174  1.269024506 -0.3582838 -0.269630556
      
      $Bup
                        Intercept dolcan -1   cpiUSA -1   dolcan -2   cpiUSA -2
      Equation dolcan  0.05798003 0.9631114 -0.02791145 -0.02695616  0.02829938
      Equation cpiUSA -0.65346050 2.0919411  1.32179865 -1.99796630 -0.31587620
      
      
      Threshold value[1] "1.2252 1.3893"
      Model TVAR with  1  thresholds
      
      $Bdown
                      dolcan -1     cpiUSA -1  dolcan -2     cpiUSA -2
      Equation dolcan 1.1452641 -0.0005511544 -0.1435657  0.0005498422
      Equation cpiUSA 0.6624965  1.4943161966 -0.5053663 -0.4946045183
      
      $Bup
                      dolcan -1    cpiUSA -1  dolcan -2    cpiUSA -2
      Equation dolcan 1.1609381 -0.007105683 -0.1655764  0.007202736
      Equation cpiUSA 0.5957397  1.314457492 -0.4870101 -0.314595092
      
      
      Threshold value[1] "1.2263"
      Model TVAR with  2  thresholds
      
      $Bdown
                        dolcan -1   cpiUSA -1  dolcan -2    cpiUSA -2
      Equation dolcan  1.19639319 0.003697079 -0.1986477 -0.003599167
      Equation cpiUSA -0.03686374 1.224304032  0.2957565 -0.226656499
      
      $Bmiddle
                      dolcan -1   cpiUSA -1  dolcan -2    cpiUSA -2
      Equation dolcan 1.1265025 0.003469589 -0.1289219 -0.003432864
      Equation cpiUSA 0.9516258 1.378481433 -0.6090764 -0.380906894
      
      $Bup
                      dolcan -1    cpiUSA -1  dolcan -2    cpiUSA -2
      Equation dolcan 1.1589176 -0.008044606 -0.1627678  0.008132882
      Equation cpiUSA 0.6281012  1.326424793 -0.5220339 -0.326551571
      
      
      Threshold value[1] "1.1439 1.2252"
      Model TVAR with  1  thresholds
      
      $Bdown
                              Trend dolcan -1     cpiUSA -1  dolcan -2     cpiUSA -2
      Equation dolcan -1.622245e-05 1.1450492 -0.0004092175 -0.1449554  0.0004693694
      Equation cpiUSA  3.096566e-03 0.7035062  1.4672230481 -0.2400976 -0.4792437363
      
      $Bup
                              Trend dolcan -1    cpiUSA -1  dolcan -2   cpiUSA -2
      Equation dolcan -0.0002216198 1.1270540 -0.009180805 -0.1574311  0.01024224
      Equation cpiUSA -0.0005662404 0.5091657  1.309155537 -0.4661987 -0.30682914
      
      
      Threshold value[1] "1.2263"
      Model TVAR with  2  thresholds
      
      $Bdown
                              Trend dolcan -1   cpiUSA -1  dolcan -2    cpiUSA -2
      Equation dolcan -6.292985e-05  1.155769 0.001521385 -0.1604489 -0.001290011
      Equation cpiUSA  3.007284e-03  1.053381 1.507577719 -0.6132165 -0.519312189
      
      $Bmiddle
                              Trend dolcan -1    cpiUSA -1   dolcan -2   cpiUSA -2
      Equation dolcan -0.0001870489 1.1844185 -0.005661468 -0.21109514  0.00657114
      Equation cpiUSA  0.0011972003 0.3427475  1.157796683 -0.04378136 -0.16328716
      
      $Bup
                              Trend dolcan -1   cpiUSA -1  dolcan -2   cpiUSA -2
      Equation dolcan -0.0001633523 0.9858677 -0.02757334 -0.0233897  0.02860914
      Equation cpiUSA  0.0035714290 2.1499439  1.31655866 -2.1653563 -0.32545538
      
      
      Threshold value[1] "1.2338 1.3893"
      Model TVECM with  1  thresholds
      
      $Bdown
                              ECT       Const         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan -0.05165772 0.005288736 -3.885903e-05   0.1675227  -0.0121845
      Equation cpiUSA  0.73785216 0.623866396 -1.062156e-03  -0.1019883   0.1395116
      
      $Bup
                               ECT       Const         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan -0.001023032 0.003613178 -3.313512e-05   0.1312604 0.002326573
      Equation cpiUSA  0.015485020 0.101293958 -3.625345e-05   0.2371179 0.620902189
      
      
      Threshold value[1] -0.1176
      Model TVECM with  2  thresholds
      
      $Bdown
                              ECT       Const         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan -0.05165772 0.005288736 -3.885903e-05   0.1675227  -0.0121845
      Equation cpiUSA  0.73785216 0.623866396 -1.062156e-03  -0.1019883   0.1395116
      
      $Bmiddle
                               ECT      Const         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan -0.002269518 0.01440822 -0.0001021736  0.08504187 0.003521756
      Equation cpiUSA -0.337573130 0.36598485 -0.0011856393  1.18867823 0.481646459
      
      $Bup
                              ECT       Const        Trend dolcan t -1  cpiUSA t -1
      Equation dolcan -0.02900588  0.01680242 3.615931e-06   0.1850201 -0.007791395
      Equation cpiUSA  0.33504178 -0.12199902 1.994943e-03  -1.5209576  0.614156675
      
      
      Threshold value[1] -0.1176  0.3236
      Model TVECM with  1  thresholds
      
      $Bdown
                              ECT        Const dolcan t -1 cpiUSA t -1
      Equation dolcan -0.03750343 -0.006852286   0.1250896 -0.01024718
      Equation cpiUSA  0.70847089  0.400683554   0.1907474  0.19818232
      
      $Bup
                              ECT        Const dolcan t -1 cpiUSA t -1
      Equation dolcan 0.006057814 -0.001555704   0.1644276 0.002927573
      Equation cpiUSA 0.058696937  0.085773731   0.2256604 0.619010979
      
      
      Threshold value[1] -0.1176
      Model TVECM with  2  thresholds
      
      $Bdown
                              ECT        Const dolcan t -1 cpiUSA t -1
      Equation dolcan -0.03750343 -0.006852286   0.1250896 -0.01024718
      Equation cpiUSA  0.70847089  0.400683554   0.1907474  0.19818232
      
      $Bmiddle
                              ECT        Const dolcan t -1 cpiUSA t -1
      Equation dolcan  0.04152116 -0.006533316   0.2580703  0.01122303
      Equation cpiUSA -0.35664491  0.184529474   0.4800999  0.27966286
      
      $Bup
                                ECT        Const dolcan t -1 cpiUSA t -1
      Equation dolcan -0.0007475249 0.0007799121   0.1052364 0.002985671
      Equation cpiUSA -0.0967219539 0.1360188778  -0.7272839 0.680080621
      
      
      Threshold value[1] -0.1205  0.1832
      Model TVECM with  1  thresholds
      
      $Bdown
                              ECT dolcan t -1 cpiUSA t -1
      Equation dolcan -0.02202667    0.117653 -0.01446972
      Equation cpiUSA -0.55070230    1.016250  0.32487254
      
      $Bup
                             ECT dolcan t -1 cpiUSA t -1
      Equation dolcan 0.00279744   0.1680402 0.001035012
      Equation cpiUSA 0.18971894  -0.2859286 0.750443210
      
      
      Threshold value[1] -0.1176
      Model TVECM with  2  thresholds
      
      $Bdown
                              ECT dolcan t -1 cpiUSA t -1
      Equation dolcan -0.02522709   0.1674051 -0.01801245
      Equation cpiUSA -0.45272804   2.2623818  0.24296512
      
      $Bmiddle
                              ECT dolcan t -1  cpiUSA t -1
      Equation dolcan -0.01545396  0.07333102 -0.009306737
      Equation cpiUSA -0.92219950 -0.54141665  0.301735719
      
      $Bup
                              ECT dolcan t -1  cpiUSA t -1
      Equation dolcan 0.003049438   0.1689193 0.0007218864
      Equation cpiUSA 0.194655833  -0.2687060 0.7443087565
      
      
      Threshold value[1] -0.2599 -0.1051
      Model TVECM with  1  thresholds
      
      $Bdown
                             ECT         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan -0.0599990 -4.729002e-05   0.2434441  -0.0116853
      Equation cpiUSA -0.5717297  4.019271e-06   1.4443120   0.1073026
      
      $Bup
                              ECT        Trend dolcan t -1  cpiUSA t -1
      Equation dolcan 0.002438058 9.075573e-08   0.1271129 0.0007793217
      Equation cpiUSA 0.188165815 5.901968e-04  -0.1180760 0.5950325178
      
      
      Threshold value[1] -0.2639
      Model TVECM with  2  thresholds
      
      $Bdown
                              ECT         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan -0.06592635 -4.893727e-05   0.1448262  -0.0165236
      Equation cpiUSA -0.74852793 -3.444909e-04   2.2377120   0.2191881
      
      $Bmiddle
                              ECT         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan 0.002513236 -3.076313e-07   0.1907503 0.002161459
      Equation cpiUSA 0.080259228  7.917784e-04  -0.4850491 0.383556172
      
      $Bup
                              ECT         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan 0.004302129 -0.0000163670  0.08396771 0.001005123
      Equation cpiUSA 0.129046087  0.0007655733 -0.67724374 0.707856779
      
      
      Threshold value[1] -0.2796  0.2373
      Model TVECM with  1  thresholds
      
      $Bdown
                              ECT      Const         Trend dolcan t -1  cpiUSA t -1
      Equation dolcan -0.01518849 0.01105132 -4.790127e-05   0.1838534 -0.004546173
      Equation cpiUSA -0.25113994 0.33469736 -8.430391e-04   0.5327942  0.337989478
                      dolcan t -2  cpiUSA t -2
      Equation dolcan  0.02939966 -0.002047247
      Equation cpiUSA -1.22933741 -0.162023742
      
      $Bup
                              ECT       Const         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan -0.01284187  0.00991402 -3.535269e-05    0.125544 0.005828206
      Equation cpiUSA  0.18220945 -0.04146424  1.147708e-03   -1.067956 0.418093241
                      dolcan t -2 cpiUSA t -2
      Equation dolcan  -0.1979610 -0.01062449
      Equation cpiUSA   0.9412034  0.30705642
      
      
      Threshold value[1] 0.2123
      Model TVECM with  2  thresholds
      
      $Bdown
                              ECT      Const         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan -0.01739221 0.01123809 -4.874173e-05   0.1570717 -0.01032538
      Equation cpiUSA -0.18905435 0.28407694 -5.662286e-04   0.4259146  0.28726509
                      dolcan t -2  cpiUSA t -2
      Equation dolcan  0.05235467  0.001209455
      Equation cpiUSA -1.18744412 -0.119469224
      
      $Bmiddle
                             ECT      Const         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan -0.3473239 0.07524714 -5.319414e-05   0.1602594  0.03364773
      Equation cpiUSA -6.1188098 2.70914734 -9.002670e-03  -1.9196425  0.43896490
                      dolcan t -2 cpiUSA t -2
      Equation dolcan   0.1326966 -0.03037366
      Equation cpiUSA  -4.9984845 -0.24567650
      
      $Bup
                             ECT        Const        Trend dolcan t -1  cpiUSA t -1
      Equation dolcan -0.0119483  0.009739977 2.161084e-05   0.1143009 -0.002826679
      Equation cpiUSA  0.3550181 -0.113831430 1.882903e-03  -1.3483856  0.368199325
                      dolcan t -2 cpiUSA t -2
      Equation dolcan  -0.2251498  -0.0106057
      Equation cpiUSA   0.4192614   0.2529322
      
      
      Threshold value[1] 0.1839 0.2161
      Model TVECM with  1  thresholds
      
      $Bdown
                               ECT       Const dolcan t -1  cpiUSA t -1 dolcan t -2
      Equation dolcan -0.002677087 0.001901141   0.1821789 -0.003349176  0.02289789
      Equation cpiUSA -0.030945320 0.173658592   0.5033226  0.359056047 -1.34376545
                        cpiUSA t -2
      Equation dolcan -0.0007533507
      Equation cpiUSA -0.1392518022
      
      $Bup
                               ECT      Const dolcan t -1 cpiUSA t -1 dolcan t -2
      Equation dolcan -0.001946964 0.00412689   0.1173725 0.005439081  -0.2105055
      Equation cpiUSA -0.171488350 0.14641214  -0.8026721 0.430726006   1.3484553
                      cpiUSA t -2
      Equation dolcan -0.01099955
      Equation cpiUSA  0.31923260
      
      
      Threshold value[1] 0.2123
      Model TVECM with  2  thresholds
      
      $Bdown
                             ECT       Const dolcan t -1 cpiUSA t -1 dolcan t -2
      Equation dolcan -0.0461568 -0.01011165   0.1582452 -0.01119553  0.01842045
      Equation cpiUSA  0.9601284  0.51899776  -0.3722760  0.20683381 -1.60201032
                       cpiUSA t -2
      Equation dolcan  0.003495803
      Equation cpiUSA -0.161519009
      
      $Bmiddle
                             ECT        Const dolcan t -1 cpiUSA t -1 dolcan t -2
      Equation dolcan  0.0204315 -0.003717539   0.1867325  0.01718006  0.01977927
      Equation cpiUSA -0.1547078  0.142118177   1.4258148  0.58946004 -1.18049183
                       cpiUSA t -2
      Equation dolcan -0.009962532
      Equation cpiUSA -0.135171254
      
      $Bup
                              ECT      Const dolcan t -1  cpiUSA t -1 dolcan t -2
      Equation dolcan -0.01764468 0.01267829   0.1202763 -0.002199864  -0.2165529
      Equation cpiUSA -0.14129491 0.14217690  -0.8277689  0.422812224   1.1682873
                      cpiUSA t -2
      Equation dolcan -0.01009317
      Equation cpiUSA  0.29758780
      
      
      Threshold value[1] -0.2136  0.2161
      Model TVECM with  1  thresholds
      
      $Bdown
                              ECT dolcan t -1 cpiUSA t -1  dolcan t -2  cpiUSA t -2
      Equation dolcan -0.01791467   0.1391735 -0.01328401  0.002174167  0.001436122
      Equation cpiUSA -0.56534507   0.6575878  0.31447878 -0.623941905 -0.088049070
      
      $Bup
                              ECT  dolcan t -1 cpiUSA t -1 dolcan t -2 cpiUSA t -2
      Equation dolcan 0.006491931  0.160034945 0.007915646 -0.06608634 -0.01065185
      Equation cpiUSA 0.154651435 -0.008908433 0.671157310 -0.91581701  0.13182613
      
      
      Threshold value[1] -0.1898
      Model TVECM with  2  thresholds
      
      $Bdown
                              ECT dolcan t -1  cpiUSA t -1 dolcan t -2  cpiUSA t -2
      Equation dolcan -0.01682907   0.2487917 -0.009049603  -0.1622036 -0.000472127
      Equation cpiUSA -0.63057527   0.7048539  0.113211738   0.7159188 -0.144705191
      
      $Bmiddle
                              ECT dolcan t -1 cpiUSA t -1 dolcan t -2 cpiUSA t -2
      Equation dolcan -0.01504519  0.03394955 -0.02606875   0.2262093  0.01536797
      Equation cpiUSA -0.70956313 -0.33245260  0.43942156  -1.4005139 -0.03286632
      
      $Bup
                              ECT dolcan t -1 cpiUSA t -1 dolcan t -2 cpiUSA t -2
      Equation dolcan 0.005762436   0.1725675  0.01049994 -0.08084231 -0.01247178
      Equation cpiUSA 0.152098278  -0.1468867  0.63963012 -1.13241543  0.16457327
      
      
      Threshold value[1] -0.2834 -0.2136
      Model TVECM with  1  thresholds
      
      $Bdown
                               ECT        Trend dolcan t -1  cpiUSA t -1 dolcan t -2
      Equation dolcan -0.001594506 5.687498e-06   0.1856852 -0.001987102  0.02346587
      Equation cpiUSA  0.160563744 7.799357e-04   0.5882712  0.415492814 -1.40904656
                        cpiUSA t -2
      Equation dolcan  0.0005909138
      Equation cpiUSA -0.0821251317
      
      $Bup
                              ECT        Trend dolcan t -1 cpiUSA t -1 dolcan t -2
      Equation dolcan 0.006034676 7.266438e-06   0.1190289 0.006094024  -0.2116016
      Equation cpiUSA 0.103260490 9.694584e-04  -1.0407069 0.416981490   0.9982535
                       cpiUSA t -2
      Equation dolcan -0.009724057
      Equation cpiUSA  0.303290461
      
      
      Threshold value[1] 0.2123
      Model TVECM with  2  thresholds
      
      $Bdown
                               ECT        Trend dolcan t -1  cpiUSA t -1 dolcan t -2
      Equation dolcan -0.002995578 6.111664e-06   0.1616106 -0.007636568  0.04661039
      Equation cpiUSA  0.174864230 8.203576e-04   0.5406502  0.355232986 -1.33264807
                       cpiUSA t -2
      Equation dolcan  0.004104931
      Equation cpiUSA -0.046277248
      
      $Bmiddle
                              ECT         Trend dolcan t -1 cpiUSA t -1 dolcan t -2
      Equation dolcan -0.05993876  6.599656e-05   0.1751623  0.04194837   0.2039924
      Equation cpiUSA  4.22801254 -4.711408e-03  -1.3830880  0.73781541  -2.4316002
                      cpiUSA t -2
      Equation dolcan -0.03413151
      Equation cpiUSA -0.38097151
      
      $Bup
                              ECT        Trend dolcan t -1  cpiUSA t -1 dolcan t -2
      Equation dolcan 0.006794794 6.522883e-05   0.1066869 -0.002756671  -0.2406763
      Equation cpiUSA 0.135966975 1.373138e-03  -1.2593996  0.367381138   0.6007203
                      cpiUSA t -2
      Equation dolcan -0.01007585
      Equation cpiUSA  0.24673984
      
      
      Threshold value[1] 0.1839 0.2161
      [[1]]
                       Intercept         Trend  dolcan -1   cpiUSA -1
      Equation dolcan 0.01044232 -2.172806e-05  0.9871079 0.000145502
      Equation cpiUSA 0.37899558 -7.098272e-04 -0.1395487 1.002279482
      
      [[2]]
                       Intercept  dolcan -1    cpiUSA -1
      Equation dolcan 0.01278017  0.9870607 6.443813e-05
      Equation cpiUSA 0.45537005 -0.1410917 9.996312e-01
      
      [[3]]
                      dolcan -1    cpiUSA -1
      Equation dolcan 1.0004113 1.200969e-05
      Equation cpiUSA 0.3346056 9.977632e-01
      
      [[4]]
                              Trend  dolcan -1    cpiUSA -1
      Equation dolcan -0.0000766771 0.99195961 0.0003319243
      Equation cpiUSA -0.0027041580 0.03653954 1.0090455285
      
      [[5]]
                       Intercept         Trend dolcan -1    cpiUSA -1  dolcan -2
      Equation dolcan 0.01360877 -2.594234e-05 1.1503909 -0.003742394 -0.1660825
      Equation cpiUSA 0.28320638  3.472133e-04 0.5331984  1.457512609 -0.6049136
                         cpiUSA -2
      Equation dolcan  0.003911853
      Equation cpiUSA -0.459101557
      
      [[6]]
                       Intercept dolcan -1    cpiUSA -1  dolcan -2    cpiUSA -2
      Equation dolcan 0.01636855 1.1506724 -0.003709029 -0.1664027  0.003781728
      Equation cpiUSA 0.24626928 0.5294298  1.457066046 -0.6006284 -0.457359952
      
      [[7]]
                      dolcan -1   cpiUSA -1  dolcan -2    cpiUSA -2
      Equation dolcan 1.1609666 -0.00194731 -0.1602002  0.001956871
      Equation cpiUSA 0.6843077  1.48357158 -0.5073094 -0.484815406
      
      [[8]]
                              Trend dolcan -1    cpiUSA -1  dolcan -2    cpiUSA -2
      Equation dolcan -9.611456e-05 1.1534930 -0.003171405 -0.1628884  0.003578896
      Equation cpiUSA -1.113111e-03 0.5977555  1.469395219 -0.5384425 -0.466030590
      
      [[9]]
                              ECT  Intercept         Trend dolcan -1    cpiUSA -1
      Equation dolcan -0.01506974 0.01098652 -4.515916e-05 0.1655667 -0.003976831
      Equation cpiUSA -0.05033679 0.19306592 -3.133701e-04 0.5871833  0.456867918
      
      [[10]]
                               ECT   Intercept dolcan -1   cpiUSA -1
      Equation dolcan -0.001106223 0.002123837 0.1595354 -0.00300748
      Equation cpiUSA  0.046559351 0.131565673 0.5453308  0.46359448
      
      [[11]]
                                ECT dolcan -1  cpiUSA -1
      Equation dolcan -0.0006768995 0.1681196 0.00242903
      Equation cpiUSA  0.0731547016 1.0770919 0.80037086
      
      [[12]]
                              ECT        Trend dolcan -1    cpiUSA -1
      Equation dolcan 0.001372394 7.603074e-06 0.1610732 -0.001227277
      Equation cpiUSA 0.238600556 6.138198e-04 0.5082187  0.505185783
      
      [[13]]
                              ECT  Intercept         Trend dolcan -1    cpiUSA -1
      Equation dolcan -0.01524527 0.01156071 -4.607701e-05 0.1678529 -0.002542136
      Equation cpiUSA -0.03983784 0.18485483 -2.787031e-04 0.6918440  0.453852695
                        dolcan -2   cpiUSA -2
      Equation dolcan -0.02196011 -0.00308873
      Equation cpiUSA -0.63011337  0.01100919
      
      [[14]]
                               ECT   Intercept dolcan -1    cpiUSA -1   dolcan -2
      Equation dolcan -0.001031649 0.002457958 0.1635127 -0.001892005 -0.03033589
      Equation cpiUSA  0.046135183 0.129795583 0.6655917  0.457785100 -0.68077540
                         cpiUSA -2
      Equation dolcan -0.002298449
      Equation cpiUSA  0.015789314
      
      [[15]]
                                ECT dolcan -1  cpiUSA -1   dolcan -2   cpiUSA -2
      Equation dolcan -0.0007525866  0.172287 0.00152943 -0.02518247 0.001222138
      Equation cpiUSA  0.0608714297  1.128933 0.63845831 -0.40864241 0.201698403
      
      [[16]]
                              ECT        Trend dolcan -1     cpiUSA -1   dolcan -2
      Equation dolcan 0.001657334 8.424941e-06 0.1655519 -0.0006685009 -0.03020354
      Equation cpiUSA 0.230433428 5.927787e-04 0.6550524  0.4838119692 -0.76192515
                         cpiUSA -2
      Equation dolcan -0.001017307
      Equation cpiUSA  0.044131076
      
      [[17]]
      [1] "1.1439"
      
      [[18]]
      [1] "1.1439 1.3893"
      
      [[19]]
      [1] "1.2252"
      
      [[20]]
      [1] "1.1439 1.2252"
      
      [[21]]
      [1] "1.2252"
      
      [[22]]
      [1] "1.1439 1.2252"
      
      [[23]]
      [1] "1.2252"
      
      [[24]]
      [1] "1.1439 1.2252"
      
      [[25]]
      [1] "1.1439"
      
      [[26]]
      [1] "1.1439 1.3893"
      
      [[27]]
      [1] "1.2252"
      
      [[28]]
      [1] "1.2252 1.3893"
      
      [[29]]
      [1] "1.2263"
      
      [[30]]
      [1] "1.1439 1.2252"
      
      [[31]]
      [1] "1.2263"
      
      [[32]]
      [1] "1.2338 1.3893"
      
      [[33]]
      [1] -0.1176
      
      [[34]]
      [1] -0.1176  0.3236
      
      [[35]]
      [1] -0.1176
      
      [[36]]
      [1] -0.1205  0.1832
      
      [[37]]
      [1] -0.1176
      
      [[38]]
      [1] -0.2599 -0.1051
      
      [[39]]
      [1] -0.2639
      
      [[40]]
      [1] -0.2796  0.2373
      
      [[41]]
      [1] 0.2123
      
      [[42]]
      [1] 0.1839 0.2161
      
      [[43]]
      [1] 0.2123
      
      [[44]]
      [1] -0.2136  0.2161
      
      [[45]]
      [1] -0.1898
      
      [[46]]
      [1] -0.2834 -0.2136
      
      [[47]]
      [1] 0.2123
      
      [[48]]
      [1] 0.1839 0.2161
      

---

    Code
      sapply(mods, summary)
    Output
      [[1]]
      #############
      ###Model VAR 
      #############
      Full sample size: 324 	End sample size: 323
      Number of variables: 2 	Number of estimated slope parameters 8
      AIC -3918.837 	BIC -3888.616 	SSR 10.04486
      
                      Intercept          Trend                 dolcan -1          
      Equation dolcan 0.0104(0.0125)     -2.2e-05(8.7e-05)     0.9871(0.0089)***  
      Equation cpiUSA 0.3790(0.1702)*    -0.0007(0.0012)       -0.1395(0.1216)    
                      cpiUSA -1         
      Equation dolcan 0.0001(0.0003)    
      Equation cpiUSA 1.0023(0.0045)*** 
      
      [[2]]
      #############
      ###Model VAR 
      #############
      Full sample size: 324 	End sample size: 323
      Number of variables: 2 	Number of estimated slope parameters 6
      AIC -3922.4 	BIC -3899.734 	SSR 10.05601
      
                      Intercept          dolcan -1           cpiUSA -1           
      Equation dolcan 0.0128(0.0082)     0.9871(0.0089)***   6.4e-05(5.1e-05)    
      Equation cpiUSA 0.4554(0.1120)***  -0.1411(0.1215)     0.9996(0.0007)***   
      
      [[3]]
      #############
      ###Model VAR 
      #############
      Full sample size: 324 	End sample size: 323
      Number of variables: 2 	Number of estimated slope parameters 4
      AIC -3907.237 	BIC -3892.126 	SSR 10.57293
      
                      dolcan -1         cpiUSA -1           
      Equation dolcan 1.0004(0.0024)*** 1.2e-05(3.9e-05)    
      Equation cpiUSA 0.3346(0.0333)*** 0.9978(0.0005)***   
      
      [[4]]
      #############
      ###Model VAR 
      #############
      Full sample size: 324 	End sample size: 323
      Number of variables: 2 	Number of estimated slope parameters 6
      AIC -3916.979 	BIC -3894.313 	SSR 10.20019
      
                      Trend                 dolcan -1          cpiUSA -1         
      Equation dolcan -7.7e-05(5.7e-05)     0.9920(0.0068)***  0.0003(0.0002)    
      Equation cpiUSA -0.0027(0.0008)***    0.0365(0.0929)     1.0090(0.0033)*** 
      
      [[5]]
      #############
      ###Model VAR 
      #############
      Full sample size: 324 	End sample size: 322
      Number of variables: 2 	Number of estimated slope parameters 12
      AIC -3982.323 	BIC -3937.029 	SSR 7.949829
      
                      Intercept          Trend                 dolcan -1         
      Equation dolcan 0.0136(0.0125)     -2.6e-05(8.7e-05)     1.1504(0.0554)*** 
      Equation cpiUSA 0.2832(0.1537).    0.0003(0.0011)        0.5332(0.6837)    
                      cpiUSA -1           dolcan -2           cpiUSA -2         
      Equation dolcan -0.0037(0.0041)     -0.1661(0.0554)**   0.0039(0.0041)    
      Equation cpiUSA 1.4575(0.0502)***   -0.6049(0.6833)     -0.4591(0.0504)***
      
      [[6]]
      #############
      ###Model VAR 
      #############
      Full sample size: 324 	End sample size: 322
      Number of variables: 2 	Number of estimated slope parameters 10
      AIC -3986.131 	BIC -3948.385 	SSR 7.952476
      
                      Intercept         dolcan -1          cpiUSA -1          
      Equation dolcan 0.0164(0.0084).   1.1507(0.0553)***  -0.0037(0.0041)    
      Equation cpiUSA 0.2463(0.1031)*   0.5294(0.6827)     1.4571(0.0501)***  
                      dolcan -2           cpiUSA -2         
      Equation dolcan -0.1664(0.0553)**   0.0038(0.0041)    
      Equation cpiUSA -0.6006(0.6822)     -0.4574(0.0501)***
      
      [[7]]
      #############
      ###Model VAR 
      #############
      Full sample size: 324 	End sample size: 322
      Number of variables: 2 	Number of estimated slope parameters 8
      AIC -3980.289 	BIC -3950.092 	SSR 8.095284
      
                      dolcan -1          cpiUSA -1           dolcan -2          
      Equation dolcan 1.1610(0.0553)***  -0.0019(0.0040)     -0.1602(0.0555)**  
      Equation cpiUSA 0.6843(0.6846)     1.4836(0.0492)***   -0.5073(0.6861)    
                      cpiUSA -2         
      Equation dolcan 0.0020(0.0040)    
      Equation cpiUSA -0.4848(0.0491)***
      
      [[8]]
      #############
      ###Model VAR 
      #############
      Full sample size: 324 	End sample size: 322
      Number of variables: 2 	Number of estimated slope parameters 10
      AIC -3981.552 	BIC -3943.806 	SSR 8.034911
      
                      Trend                 dolcan -1          cpiUSA -1          
      Equation dolcan -9.6e-05(5.8e-05)     1.1535(0.0554)***  -0.0032(0.0040)    
      Equation cpiUSA -0.0011(0.0007)       0.5978(0.6854)     1.4694(0.0499)***  
                      dolcan -2           cpiUSA -2         
      Equation dolcan -0.1629(0.0553)**   0.0036(0.0041)    
      Equation cpiUSA -0.5384(0.6849)     -0.4660(0.0505)***
      
      [[9]]
      #############
      ###Model VECM 
      #############
      Full sample size: 324 	End sample size: 322
      Number of variables: 2 	Number of estimated slope parameters 10
      AIC -3983.832 	BIC -3942.312 	SSR 7.960313
      Cointegrating vector (estimated by 2OLS):
         dolcan     cpiUSA
      r1      1 -0.0158527
      
      
                      ECT                 Intercept         Trend               
      Equation dolcan -0.0151(0.0084).    0.0110(0.0053)*   -4.5e-05(2.6e-05).  
      Equation cpiUSA -0.0503(0.1041)     0.1931(0.0650)**  -0.0003(0.0003)     
                      dolcan -1          cpiUSA -1          
      Equation dolcan 0.1656(0.0553)**   -0.0040(0.0041)    
      Equation cpiUSA 0.5872(0.6821)     0.4569(0.0503)***  
      
      [[10]]
      #############
      ###Model VECM 
      #############
      Full sample size: 324 	End sample size: 322
      Number of variables: 2 	Number of estimated slope parameters 8
      AIC -3983.729 	BIC -3949.758 	SSR 7.984527
      Cointegrating vector (estimated by 2OLS):
         dolcan     cpiUSA
      r1      1 -0.0158527
      
      
                      ECT                 Intercept         dolcan -1         
      Equation dolcan -0.0011(0.0025)     0.0021(0.0013).   0.1595(0.0553)**  
      Equation cpiUSA 0.0466(0.0310)      0.1316(0.0156)*** 0.5453(0.6807)    
                      cpiUSA -1          
      Equation dolcan -0.0030(0.0040)    
      Equation cpiUSA 0.4636(0.0498)***  
      
      [[11]]
      #############
      ###Model VECM 
      #############
      Full sample size: 324 	End sample size: 322
      Number of variables: 2 	Number of estimated slope parameters 6
      AIC -3919.627 	BIC -3893.205 	SSR 9.764497
      Cointegrating vector (estimated by 2OLS):
         dolcan     cpiUSA
      r1      1 -0.0158527
      
      
                      ECT                 dolcan -1          cpiUSA -1         
      Equation dolcan -0.0007(0.0025)     0.1681(0.0552)**   0.0024(0.0024)    
      Equation cpiUSA 0.0732(0.0340)*     1.0771(0.7488)     0.8004(0.0329)*** 
      
      [[12]]
      #############
      ###Model VECM 
      #############
      Full sample size: 324 	End sample size: 322
      Number of variables: 2 	Number of estimated slope parameters 8
      AIC -3974.33 	BIC -3940.359 	SSR 8.181358
      Cointegrating vector (estimated by 2OLS):
         dolcan     cpiUSA
      r1      1 -0.0158527
      
      
                      ECT                Trend                dolcan -1         
      Equation dolcan 0.0014(0.0030)     7.6e-06(6.3e-06)     0.1611(0.0555)**  
      Equation cpiUSA 0.2386(0.0376)***  0.0006(7.8e-05)***   0.5082(0.6900)    
                      cpiUSA -1          
      Equation dolcan -0.0012(0.0039)    
      Equation cpiUSA 0.5052(0.0481)***  
      
      [[13]]
      #############
      ###Model VECM 
      #############
      Full sample size: 324 	End sample size: 321
      Number of variables: 2 	Number of estimated slope parameters 14
      AIC -3962.912 	BIC -3906.34 	SSR 7.937843
      Cointegrating vector (estimated by 2OLS):
         dolcan     cpiUSA
      r1      1 -0.0158527
      
      
                      ECT                 Intercept         Trend               
      Equation dolcan -0.0152(0.0085).    0.0116(0.0054)*   -4.6e-05(2.6e-05).  
      Equation cpiUSA -0.0398(0.1052)     0.1849(0.0667)**  -0.0003(0.0003)     
                      dolcan -1          cpiUSA -1           dolcan -2          
      Equation dolcan 0.1679(0.0562)**   -0.0025(0.0046)     -0.0220(0.0564)    
      Equation cpiUSA 0.6918(0.6934)     0.4539(0.0565)***   -0.6301(0.6962)    
                      cpiUSA -2          
      Equation dolcan -0.0031(0.0046)    
      Equation cpiUSA 0.0110(0.0567)     
      
      [[14]]
      #############
      ###Model VECM 
      #############
      Full sample size: 324 	End sample size: 321
      Number of variables: 2 	Number of estimated slope parameters 12
      AIC -3962.974 	BIC -3913.945 	SSR 7.956742
      Cointegrating vector (estimated by 2OLS):
         dolcan     cpiUSA
      r1      1 -0.0158527
      
      
                      ECT                 Intercept         dolcan -1         
      Equation dolcan -0.0010(0.0025)     0.0025(0.0014).   0.1635(0.0563)**  
      Equation cpiUSA 0.0461(0.0313)      0.1298(0.0175)*** 0.6656(0.6924)    
                      cpiUSA -1           dolcan -2           cpiUSA -2          
      Equation dolcan -0.0019(0.0046)     -0.0303(0.0564)     -0.0023(0.0046)    
      Equation cpiUSA 0.4578(0.0563)***   -0.6808(0.6933)     0.0158(0.0564)     
      
      [[15]]
      #############
      ###Model VECM 
      #############
      Full sample size: 324 	End sample size: 321
      Number of variables: 2 	Number of estimated slope parameters 10
      AIC -3911.881 	BIC -3870.395 	SSR 9.344372
      Cointegrating vector (estimated by 2OLS):
         dolcan     cpiUSA
      r1      1 -0.0158527
      
      
                      ECT                 dolcan -1          cpiUSA -1         
      Equation dolcan -0.0008(0.0025)     0.1723(0.0563)**   0.0015(0.0041)    
      Equation cpiUSA 0.0609(0.0338).     1.1289(0.7465)     0.6385(0.0550)*** 
                      dolcan -2           cpiUSA -2         
      Equation dolcan -0.0252(0.0565)     0.0012(0.0041)    
      Equation cpiUSA -0.4086(0.7495)     0.2017(0.0547)*** 
      
      [[16]]
      #############
      ###Model VECM 
      #############
      Full sample size: 324 	End sample size: 321
      Number of variables: 2 	Number of estimated slope parameters 12
      AIC -3954.238 	BIC -3905.209 	SSR 8.131732
      Cointegrating vector (estimated by 2OLS):
         dolcan     cpiUSA
      r1      1 -0.0158527
      
      
                      ECT                Trend                dolcan -1         
      Equation dolcan 0.0017(0.0032)     8.4e-06(7.0e-06)     0.1656(0.0565)**  
      Equation cpiUSA 0.2304(0.0400)***  0.0006(8.6e-05)***   0.6551(0.7006)    
                      cpiUSA -1           dolcan -2           cpiUSA -2          
      Equation dolcan -0.0007(0.0045)     -0.0302(0.0566)     -0.0010(0.0045)    
      Equation cpiUSA 0.4838(0.0561)***   -0.7619(0.7019)     0.0441(0.0560)     
      
      [[17]]
      Model TVAR with  1  thresholds
      
      Full sample size: 324 	End sample size: 323
      Number of variables: 2 	Number of estimated parameters: 16 + 1
      AIC -3974.104 	BIC -3909.884 	 SSR 8.161472 
      
      [[1]]
                      Intercept        Trend            dolcan -1        
      Equation dolcan 0.0821(0.0789)   0.0007(0.0004)   0.9823(0.0483)***
      Equation cpiUSA -0.5559(0.9804)  -0.0018(0.0056)  0.7181(0.6005)   
                      cpiUSA -1        
      Equation dolcan -0.0024(0.0015)  
      Equation cpiUSA 1.0032(0.0191)***
      
      [[2]]
                      Intercept         Trend             dolcan -1         
      Equation dolcan -0.0135(0.0215)   -0.0002(0.0001)   0.9915(0.0107)*** 
      Equation cpiUSA 2.2021(0.2674)*** 0.0073(0.0015)*** -0.7468(0.1326)***
                      cpiUSA -1        
      Equation dolcan 0.0007(0.0005)   
      Equation cpiUSA 0.9702(0.0059)***
      
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
      
      Threshold value: 1.1439
      Percentage of Observations in each regime: 22.3% 77.7% 
      
      [[18]]
      Model TVAR with  2  thresholds
      
      Full sample size: 324 	End sample size: 323
      Number of variables: 2 	Number of estimated parameters: 24 + 2
      AIC -3974.97 	BIC -3876.751 	 SSR 7.794644 
      
      [[1]]
                      Intercept        Trend            dolcan -1        
      Equation dolcan 0.0821(0.0790)   0.0007(0.0004)   0.9823(0.0484)***
      Equation cpiUSA -0.5559(0.9641)  -0.0018(0.0055)  0.7181(0.5905)   
                      cpiUSA -1        
      Equation dolcan -0.0024(0.0015)  
      Equation cpiUSA 1.0032(0.0188)***
      
      [[2]]
                      Intercept         Trend             dolcan -1         
      Equation dolcan -0.0213(0.0222)   -0.0002(0.0001)   0.9952(0.0124)*** 
      Equation cpiUSA 2.0413(0.2707)*** 0.0064(0.0016)*** -0.6686(0.1518)***
                      cpiUSA -1        
      Equation dolcan 0.0008(0.0005)   
      Equation cpiUSA 0.9732(0.0063)***
      
      [[3]]
                      Intercept        Trend            dolcan -1        
      Equation dolcan 0.3579(0.2290)   0.0014(0.0010)   0.8674(0.0769)***
      Equation cpiUSA -1.2736(2.7959)  -0.0035(0.0126)  0.1600(0.9391)   
                      cpiUSA -1        
      Equation dolcan -0.0056(0.0042)  
      Equation cpiUSA 1.0215(0.0519)***
      
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
      
      Threshold value: 1.1439 1.3893
      Percentage of Observations in each regime: 22.3% 66.9% 10.8% 
      
      [[19]]
      Model TVAR with  1  thresholds
      
      Full sample size: 324 	End sample size: 323
      Number of variables: 2 	Number of estimated parameters: 12 + 1
      AIC -3971.026 	BIC -3921.917 	 SSR 8.348522 
      
      [[1]]
                      Intercept          dolcan -1         cpiUSA -1        
      Equation dolcan 0.0077(0.0171)     0.9940(0.0178)*** 1.4e-05(6.9e-05) 
      Equation cpiUSA -1.0125(0.2138)*** 1.3416(0.2218)*** 0.9969(0.0009)***
      
      [[2]]
                      Intercept       dolcan -1         cpiUSA -1        
      Equation dolcan 0.0305(0.0196)  0.9689(0.0175)*** 0.0001(8.2e-05). 
      Equation cpiUSA 0.4650(0.2452). -0.2413(0.2179)   1.0007(0.0010)***
      
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
      
      Threshold value: 1.2252
      Percentage of Observations in each regime: 50.5% 49.5% 
      
      [[20]]
      Model TVAR with  2  thresholds
      
      Full sample size: 324 	End sample size: 323
      Number of variables: 2 	Number of estimated parameters: 18 + 2
      AIC -3976.778 	BIC -3901.225 	 SSR 7.915594 
      
      [[1]]
                      Intercept        dolcan -1         cpiUSA -1        
      Equation dolcan -0.0298(0.0348)  1.0311(0.0373)*** 8.4e-06(0.0002)  
      Equation cpiUSA -0.2707(0.4248)  0.5938(0.4552)    0.9971(0.0019)***
      
      [[2]]
                      Intercept       dolcan -1         cpiUSA -1        
      Equation dolcan 0.0672(0.0777)  0.9434(0.0654)*** 1.1e-05(8.1e-05) 
      Equation cpiUSA 1.9185(0.9477)* -1.0550(0.7981)   0.9957(0.0010)***
      
      [[3]]
                      Intercept       dolcan -1         cpiUSA -1        
      Equation dolcan 0.0305(0.0197)  0.9689(0.0175)*** 0.0001(8.2e-05). 
      Equation cpiUSA 0.4650(0.2399). -0.2413(0.2131)   1.0007(0.0010)***
      
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
      
      Threshold value: 1.1439 1.2252
      Percentage of Observations in each regime: 22.3% 28.2% 49.5% 
      
      [[21]]
      Model TVAR with  1  thresholds
      
      Full sample size: 324 	End sample size: 323
      Number of variables: 2 	Number of estimated parameters: 8 + 1
      AIC -3950.762 	BIC -3916.763 	 SSR 9.029764 
      
      [[1]]
                      dolcan -1         cpiUSA -1        
      Equation dolcan 1.0018(0.0029)*** -6e-06(5.4e-05)  
      Equation cpiUSA 0.3056(0.0377)*** 0.9995(0.0007)***
      
      [[2]]
                      dolcan -1         cpiUSA -1        
      Equation dolcan 0.9949(0.0049)*** 9.2e-05(7.4e-05) 
      Equation cpiUSA 0.1555(0.0630)*   0.9999(0.0010)***
      
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
      
      Threshold value: 1.2252
      Percentage of Observations in each regime: 50.5% 49.5% 
      
      [[22]]
      Model TVAR with  2  thresholds
      
      Full sample size: 324 	End sample size: 323
      Number of variables: 2 	Number of estimated parameters: 12 + 2
      AIC -3975.972 	BIC -3923.085 	 SSR 8.123147 
      
      [[1]]
                      dolcan -1         cpiUSA -1        
      Equation dolcan 0.9995(0.0051)*** 8.1e-05(0.0001)  
      Equation cpiUSA 0.3063(0.0627)*** 0.9977(0.0016)***
      
      [[2]]
                      dolcan -1         cpiUSA -1        
      Equation dolcan 0.9998(0.0049)*** 1.9e-05(8.1e-05) 
      Equation cpiUSA 0.5560(0.0606)*** 0.9960(0.0010)***
      
      [[3]]
                      dolcan -1         cpiUSA -1        
      Equation dolcan 0.9949(0.0049)*** 9.2e-05(7.4e-05) 
      Equation cpiUSA 0.1555(0.0599)**  0.9999(0.0009)***
      
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
      
      Threshold value: 1.1439 1.2252
      Percentage of Observations in each regime: 22.3% 28.2% 49.5% 
      
      [[23]]
      Model TVAR with  1  thresholds
      
      Full sample size: 324 	End sample size: 323
      Number of variables: 2 	Number of estimated parameters: 12 + 1
      AIC -3954.841 	BIC -3905.732 	 SSR 8.842832 
      
      [[1]]
                      Trend             dolcan -1         cpiUSA -1        
      Equation dolcan  2.1e-05(0.0001)  1.0038(0.0130)*** -8.5e-05(0.0005) 
      Equation cpiUSA 0.0039(0.0017)*   0.6855(0.1673)*** 0.9846(0.0064)***
      
      [[2]]
                      Trend             dolcan -1         cpiUSA -1        
      Equation dolcan -0.0002(9.4e-05)* 0.9703(0.0121)*** 0.0010(0.0004)*  
      Equation cpiUSA -0.0014(0.0012)   -0.0072(0.1564)   1.0059(0.0054)***
      
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
      
      Threshold value: 1.2252
      Percentage of Observations in each regime: 50.5% 49.5% 
      
      [[24]]
      Model TVAR with  2  thresholds
      
      Full sample size: 324 	End sample size: 323
      Number of variables: 2 	Number of estimated parameters: 18 + 2
      AIC -3974.302 	BIC -3898.749 	 SSR 8.06911 
      
      [[1]]
                      Trend           dolcan -1         cpiUSA -1        
      Equation dolcan 0.0003(0.0002)  1.0280(0.0201)*** -0.0010(0.0008)  
      Equation cpiUSA 0.0010(0.0024)  0.4086(0.2491)    0.9938(0.0094)***
      
      [[2]]
                      Trend            dolcan -1         cpiUSA -1        
      Equation dolcan -0.0002(0.0002)  0.9847(0.0175)*** 0.0007(0.0007)   
      Equation cpiUSA 0.0017(0.0024)   0.7017(0.2173)**  0.9899(0.0088)***
      
      [[3]]
                      Trend             dolcan -1         cpiUSA -1        
      Equation dolcan -0.0002(9.4e-05)* 0.9703(0.0121)*** 0.0010(0.0004)*  
      Equation cpiUSA -0.0014(0.0012)   -0.0072(0.1501)   1.0059(0.0051)***
      
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
      
      Threshold value: 1.1439 1.2252
      Percentage of Observations in each regime: 22.3% 28.2% 49.5% 
      
      [[25]]
      Model TVAR with  1  thresholds
      
      Full sample size: 324 	End sample size: 322
      Number of variables: 2 	Number of estimated parameters: 24 + 1
      AIC -3999.979 	BIC -3905.615 	 SSR 7.045296 
      
      [[1]]
                      Intercept        Trend            dolcan -1        
      Equation dolcan 0.0923(0.0822)   0.0007(0.0005)   1.1164(0.1652)***
      Equation cpiUSA -0.3846(0.9608)  -0.0006(0.0055)  -0.0579(1.9317)  
                      cpiUSA -1         dolcan -2        cpiUSA -2       
      Equation dolcan 0.0057(0.0160)    -0.1477(0.1670)  -0.0081(0.0163) 
      Equation cpiUSA 1.2098(0.1872)*** 0.6665(1.9529)   -0.2107(0.1901) 
      
      [[2]]
                      Intercept         Trend             dolcan -1        
      Equation dolcan 0.0006(0.0231)    -0.0001(0.0001)   1.1410(0.0596)***
      Equation cpiUSA 1.6399(0.2701)*** 0.0060(0.0014)*** 0.5722(0.6974)   
                      cpiUSA -1         dolcan -2         cpiUSA -2         
      Equation dolcan -0.0035(0.0047)   -0.1549(0.0597)** 0.0040(0.0046)    
      Equation cpiUSA 1.3343(0.0546)*** -1.1085(0.6980)   -0.3583(0.0533)***
      
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
      
      Threshold value: 1.1439
      Percentage of Observations in each regime: 22% 78% 
      
      [[26]]
      Model TVAR with  2  thresholds
      
      Full sample size: 324 	End sample size: 322
      Number of variables: 2 	Number of estimated parameters: 36 + 2
      AIC -3995.99 	BIC -3852.557 	 SSR 6.794512 
      
      [[1]]
                      Intercept        Trend            dolcan -1        
      Equation dolcan 0.0923(0.0817)   0.0007(0.0005)   1.1164(0.1642)***
      Equation cpiUSA -0.3846(0.9528)  -0.0006(0.0054)  -0.0579(1.9156)  
                      cpiUSA -1         dolcan -2        cpiUSA -2       
      Equation dolcan 0.0057(0.0159)    -0.1477(0.1660)  -0.0081(0.0162) 
      Equation cpiUSA 1.2098(0.1856)*** 0.6665(1.9366)   -0.2107(0.1885) 
      
      [[2]]
                      Intercept         Trend             dolcan -1        
      Equation dolcan -0.0154(0.0238)   -0.0002(0.0001)   1.1806(0.0687)***
      Equation cpiUSA 1.5002(0.2773)*** 0.0054(0.0016)*** 0.5086(0.8014)   
                      cpiUSA -1         dolcan -2         cpiUSA -2         
      Equation dolcan 0.0010(0.0052)    -0.1865(0.0681)** -0.0003(0.0051)   
      Equation cpiUSA 1.3202(0.0602)*** -0.9522(0.7944)   -0.3425(0.0589)***
      
      [[3]]
                      Intercept       Trend           dolcan -1        
      Equation dolcan 0.3044(0.2443)  0.0012(0.0011)  0.9527(0.1378)***
      Equation cpiUSA 1.2815(2.8498)  0.0092(0.0130)  2.0105(1.6081)   
                      cpiUSA -1         dolcan -2        cpiUSA -2       
      Equation dolcan -0.0297(0.0111)** -0.0770(0.1330)  0.0253(0.0114)* 
      Equation cpiUSA 1.3074(0.1297)*** -2.3912(1.5514)  -0.3393(0.1329)*
      
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
      
      Threshold value: 1.1439 1.3893
      Percentage of Observations in each regime: 22% 67.1% 10.9% 
      
      [[27]]
      Model TVAR with  1  thresholds
      
      Full sample size: 324 	End sample size: 322
      Number of variables: 2 	Number of estimated parameters: 20 + 1
      AIC -3997.421 	BIC -3918.155 	 SSR 7.278526 
      
      [[1]]
                      Intercept         dolcan -1         cpiUSA -1        
      Equation dolcan 0.0140(0.0181)    1.1438(0.0902)*** 0.0026(0.0065)   
      Equation cpiUSA -0.6473(0.2154)** 0.8354(1.0711)    1.3948(0.0772)***
                      dolcan -2        cpiUSA -2         
      Equation dolcan -0.1574(0.0911). -0.0025(0.0065)   
      Equation cpiUSA 0.0150(1.0824)   -0.3968(0.0770)***
      
      [[2]]
                      Intercept       dolcan -1         cpiUSA -1        
      Equation dolcan 0.0402(0.0198)* 1.1335(0.0713)*** -0.0102(0.0060). 
      Equation cpiUSA 0.2830(0.2353)  0.4492(0.8467)    1.3110(0.0713)***
                      dolcan -2        cpiUSA -2         
      Equation dolcan -0.1713(0.0702)* 0.0104(0.0060).   
      Equation cpiUSA -0.5822(0.8334)  -0.3106(0.0714)***
      
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
      
      Threshold value: 1.2252
      Percentage of Observations in each regime: 50.3% 49.7% 
      
      [[28]]
      Model TVAR with  2  thresholds
      
      Full sample size: 324 	End sample size: 322
      Number of variables: 2 	Number of estimated parameters: 30 + 2
      AIC -3995.703 	BIC -3874.917 	 SSR 6.972166 
      
      [[1]]
                      Intercept         dolcan -1         cpiUSA -1        
      Equation dolcan 0.0140(0.0181)    1.1438(0.0900)*** 0.0026(0.0065)   
      Equation cpiUSA -0.6473(0.2125)** 0.8354(1.0568)    1.3948(0.0762)***
                      dolcan -2        cpiUSA -2         
      Equation dolcan -0.1574(0.0910). -0.0025(0.0065)   
      Equation cpiUSA 0.0150(1.0678)   -0.3968(0.0759)***
      
      [[2]]
                      Intercept        dolcan -1         cpiUSA -1        
      Equation dolcan 0.0233(0.0295)   1.2136(0.0890)*** -0.0041(0.0072)  
      Equation cpiUSA -0.0503(0.3461)  0.5533(1.0451)    1.2690(0.0850)***
                      dolcan -2         cpiUSA -2        
      Equation dolcan -0.2374(0.0858)** 0.0042(0.0072)   
      Equation cpiUSA -0.3583(1.0068)   -0.2696(0.0851)**
      
      [[3]]
                      Intercept        dolcan -1         cpiUSA -1        
      Equation dolcan 0.0580(0.0681)   0.9631(0.1377)*** -0.0279(0.0110)* 
      Equation cpiUSA -0.6535(0.7991)  2.0919(1.6169)    1.3218(0.1292)***
                      dolcan -2        cpiUSA -2       
      Equation dolcan -0.0270(0.1244)  0.0283(0.0110)* 
      Equation cpiUSA -1.9980(1.4599)  -0.3159(0.1297)*
      
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
      
      Threshold value: 1.2252 1.3893
      Percentage of Observations in each regime: 50.3% 38.8% 10.9% 
      
      [[29]]
      Model TVAR with  1  thresholds
      
      Full sample size: 324 	End sample size: 322
      Number of variables: 2 	Number of estimated parameters: 16 + 1
      AIC -3990.483 	BIC -3926.315 	 SSR 7.512746 
      
      [[1]]
                      dolcan -1         cpiUSA -1         dolcan -2       
      Equation dolcan 1.1453(0.0902)*** -0.0006(0.0059)   -0.1436(0.0905) 
      Equation cpiUSA 0.6625(1.0808)    1.4943(0.0712)*** -0.5054(1.0845) 
                      cpiUSA -2         
      Equation dolcan 0.0005(0.0059)    
      Equation cpiUSA -0.4946(0.0712)***
      
      [[2]]
                      dolcan -1         cpiUSA -1         dolcan -2       
      Equation dolcan 1.1609(0.0705)*** -0.0071(0.0060)   -0.1656(0.0704)*
      Equation cpiUSA 0.5957(0.8445)    1.3145(0.0719)*** -0.4870(0.8432) 
                      cpiUSA -2         
      Equation dolcan 0.0072(0.0060)    
      Equation cpiUSA -0.3146(0.0719)***
      
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
      
      Threshold value: 1.2263
      Percentage of Observations in each regime: 50.9% 49.1% 
      
      [[30]]
      Model TVAR with  2  thresholds
      
      Full sample size: 324 	End sample size: 322
      Number of variables: 2 	Number of estimated parameters: 24 + 2
      AIC -3990.169 	BIC -3892.031 	 SSR 7.141784 
      
      [[1]]
                      dolcan -1         cpiUSA -1         dolcan -2       
      Equation dolcan 1.1964(0.1579)*** 0.0037(0.0159)    -0.1986(0.1590) 
      Equation cpiUSA -0.0369(1.8495)   1.2243(0.1863)*** 0.2958(1.8622)  
                      cpiUSA -2       
      Equation dolcan -0.0036(0.0159) 
      Equation cpiUSA -0.2267(0.1863) 
      
      [[2]]
                      dolcan -1         cpiUSA -1         dolcan -2       
      Equation dolcan 1.1265(0.1110)*** 0.0035(0.0075)    -0.1289(0.1112) 
      Equation cpiUSA 0.9516(1.2997)    1.3785(0.0879)*** -0.6091(1.3020) 
                      cpiUSA -2         
      Equation dolcan -0.0034(0.0075)   
      Equation cpiUSA -0.3809(0.0876)***
      
      [[3]]
                      dolcan -1         cpiUSA -1         dolcan -2       
      Equation dolcan 1.1589(0.0707)*** -0.0080(0.0060)   -0.1628(0.0706)*
      Equation cpiUSA 0.6281(0.8282)    1.3264(0.0697)*** -0.5220(0.8266) 
                      cpiUSA -2         
      Equation dolcan 0.0081(0.0060)    
      Equation cpiUSA -0.3266(0.0697)***
      
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
      
      Threshold value: 1.1439 1.2252
      Percentage of Observations in each regime: 22% 28.3% 49.7% 
      
      [[31]]
      Model TVAR with  1  thresholds
      
      Full sample size: 324 	End sample size: 322
      Number of variables: 2 	Number of estimated parameters: 20 + 1
      AIC -3992.26 	BIC -3912.995 	 SSR 7.415455 
      
      [[1]]
                      Trend             dolcan -1         cpiUSA -1        
      Equation dolcan -1.6e-05(0.0001)  1.1450(0.0897)*** -0.0004(0.0060)  
      Equation cpiUSA 0.0031(0.0016).   0.7035(1.0774)    1.4672(0.0723)***
                      dolcan -2        cpiUSA -2         
      Equation dolcan -0.1450(0.0907)  0.0005(0.0059)    
      Equation cpiUSA -0.2401(1.0893)  -0.4792(0.0714)***
      
      [[2]]
                      Trend             dolcan -1         cpiUSA -1        
      Equation dolcan -0.0002(9.6e-05)* 1.1271(0.0716)*** -0.0092(0.0060)  
      Equation cpiUSA -0.0006(0.0012)   0.5092(0.8599)    1.3092(0.0725)***
                      dolcan -2        cpiUSA -2         
      Equation dolcan -0.1574(0.0701)* 0.0102(0.0061).   
      Equation cpiUSA -0.4662(0.8415)  -0.3068(0.0734)***
      
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
      
      Threshold value: 1.2263
      Percentage of Observations in each regime: 50.9% 49.1% 
      
      [[32]]
      Model TVAR with  2  thresholds
      
      Full sample size: 324 	End sample size: 322
      Number of variables: 2 	Number of estimated parameters: 30 + 2
      AIC -3994.773 	BIC -3873.988 	 SSR 7.025195 
      
      [[1]]
                      Trend             dolcan -1         cpiUSA -1        
      Equation dolcan -6.3e-05(0.0001)  1.1558(0.0889)*** 0.0015(0.0058)   
      Equation cpiUSA 0.0030(0.0015).   1.0534(1.0500)    1.5076(0.0686)***
                      dolcan -2        cpiUSA -2         
      Equation dolcan -0.1604(0.0899). -0.0013(0.0057)   
      Equation cpiUSA -0.6132(1.0612)  -0.5193(0.0677)***
      
      [[2]]
                      Trend            dolcan -1         cpiUSA -1        
      Equation dolcan -0.0002(0.0002)  1.1844(0.0926)*** -0.0057(0.0079)  
      Equation cpiUSA 0.0012(0.0020)   0.3427(1.0931)    1.1578(0.0929)***
                      dolcan -2        cpiUSA -2       
      Equation dolcan -0.2111(0.0877)* 0.0066(0.0079)  
      Equation cpiUSA -0.0438(1.0356)  -0.1633(0.0937).
      
      [[3]]
                      Trend            dolcan -1         cpiUSA -1        
      Equation dolcan -0.0002(0.0003)  0.9859(0.1352)*** -0.0276(0.0110)* 
      Equation cpiUSA 0.0036(0.0037)   2.1499(1.5968)    1.3166(0.1297)***
                      dolcan -2        cpiUSA -2       
      Equation dolcan -0.0234(0.1258)  0.0286(0.0111)* 
      Equation cpiUSA -2.1654(1.4854)  -0.3255(0.1308)*
      
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
      
      Threshold value: 1.2338 1.3893
      Percentage of Observations in each regime: 54.7% 34.5% 10.9% 
      
      [[33]]
      #############
      ###Model TVECM
      #############
      Full sample size: 324 	End sample size: 322
      Number of variables: 2 	Number of estimated parameters 20
      AIC -4002.429 	BIC -3923.163 	SSR 7.298715 
      
      
      Cointegrating vector: (1, - 0.01550732 )
      $Bdown
                      ECT              Const              Trend            
      Equation dolcan -0.0517(0.0694). 0.0053(0.6266)     -3.9e-05(0.2682) 
      Equation cpiUSA 0.7379(0.0310)*  0.6239(2.6e-06)*** -0.0011(0.0121)* 
                      dolcan t -1      cpiUSA t -1     
      Equation dolcan 0.1675(0.0457)*  -0.0122(0.0750).
      Equation cpiUSA -0.1020(0.9190)  0.1395(0.0894). 
      
      $Bup
                      ECT              Const           Trend            
      Equation dolcan -0.0010(0.9412)  0.0036(0.6636)  -3.3e-05(0.3894) 
      Equation cpiUSA 0.0155(0.9259)   0.1013(0.3101)  -3.6e-05(0.9374) 
                      dolcan t -1     cpiUSA t -1       
      Equation dolcan 0.1313(0.0797). 0.0023(0.6488)    
      Equation cpiUSA 0.2371(0.7916)  0.6209(4.7e-21)***
      
      
      Threshold
      Values: -0.1176
      Percentage of Observations in each regime 36.3% 63.7% 
      
      [[34]]
      #############
      ###Model TVECM
      #############
      Full sample size: 324 	End sample size: 322
      Number of variables: 2 	Number of estimated parameters 30
      AIC -4001.215 	BIC -3880.429 	SSR 7.00715 
      
      
      Cointegrating vector: (1, - 0.01550732 )
      $Bdown
                      ECT              Const            Trend            
      Equation dolcan -0.0517(0.0684). 0.0053(0.6253)   -3.9e-05(0.2665) 
      Equation cpiUSA 0.7379(0.0290)*  0.6239(2e-06)*** -0.0011(0.0111)* 
                      dolcan t -1      cpiUSA t -1     
      Equation dolcan 0.1675(0.0450)*  -0.0122(0.0739).
      Equation cpiUSA -0.1020(0.9180)  0.1395(0.0855). 
      
      $Bmiddle
                      ECT              Const           Trend           
      Equation dolcan -0.0023(0.8928)  0.0144(0.2359)  -0.0001(0.0823).
      Equation cpiUSA -0.3376(0.0929). 0.3660(0.0117)* -0.0012(0.0902).
                      dolcan t -1     cpiUSA t -1       
      Equation dolcan 0.0850(0.3605)  0.0035(0.5998)    
      Equation cpiUSA 1.1887(0.2830)  0.4816(4.6e-09)***
      
      $Bup
                      ECT              Const            Trend           
      Equation dolcan -0.0290(0.4313)  0.0168(0.3808)   3.6e-06(0.9649) 
      Equation cpiUSA 0.3350(0.4450)   -0.1220(0.5927)  0.0020(0.0423)* 
                      dolcan t -1      cpiUSA t -1       
      Equation dolcan 0.1850(0.1512)   -0.0078(0.4642)   
      Equation cpiUSA -1.5210(0.3211)  0.6142(1.9e-06)***
      
      
      Threshold
      Values: -0.1176 0.3236
      Percentage of Observations in each regime 36.3% 36.3% 27.3% 
      
      [[35]]
      #############
      ###Model TVECM
      #############
      Full sample size: 324 	End sample size: 322
      Number of variables: 2 	Number of estimated parameters 16
      AIC -4001.207 	BIC -3937.04 	SSR 7.41122 
      
      
      Cointegrating vector: (1, - 0.01665858 )
      $Bdown
                      ECT              Const              dolcan t -1    
      Equation dolcan -0.0375(0.0412)* -0.0069(0.2562)    0.1251(0.1130) 
      Equation cpiUSA 0.7085(0.0014)** 0.4007(6.7e-08)*** 0.1907(0.8404) 
                      cpiUSA t -1     
      Equation dolcan -0.0102(0.1122) 
      Equation cpiUSA 0.1982(0.0109)* 
      
      $Bup
                      ECT             Const            dolcan t -1    
      Equation dolcan 0.0061(0.3280)  -0.0016(0.4998)  0.1644(0.0365)*
      Equation cpiUSA 0.0587(0.4304)  0.0858(0.0021)** 0.2257(0.8107) 
                      cpiUSA t -1       
      Equation dolcan 0.0029(0.5748)    
      Equation cpiUSA 0.6190(3.3e-20)***
      
      
      Threshold
      Values: -0.1176
      Percentage of Observations in each regime 41% 59% 
      
      [[36]]
      #############
      ###Model TVECM
      #############
      Full sample size: 324 	End sample size: 322
      Number of variables: 2 	Number of estimated parameters 24
      AIC -4000.627 	BIC -3902.489 	SSR 7.128824 
      
      
      Cointegrating vector: (1, - 0.01665858 )
      $Bdown
                      ECT              Const              dolcan t -1    
      Equation dolcan -0.0375(0.0409)* -0.0069(0.2556)    0.1251(0.1125) 
      Equation cpiUSA 0.7085(0.0012)** 0.4007(4.6e-08)*** 0.1907(0.8383) 
                      cpiUSA t -1     
      Equation dolcan -0.0102(0.1117) 
      Equation cpiUSA 0.1982(0.0099)**
      
      $Bmiddle
                      ECT              Const             dolcan t -1    
      Equation dolcan 0.0415(0.0375)*  -0.0065(0.1015)   0.2581(0.0660).
      Equation cpiUSA -0.3566(0.1320)  0.1845(0.0001)*** 0.4801(0.7729) 
                      cpiUSA t -1    
      Equation dolcan 0.0112(0.3200) 
      Equation cpiUSA 0.2797(0.0375)*
      
      $Bup
                      ECT              Const           dolcan t -1     
      Equation dolcan -0.0007(0.9464)  0.0008(0.8684)  0.1052(0.2848)  
      Equation cpiUSA -0.0967(0.4647)  0.1360(0.0155)* -0.7273(0.5337) 
                      cpiUSA t -1       
      Equation dolcan 0.0030(0.6310)    
      Equation cpiUSA 0.6801(4.8e-18)***
      
      
      Threshold
      Values: -0.1205 0.1832
      Percentage of Observations in each regime 41% 20.2% 38.8% 
      
      [[37]]
      #############
      ###Model TVECM
      #############
      Full sample size: 324 	End sample size: 322
      Number of variables: 2 	Number of estimated parameters 12
      AIC -3974.743 	BIC -3925.674 	SSR 8.235945 
      
      
      Cointegrating vector: (1, - 0.01616518 )
      $Bdown
                      ECT                 dolcan t -1     cpiUSA t -1       
      Equation dolcan -0.0220(0.0014)**   0.1177(0.1375)  -0.0145(0.0175)*  
      Equation cpiUSA -0.5507(7.7e-10)*** 1.0162(0.3105)  0.3249(3.0e-05)***
      
      $Bup
                      ECT               dolcan t -1      cpiUSA t -1       
      Equation dolcan 0.0028(0.5042)    0.1680(0.0297)*  0.0010(0.8048)    
      Equation cpiUSA 0.1897(0.0004)*** -0.2859(0.7693)  0.7504(1.3e-35)***
      
      
      Threshold
      Values: -0.1176
      Percentage of Observations in each regime 39.4% 60.6% 
      
      [[38]]
      #############
      ###Model TVECM
      #############
      Full sample size: 324 	End sample size: 322
      Number of variables: 2 	Number of estimated parameters 18
      AIC -3980.323 	BIC -3904.831 	SSR 7.768609 
      
      
      Cointegrating vector: (1, - 0.01616518 )
      $Bdown
                      ECT                dolcan t -1     cpiUSA t -1     
      Equation dolcan -0.0252(0.0075)**  0.1674(0.1983)  -0.0180(0.0752).
      Equation cpiUSA -0.4527(0.0001)*** 2.2624(0.1578)  0.2430(0.0513). 
      
      $Bmiddle
                      ECT               dolcan t -1      cpiUSA t -1     
      Equation dolcan -0.0155(0.1904)   0.0733(0.4738)   -0.0093(0.2401) 
      Equation cpiUSA -0.9222(7e-10)*** -0.5414(0.6673)  0.3017(0.0021)**
      
      $Bup
                      ECT               dolcan t -1      cpiUSA t -1       
      Equation dolcan 0.0030(0.4741)    0.1689(0.0296)*  0.0007(0.8663)    
      Equation cpiUSA 0.1947(0.0002)*** -0.2687(0.7777)  0.7443(2.3e-35)***
      
      
      Threshold
      Values: -0.2599 -0.1051
      Percentage of Observations in each regime 20.5% 19.3% 60.2% 
      
      [[39]]
      #############
      ###Model TVECM
      #############
      Full sample size: 324 	End sample size: 322
      Number of variables: 2 	Number of estimated parameters 16
      AIC -3986.448 	BIC -3922.28 	SSR 7.674577 
      
      
      Cointegrating vector: (1, - 0.01626386 )
      $Bdown
                      ECT              Trend             dolcan t -1    
      Equation dolcan -0.0600(0.1237)  -4.7e-05(0.2742)  0.2434(0.0454)*
      Equation cpiUSA -0.5717(0.2277)   4.0e-06(0.9939)  1.4443(0.3282) 
                      cpiUSA t -1     
      Equation dolcan -0.0117(0.2077) 
      Equation cpiUSA 0.1073(0.3416)  
      
      $Bup
                      ECT                Trend              dolcan t -1     
      Equation dolcan 0.0024(0.4726)     9.1e-08(0.9908)    0.1271(0.0430)* 
      Equation cpiUSA 0.1882(7.4e-06)*** 0.0006(1.9e-09)*** -0.1181(0.8768) 
                      cpiUSA t -1       
      Equation dolcan 0.0008(0.8562)    
      Equation cpiUSA 0.5950(2.4e-25)***
      
      
      Threshold
      Values: -0.2639
      Percentage of Observations in each regime 22.7% 77.3% 
      
      [[40]]
      #############
      ###Model TVECM
      #############
      Full sample size: 324 	End sample size: 322
      Number of variables: 2 	Number of estimated parameters 24
      AIC -3985.205 	BIC -3887.067 	SSR 7.313595 
      
      
      Cointegrating vector: (1, - 0.01626386 )
      $Bdown
                      ECT              Trend             dolcan t -1    
      Equation dolcan -0.0659(0.1087)  -4.9e-05(0.2894)  0.1448(0.2877) 
      Equation cpiUSA -0.7485(0.1259)  -0.0003(0.5308)   2.2377(0.1679) 
                      cpiUSA t -1     
      Equation dolcan -0.0165(0.1101) 
      Equation cpiUSA 0.2192(0.0752). 
      
      $Bmiddle
                      ECT             Trend              dolcan t -1     
      Equation dolcan 0.0025(0.6581)  -3.1e-07(0.9718)   0.1908(0.0112)* 
      Equation cpiUSA 0.0803(0.2357)  0.0008(2.5e-13)*** -0.4850(0.5863) 
                      cpiUSA t -1       
      Equation dolcan 0.0022(0.7034)    
      Equation cpiUSA 0.3836(3.1e-08)***
      
      $Bup
                      ECT             Trend             dolcan t -1     
      Equation dolcan 0.0043(0.4257)  -1.6e-05(0.6378)  0.0840(0.4470)  
      Equation cpiUSA 0.1290(0.0453)* 0.0008(0.0651).   -0.6772(0.6064) 
                      cpiUSA t -1       
      Equation dolcan 0.0010(0.8868)    
      Equation cpiUSA 0.7079(1.4e-15)***
      
      
      Threshold
      Values: -0.2796 0.2373
      Percentage of Observations in each regime 18.9% 48.1% 32.9% 
      
      [[41]]
      #############
      ###Model TVECM
      #############
      Full sample size: 324 	End sample size: 321
      Number of variables: 2 	Number of estimated parameters 28
      AIC -3978.469 	BIC -3869.098 	SSR 7.026082 
      
      
      Cointegrating vector: (1, - 0.01662568 )
      $Bdown
                      ECT              Const              Trend            
      Equation dolcan -0.0152(0.1067)  0.0111(0.0875).    -4.8e-05(0.1414) 
      Equation cpiUSA -0.2511(0.0227)* 0.3347(1.2e-05)*** -0.0008(0.0270)* 
                      dolcan t -1      cpiUSA t -1        dolcan t -2     
      Equation dolcan 0.1839(0.0057)** -0.0045(0.4014)    0.0294(0.6485)  
      Equation cpiUSA 0.5328(0.4903)   0.3380(1.7e-07)*** -1.2293(0.1033) 
                      cpiUSA t -2     
      Equation dolcan -0.0020(0.7027) 
      Equation cpiUSA -0.1620(0.0101)*
      
      $Bup
                      ECT              Const            Trend            
      Equation dolcan -0.0128(0.6635)  0.0099(0.5146)   -3.5e-05(0.6501) 
      Equation cpiUSA 0.1822(0.5972)   -0.0415(0.8154)  0.0011(0.2079)   
                      dolcan t -1      cpiUSA t -1       dolcan t -2     
      Equation dolcan 0.1255(0.2732)   0.0058(0.5859)    -0.1980(0.1065) 
      Equation cpiUSA -1.0680(0.4246)  0.4181(0.0009)*** 0.9412(0.5104)  
                      cpiUSA t -2     
      Equation dolcan -0.0106(0.3459) 
      Equation cpiUSA 0.3071(0.0201)* 
      
      
      Threshold
      Values: 0.2123
      Percentage of Observations in each regime 67.6% 32.4% 
      
      [[42]]
      #############
      ###Model TVECM
      #############
      Full sample size: 324 	End sample size: 321
      Number of variables: 2 	Number of estimated parameters 42
      AIC -3987.609 	BIC -3821.666 	SSR 6.526916 
      
      
      Cointegrating vector: (1, - 0.01662568 )
      $Bdown
                      ECT              Const             Trend            
      Equation dolcan -0.0174(0.0938). 0.0112(0.0886).   -4.9e-05(0.1480) 
      Equation cpiUSA -0.1891(0.1132)  0.2841(0.0002)*** -0.0006(0.1441)  
                      dolcan t -1     cpiUSA t -1        dolcan t -2     
      Equation dolcan 0.1571(0.0241)* -0.0103(0.0746).   0.0524(0.4485)  
      Equation cpiUSA 0.4259(0.5934)  0.2873(2.1e-05)*** -1.1874(0.1356) 
                      cpiUSA t -2     
      Equation dolcan 0.0012(0.8347)  
      Equation cpiUSA -0.1195(0.0739).
      
      $Bmiddle
                      ECT              Const             Trend              
      Equation dolcan -0.3473(0.1995)  0.0752(0.2544)    -5.3e-05(0.7407)   
      Equation cpiUSA -6.1188(0.0498)* 2.7091(0.0004)*** -0.0090(1.8e-06)***
                      dolcan t -1      cpiUSA t -1     dolcan t -2     
      Equation dolcan 0.1603(0.5168)   0.0336(0.0424)* 0.1327(0.6105)  
      Equation cpiUSA -1.9196(0.4996)  0.4390(0.0215)* -4.9985(0.0960).
                      cpiUSA t -2     
      Equation dolcan -0.0304(0.0602).
      Equation cpiUSA -0.2457(0.1857) 
      
      $Bup
                      ECT              Const            Trend           
      Equation dolcan -0.0119(0.6924)  0.0097(0.5270)   2.2e-05(0.8067) 
      Equation cpiUSA 0.3550(0.3071)   -0.1138(0.5204)  0.0019(0.0646). 
                      dolcan t -1      cpiUSA t -1      dolcan t -2     
      Equation dolcan 0.1143(0.3207)   -0.0028(0.8039)  -0.2251(0.0693).
      Equation cpiUSA -1.3484(0.3085)  0.3682(0.0052)** 0.4193(0.7681)  
                      cpiUSA t -2     
      Equation dolcan -0.0106(0.3633) 
      Equation cpiUSA 0.2529(0.0600). 
      
      
      Threshold
      Values: 0.1839 0.2161
      Percentage of Observations in each regime 60.7% 7.8% 31.5% 
      
      [[43]]
      #############
      ###Model TVECM
      #############
      Full sample size: 324 	End sample size: 321
      Number of variables: 2 	Number of estimated parameters 24
      AIC -3977.099 	BIC -3882.813 	SSR 7.174833 
      
      
      Cointegrating vector: (1, - 0.01662568 )
      $Bdown
                      ECT              Const              dolcan t -1     
      Equation dolcan -0.0027(0.5060)  0.0019(0.2786)     0.1822(0.0062)**
      Equation cpiUSA -0.0309(0.5132)  0.1737(1.3e-15)*** 0.5033(0.5176)  
                      cpiUSA t -1        dolcan t -2      cpiUSA t -2     
      Equation dolcan -0.0033(0.5319)    0.0229(0.7221)   -0.0008(0.8869) 
      Equation cpiUSA 0.3591(2.7e-08)*** -1.3438(0.0765). -0.1393(0.0259)*
      
      $Bup
                      ECT              Const           dolcan t -1     
      Equation dolcan -0.0019(0.9097)  0.0041(0.6184)  0.1174(0.2998)  
      Equation cpiUSA -0.1715(0.3958)  0.1464(0.1334)  -0.8027(0.5463) 
                      cpiUSA t -1       dolcan t -2      cpiUSA t -2     
      Equation dolcan 0.0054(0.6102)    -0.2105(0.0784). -0.0110(0.3281) 
      Equation cpiUSA 0.4307(0.0007)*** 1.3485(0.3367)   0.3192(0.0162)* 
      
      
      Threshold
      Values: 0.2123
      Percentage of Observations in each regime 67.6% 32.4% 
      
      [[44]]
      #############
      ###Model TVECM
      #############
      Full sample size: 324 	End sample size: 321
      Number of variables: 2 	Number of estimated parameters 36
      AIC -3977.753 	BIC -3834.438 	SSR 6.958434 
      
      
      Cointegrating vector: (1, - 0.01662568 )
      $Bdown
                      ECT              Const              dolcan t -1     
      Equation dolcan -0.0462(0.0843). -0.0101(0.2528)    0.1582(0.0625). 
      Equation cpiUSA 0.9601(0.0026)** 0.5190(1.3e-06)*** -0.3723(0.7116) 
                      cpiUSA t -1      dolcan t -2      cpiUSA t -2     
      Equation dolcan -0.0112(0.0970). 0.0184(0.8242)   0.0035(0.6289)  
      Equation cpiUSA 0.2068(0.0102)*  -1.6020(0.1049)  -0.1615(0.0610).
      
      $Bmiddle
                      ECT              Const              dolcan t -1    
      Equation dolcan 0.0204(0.0670).  -0.0037(0.2006)    0.1867(0.0867).
      Equation cpiUSA -0.1547(0.2425)  0.1421(4.8e-05)*** 1.4258(0.2706) 
                      cpiUSA t -1        dolcan t -2      cpiUSA t -2     
      Equation dolcan 0.0172(0.0432)*    0.0198(0.8526)   -0.0100(0.2130) 
      Equation cpiUSA 0.5895(1.2e-08)*** -1.1805(0.3511)  -0.1352(0.1553) 
      
      $Bup
                      ECT              Const           dolcan t -1     
      Equation dolcan -0.0176(0.3564)  0.0127(0.1859)  0.1203(0.2822)  
      Equation cpiUSA -0.1413(0.5343)  0.1422(0.2120)  -0.8278(0.5333) 
                      cpiUSA t -1      dolcan t -2      cpiUSA t -2     
      Equation dolcan -0.0022(0.8419)  -0.2166(0.0668). -0.0101(0.3764) 
      Equation cpiUSA 0.4228(0.0014)** 1.1683(0.4045)   0.2976(0.0288)* 
      
      
      Threshold
      Values: -0.2136 0.2161
      Percentage of Observations in each regime 37.1% 31.5% 31.5% 
      
      [[45]]
      #############
      ###Model TVECM
      #############
      Full sample size: 324 	End sample size: 321
      Number of variables: 2 	Number of estimated parameters 20
      AIC -3951.727 	BIC -3872.527 	SSR 8.129782 
      
      
      Cointegrating vector: (1, - 0.01649411 )
      $Bdown
                      ECT                 dolcan t -1     cpiUSA t -1      
      Equation dolcan -0.0179(0.0204)*    0.1392(0.0950). -0.0133(0.0418)* 
      Equation cpiUSA -0.5653(1.5e-08)*** 0.6576(0.5319)  0.3145(0.0002)***
                      dolcan t -2      cpiUSA t -2     
      Equation dolcan 0.0022(0.9787)   0.0014(0.8397)  
      Equation cpiUSA -0.6239(0.5450)  -0.0880(0.3272) 
      
      $Bup
                      ECT              dolcan t -1      cpiUSA t -1       
      Equation dolcan 0.0065(0.1213)   0.1600(0.0381)*  0.0079(0.2007)    
      Equation cpiUSA 0.1547(0.0037)** -0.0089(0.9927)  0.6712(4.1e-16)***
                      dolcan t -2      cpiUSA t -2     
      Equation dolcan -0.0661(0.3952)  -0.0107(0.0767).
      Equation cpiUSA -0.9158(0.3516)  0.1318(0.0831). 
      
      
      Threshold
      Values: -0.1898
      Percentage of Observations in each regime 37.4% 62.6% 
      
      [[46]]
      #############
      ###Model TVECM
      #############
      Full sample size: 324 	End sample size: 321
      Number of variables: 2 	Number of estimated parameters 30
      AIC -3957.032 	BIC -3836.346 	SSR 7.778024 
      
      
      Cointegrating vector: (1, - 0.01649411 )
      $Bdown
                      ECT                 dolcan t -1     cpiUSA t -1     
      Equation dolcan -0.0168(0.0755).    0.2488(0.0420)* -0.0090(0.3383) 
      Equation cpiUSA -0.6306(2.3e-07)*** 0.7049(0.6472)  0.1132(0.3429)  
                      dolcan t -2      cpiUSA t -2     
      Equation dolcan -0.1622(0.1711)  -0.0005(0.9569) 
      Equation cpiUSA 0.7159(0.6319)   -0.1447(0.1901) 
      
      $Bmiddle
                      ECT                dolcan t -1      cpiUSA t -1     
      Equation dolcan -0.0150(0.3418)    0.0339(0.7750)   -0.0261(0.0252)*
      Equation cpiUSA -0.7096(0.0004)*** -0.3325(0.8246)  0.4394(0.0029)**
                      dolcan t -2      cpiUSA t -2     
      Equation dolcan 0.2262(0.0625).  0.0154(0.2405)  
      Equation cpiUSA -1.4005(0.3601)  -0.0329(0.8423) 
      
      $Bup
                      ECT              dolcan t -1      cpiUSA t -1     
      Equation dolcan 0.0058(0.1433)   0.1726(0.0215)*  0.0105(0.0650). 
      Equation cpiUSA 0.1521(0.0024)** -0.1469(0.8764)  0.6396(4e-17)***
                      dolcan t -2      cpiUSA t -2     
      Equation dolcan -0.0808(0.2810)  -0.0125(0.0285)*
      Equation cpiUSA -1.1324(0.2319)  0.1646(0.0221)* 
      
      
      Threshold
      Values: -0.2834 -0.2136
      Percentage of Observations in each regime 23.4% 12.5% 64.2% 
      
      [[47]]
      #############
      ###Model TVECM
      #############
      Full sample size: 324 	End sample size: 321
      Number of variables: 2 	Number of estimated parameters 24
      AIC -3962.656 	BIC -3868.37 	SSR 7.476702 
      
      
      Cointegrating vector: (1, - 0.01662568 )
      $Bdown
                      ECT              Trend              dolcan t -1     
      Equation dolcan -0.0016(0.7517)  5.7e-06(0.5204)    0.1857(0.0053)**
      Equation cpiUSA 0.1606(0.0082)** 0.0008(1.7e-12)*** 0.5883(0.4589)  
                      cpiUSA t -1        dolcan t -2      cpiUSA t -2     
      Equation dolcan -0.0020(0.7033)    0.0235(0.7161)   0.0006(0.9086)  
      Equation cpiUSA 0.4155(1.3e-10)*** -1.4090(0.0692). -0.0821(0.1839) 
      
      $Bup
                      ECT             Trend            dolcan t -1     
      Equation dolcan 0.0060(0.2932)  7.3e-06(0.8643)  0.1190(0.2980)  
      Equation cpiUSA 0.1033(0.1338)  0.0010(0.0578).  -1.0407(0.4476) 
                      cpiUSA t -1      dolcan t -2      cpiUSA t -2     
      Equation dolcan 0.0061(0.5695)   -0.2116(0.0807). -0.0097(0.3857) 
      Equation cpiUSA 0.4170(0.0013)** 0.9983(0.4909)   0.3033(0.0245)* 
      
      
      Threshold
      Values: 0.2123
      Percentage of Observations in each regime 67.6% 32.4% 
      
      [[48]]
      #############
      ###Model TVECM
      #############
      Full sample size: 324 	End sample size: 321
      Number of variables: 2 	Number of estimated parameters 36
      AIC -3965.483 	BIC -3822.168 	SSR 7.116788 
      
      
      Cointegrating vector: (1, - 0.01662568 )
      $Bdown
                      ECT              Trend              dolcan t -1    
      Equation dolcan -0.0030(0.6191)  6.1e-06(0.5396)    0.1616(0.0206)*
      Equation cpiUSA 0.1749(0.0154)*  0.0008(2.8e-11)*** 0.5407(0.5141) 
                      cpiUSA t -1        dolcan t -2      cpiUSA t -2     
      Equation dolcan -0.0076(0.1712)    0.0466(0.5005)   0.0041(0.4604)  
      Equation cpiUSA 0.3552(1.7e-07)*** -1.3326(0.1068)  -0.0463(0.4851) 
      
      $Bmiddle
                      ECT               Trend             dolcan t -1     
      Equation dolcan -0.0599(0.5427)   6.6e-05(0.5900)   0.1752(0.4792)  
      Equation cpiUSA 4.2280(0.0004)*** -0.0047(0.0014)** -1.3831(0.6392) 
                      cpiUSA t -1        dolcan t -2      cpiUSA t -2     
      Equation dolcan 0.0419(0.0051)**   0.2040(0.4213)   -0.0341(0.0316)*
      Equation cpiUSA 0.7378(4.1e-05)*** -2.4316(0.4214)  -0.3810(0.0441)*
      
      $Bup
                      ECT             Trend            dolcan t -1     
      Equation dolcan 0.0068(0.2495)  6.5e-05(0.2395)  0.1067(0.3525)  
      Equation cpiUSA 0.1360(0.0537). 0.0014(0.0383)*  -1.2594(0.3573) 
                      cpiUSA t -1      dolcan t -2      cpiUSA t -2     
      Equation dolcan -0.0028(0.8092)  -0.2407(0.0483)* -0.0101(0.3878) 
      Equation cpiUSA 0.3674(0.0073)** 0.6007(0.6783)   0.2467(0.0767). 
      
      
      Threshold
      Values: 0.1839 0.2161
      Percentage of Observations in each regime 60.7% 7.8% 31.5% 
      

---

    Code
      sapply(mods, coef)
    Output
      [[1]]
                       Intercept         Trend  dolcan -1   cpiUSA -1
      Equation dolcan 0.01044232 -2.172806e-05  0.9871079 0.000145502
      Equation cpiUSA 0.37899558 -7.098272e-04 -0.1395487 1.002279482
      
      [[2]]
                       Intercept  dolcan -1    cpiUSA -1
      Equation dolcan 0.01278017  0.9870607 6.443813e-05
      Equation cpiUSA 0.45537005 -0.1410917 9.996312e-01
      
      [[3]]
                      dolcan -1    cpiUSA -1
      Equation dolcan 1.0004113 1.200969e-05
      Equation cpiUSA 0.3346056 9.977632e-01
      
      [[4]]
                              Trend  dolcan -1    cpiUSA -1
      Equation dolcan -0.0000766771 0.99195961 0.0003319243
      Equation cpiUSA -0.0027041580 0.03653954 1.0090455285
      
      [[5]]
                       Intercept         Trend dolcan -1    cpiUSA -1  dolcan -2
      Equation dolcan 0.01360877 -2.594234e-05 1.1503909 -0.003742394 -0.1660825
      Equation cpiUSA 0.28320638  3.472133e-04 0.5331984  1.457512609 -0.6049136
                         cpiUSA -2
      Equation dolcan  0.003911853
      Equation cpiUSA -0.459101557
      
      [[6]]
                       Intercept dolcan -1    cpiUSA -1  dolcan -2    cpiUSA -2
      Equation dolcan 0.01636855 1.1506724 -0.003709029 -0.1664027  0.003781728
      Equation cpiUSA 0.24626928 0.5294298  1.457066046 -0.6006284 -0.457359952
      
      [[7]]
                      dolcan -1   cpiUSA -1  dolcan -2    cpiUSA -2
      Equation dolcan 1.1609666 -0.00194731 -0.1602002  0.001956871
      Equation cpiUSA 0.6843077  1.48357158 -0.5073094 -0.484815406
      
      [[8]]
                              Trend dolcan -1    cpiUSA -1  dolcan -2    cpiUSA -2
      Equation dolcan -9.611456e-05 1.1534930 -0.003171405 -0.1628884  0.003578896
      Equation cpiUSA -1.113111e-03 0.5977555  1.469395219 -0.5384425 -0.466030590
      
      [[9]]
                              ECT  Intercept         Trend dolcan -1    cpiUSA -1
      Equation dolcan -0.01506974 0.01098652 -4.515916e-05 0.1655667 -0.003976831
      Equation cpiUSA -0.05033679 0.19306592 -3.133701e-04 0.5871833  0.456867918
      
      [[10]]
                               ECT   Intercept dolcan -1   cpiUSA -1
      Equation dolcan -0.001106223 0.002123837 0.1595354 -0.00300748
      Equation cpiUSA  0.046559351 0.131565673 0.5453308  0.46359448
      
      [[11]]
                                ECT dolcan -1  cpiUSA -1
      Equation dolcan -0.0006768995 0.1681196 0.00242903
      Equation cpiUSA  0.0731547016 1.0770919 0.80037086
      
      [[12]]
                              ECT        Trend dolcan -1    cpiUSA -1
      Equation dolcan 0.001372394 7.603074e-06 0.1610732 -0.001227277
      Equation cpiUSA 0.238600556 6.138198e-04 0.5082187  0.505185783
      
      [[13]]
                              ECT  Intercept         Trend dolcan -1    cpiUSA -1
      Equation dolcan -0.01524527 0.01156071 -4.607701e-05 0.1678529 -0.002542136
      Equation cpiUSA -0.03983784 0.18485483 -2.787031e-04 0.6918440  0.453852695
                        dolcan -2   cpiUSA -2
      Equation dolcan -0.02196011 -0.00308873
      Equation cpiUSA -0.63011337  0.01100919
      
      [[14]]
                               ECT   Intercept dolcan -1    cpiUSA -1   dolcan -2
      Equation dolcan -0.001031649 0.002457958 0.1635127 -0.001892005 -0.03033589
      Equation cpiUSA  0.046135183 0.129795583 0.6655917  0.457785100 -0.68077540
                         cpiUSA -2
      Equation dolcan -0.002298449
      Equation cpiUSA  0.015789314
      
      [[15]]
                                ECT dolcan -1  cpiUSA -1   dolcan -2   cpiUSA -2
      Equation dolcan -0.0007525866  0.172287 0.00152943 -0.02518247 0.001222138
      Equation cpiUSA  0.0608714297  1.128933 0.63845831 -0.40864241 0.201698403
      
      [[16]]
                              ECT        Trend dolcan -1     cpiUSA -1   dolcan -2
      Equation dolcan 0.001657334 8.424941e-06 0.1655519 -0.0006685009 -0.03020354
      Equation cpiUSA 0.230433428 5.927787e-04 0.6550524  0.4838119692 -0.76192515
                         cpiUSA -2
      Equation dolcan -0.001017307
      Equation cpiUSA  0.044131076
      
      [[17]]
      [[17]]$Bdown
                        Intercept         Trend dolcan -1    cpiUSA -1
      Equation dolcan  0.08210166  0.0007054584 0.9822979 -0.002408619
      Equation cpiUSA -0.55588003 -0.0017977750 0.7181117  1.003209833
      
      [[17]]$Bup
                        Intercept         Trend  dolcan -1    cpiUSA -1
      Equation dolcan -0.01354042 -0.0001653546  0.9914561 0.0007168335
      Equation cpiUSA  2.20206986  0.0073398658 -0.7468369 0.9702240258
      
      
      [[18]]
      [[18]]$Bdown
                        Intercept         Trend dolcan -1    cpiUSA -1
      Equation dolcan  0.08210166  0.0007054584 0.9822979 -0.002408619
      Equation cpiUSA -0.55588003 -0.0017977750 0.7181117  1.003209833
      
      [[18]]$Bmiddle
                      Intercept         Trend  dolcan -1    cpiUSA -1
      Equation dolcan -0.021278 -0.0002016937  0.9952333 0.0008405793
      Equation cpiUSA  2.041256  0.0064145304 -0.6685873 0.9732236367
      
      [[18]]$Bup
                       Intercept        Trend dolcan -1    cpiUSA -1
      Equation dolcan  0.3579121  0.001426282 0.8674328 -0.005605241
      Equation cpiUSA -1.2735570 -0.003504812 0.1599709  1.021491583
      
      
      [[19]]
      [[19]]$Bdown
                         Intercept dolcan -1    cpiUSA -1
      Equation dolcan  0.007691528 0.9939684 1.357332e-05
      Equation cpiUSA -1.012509544 1.3415788 9.969201e-01
      
      [[19]]$Bup
                       Intercept  dolcan -1   cpiUSA -1
      Equation dolcan 0.03045117  0.9689214 0.000145664
      Equation cpiUSA 0.46501765 -0.2412925 1.000711673
      
      
      [[20]]
      [[20]]$Bdown
                        Intercept dolcan -1    cpiUSA -1
      Equation dolcan -0.02978663 1.0310970 8.444228e-06
      Equation cpiUSA -0.27074633 0.5937531 9.970502e-01
      
      [[20]]$Bmiddle
                       Intercept dolcan -1    cpiUSA -1
      Equation dolcan 0.06718762  0.943391 1.089975e-05
      Equation cpiUSA 1.91845343 -1.055018 9.957302e-01
      
      [[20]]$Bup
                       Intercept  dolcan -1   cpiUSA -1
      Equation dolcan 0.03045117  0.9689214 0.000145664
      Equation cpiUSA 0.46501765 -0.2412925 1.000711673
      
      
      [[21]]
      [[21]]$Bdown
                      dolcan -1     cpiUSA -1
      Equation dolcan 1.0018378 -5.982604e-06
      Equation cpiUSA 0.3056463  9.994945e-01
      
      [[21]]$Bup
                      dolcan -1    cpiUSA -1
      Equation dolcan 0.9949031 9.159023e-05
      Equation cpiUSA 0.1554722 9.998859e-01
      
      
      [[22]]
      [[22]]$Bdown
                      dolcan -1    cpiUSA -1
      Equation dolcan 0.9994742 8.089218e-05
      Equation cpiUSA 0.3063167 9.977088e-01
      
      [[22]]$Bmiddle
                      dolcan -1    cpiUSA -1
      Equation dolcan 0.9998115 1.872682e-05
      Equation cpiUSA 0.5559962 9.959537e-01
      
      [[22]]$Bup
                      dolcan -1    cpiUSA -1
      Equation dolcan 0.9949031 9.159023e-05
      Equation cpiUSA 0.1554722 9.998859e-01
      
      
      [[23]]
      [[23]]$Bdown
                             Trend dolcan -1    cpiUSA -1
      Equation dolcan 2.085515e-05 1.0038460 -0.000084888
      Equation cpiUSA 3.944688e-03 0.6854827  0.984569738
      
      [[23]]$Bup
                              Trend    dolcan -1   cpiUSA -1
      Equation dolcan -0.0002075003  0.970279328 0.000997695
      Equation cpiUSA -0.0013705132 -0.007164849 1.005870624
      
      
      [[24]]
      [[24]]$Bdown
                             Trend dolcan -1    cpiUSA -1
      Equation dolcan 0.0002877063 1.0280167 -0.001015842
      Equation cpiUSA 0.0010306707 0.4085665  0.993779851
      
      [[24]]$Bmiddle
                              Trend dolcan -1    cpiUSA -1
      Equation dolcan -0.0001727056 0.9846775 0.0006501501
      Equation cpiUSA  0.0016630729 0.7017294 0.9898733656
      
      [[24]]$Bup
                              Trend    dolcan -1   cpiUSA -1
      Equation dolcan -0.0002075003  0.970279328 0.000997695
      Equation cpiUSA -0.0013705132 -0.007164849 1.005870624
      
      
      [[25]]
      [[25]]$Bdown
                        Intercept         Trend   dolcan -1   cpiUSA -1  dolcan -2
      Equation dolcan  0.09231288  0.0006969576  1.11640771 0.005716144 -0.1476973
      Equation cpiUSA -0.38461556 -0.0006441612 -0.05789139 1.209808395  0.6664581
                         cpiUSA -2
      Equation dolcan -0.008053514
      Equation cpiUSA -0.210683537
      
      [[25]]$Bup
                         Intercept         Trend dolcan -1    cpiUSA -1  dolcan -2
      Equation dolcan 0.0006464748 -0.0001118549  1.141001 -0.003541846 -0.1549285
      Equation cpiUSA 1.6398792639  0.0059730476  0.572203  1.334289684 -1.1084546
                        cpiUSA -2
      Equation dolcan  0.00404909
      Equation cpiUSA -0.35828351
      
      
      [[26]]
      [[26]]$Bdown
                        Intercept         Trend   dolcan -1   cpiUSA -1  dolcan -2
      Equation dolcan  0.09231288  0.0006969576  1.11640771 0.005716144 -0.1476973
      Equation cpiUSA -0.38461556 -0.0006441612 -0.05789139 1.209808395  0.6664581
                         cpiUSA -2
      Equation dolcan -0.008053514
      Equation cpiUSA -0.210683537
      
      [[26]]$Bmiddle
                        Intercept        Trend dolcan -1    cpiUSA -1  dolcan -2
      Equation dolcan -0.01537662 -0.000157434 1.1806500 0.0009972436 -0.1865112
      Equation cpiUSA  1.50018865  0.005426318 0.5086428 1.3201624875 -0.9522052
                          cpiUSA -2
      Equation dolcan -0.0003233887
      Equation cpiUSA -0.3424836706
      
      [[26]]$Bup
                      Intercept       Trend dolcan -1   cpiUSA -1   dolcan -2
      Equation dolcan 0.3044483 0.001169578 0.9527351 -0.02974125 -0.07704594
      Equation cpiUSA 1.2815143 0.009182135 2.0104789  1.30743327 -2.39121154
                        cpiUSA -2
      Equation dolcan  0.02531077
      Equation cpiUSA -0.33933921
      
      
      [[27]]
      [[27]]$Bdown
                        Intercept dolcan -1   cpiUSA -1   dolcan -2   cpiUSA -2
      Equation dolcan  0.01395885 1.1438210 0.002560297 -0.15743623 -0.00252461
      Equation cpiUSA -0.64727526 0.8354151 1.394781489  0.01498696 -0.39677145
      
      [[27]]$Bup
                       Intercept dolcan -1   cpiUSA -1  dolcan -2   cpiUSA -2
      Equation dolcan 0.04017077 1.1335259 -0.01023868 -0.1713070  0.01039812
      Equation cpiUSA 0.28300510 0.4492155  1.31096745 -0.5821929 -0.31059283
      
      
      [[28]]
      [[28]]$Bdown
                        Intercept dolcan -1   cpiUSA -1   dolcan -2   cpiUSA -2
      Equation dolcan  0.01395885 1.1438210 0.002560297 -0.15743623 -0.00252461
      Equation cpiUSA -0.64727526 0.8354151 1.394781489  0.01498696 -0.39677145
      
      [[28]]$Bmiddle
                        Intercept dolcan -1    cpiUSA -1  dolcan -2    cpiUSA -2
      Equation dolcan  0.02333395 1.2135735 -0.004115204 -0.2373632  0.004240793
      Equation cpiUSA -0.05032674 0.5533174  1.269024506 -0.3582838 -0.269630556
      
      [[28]]$Bup
                        Intercept dolcan -1   cpiUSA -1   dolcan -2   cpiUSA -2
      Equation dolcan  0.05798003 0.9631114 -0.02791145 -0.02695616  0.02829938
      Equation cpiUSA -0.65346050 2.0919411  1.32179865 -1.99796630 -0.31587620
      
      
      [[29]]
      [[29]]$Bdown
                      dolcan -1     cpiUSA -1  dolcan -2     cpiUSA -2
      Equation dolcan 1.1452641 -0.0005511544 -0.1435657  0.0005498422
      Equation cpiUSA 0.6624965  1.4943161966 -0.5053663 -0.4946045183
      
      [[29]]$Bup
                      dolcan -1    cpiUSA -1  dolcan -2    cpiUSA -2
      Equation dolcan 1.1609381 -0.007105683 -0.1655764  0.007202736
      Equation cpiUSA 0.5957397  1.314457492 -0.4870101 -0.314595092
      
      
      [[30]]
      [[30]]$Bdown
                        dolcan -1   cpiUSA -1  dolcan -2    cpiUSA -2
      Equation dolcan  1.19639319 0.003697079 -0.1986477 -0.003599167
      Equation cpiUSA -0.03686374 1.224304032  0.2957565 -0.226656499
      
      [[30]]$Bmiddle
                      dolcan -1   cpiUSA -1  dolcan -2    cpiUSA -2
      Equation dolcan 1.1265025 0.003469589 -0.1289219 -0.003432864
      Equation cpiUSA 0.9516258 1.378481433 -0.6090764 -0.380906894
      
      [[30]]$Bup
                      dolcan -1    cpiUSA -1  dolcan -2    cpiUSA -2
      Equation dolcan 1.1589176 -0.008044606 -0.1627678  0.008132882
      Equation cpiUSA 0.6281012  1.326424793 -0.5220339 -0.326551571
      
      
      [[31]]
      [[31]]$Bdown
                              Trend dolcan -1     cpiUSA -1  dolcan -2     cpiUSA -2
      Equation dolcan -1.622245e-05 1.1450492 -0.0004092175 -0.1449554  0.0004693694
      Equation cpiUSA  3.096566e-03 0.7035062  1.4672230481 -0.2400976 -0.4792437363
      
      [[31]]$Bup
                              Trend dolcan -1    cpiUSA -1  dolcan -2   cpiUSA -2
      Equation dolcan -0.0002216198 1.1270540 -0.009180805 -0.1574311  0.01024224
      Equation cpiUSA -0.0005662404 0.5091657  1.309155537 -0.4661987 -0.30682914
      
      
      [[32]]
      [[32]]$Bdown
                              Trend dolcan -1   cpiUSA -1  dolcan -2    cpiUSA -2
      Equation dolcan -6.292985e-05  1.155769 0.001521385 -0.1604489 -0.001290011
      Equation cpiUSA  3.007284e-03  1.053381 1.507577719 -0.6132165 -0.519312189
      
      [[32]]$Bmiddle
                              Trend dolcan -1    cpiUSA -1   dolcan -2   cpiUSA -2
      Equation dolcan -0.0001870489 1.1844185 -0.005661468 -0.21109514  0.00657114
      Equation cpiUSA  0.0011972003 0.3427475  1.157796683 -0.04378136 -0.16328716
      
      [[32]]$Bup
                              Trend dolcan -1   cpiUSA -1  dolcan -2   cpiUSA -2
      Equation dolcan -0.0001633523 0.9858677 -0.02757334 -0.0233897  0.02860914
      Equation cpiUSA  0.0035714290 2.1499439  1.31655866 -2.1653563 -0.32545538
      
      
      [[33]]
      [[33]]$Bdown
                              ECT       Const         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan -0.05165772 0.005288736 -3.885903e-05   0.1675227  -0.0121845
      Equation cpiUSA  0.73785216 0.623866396 -1.062156e-03  -0.1019883   0.1395116
      
      [[33]]$Bup
                               ECT       Const         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan -0.001023032 0.003613178 -3.313512e-05   0.1312604 0.002326573
      Equation cpiUSA  0.015485020 0.101293958 -3.625345e-05   0.2371179 0.620902189
      
      
      [[34]]
      [[34]]$Bdown
                              ECT       Const         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan -0.05165772 0.005288736 -3.885903e-05   0.1675227  -0.0121845
      Equation cpiUSA  0.73785216 0.623866396 -1.062156e-03  -0.1019883   0.1395116
      
      [[34]]$Bmiddle
                               ECT      Const         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan -0.002269518 0.01440822 -0.0001021736  0.08504187 0.003521756
      Equation cpiUSA -0.337573130 0.36598485 -0.0011856393  1.18867823 0.481646459
      
      [[34]]$Bup
                              ECT       Const        Trend dolcan t -1  cpiUSA t -1
      Equation dolcan -0.02900588  0.01680242 3.615931e-06   0.1850201 -0.007791395
      Equation cpiUSA  0.33504178 -0.12199902 1.994943e-03  -1.5209576  0.614156675
      
      
      [[35]]
      [[35]]$Bdown
                              ECT        Const dolcan t -1 cpiUSA t -1
      Equation dolcan -0.03750343 -0.006852286   0.1250896 -0.01024718
      Equation cpiUSA  0.70847089  0.400683554   0.1907474  0.19818232
      
      [[35]]$Bup
                              ECT        Const dolcan t -1 cpiUSA t -1
      Equation dolcan 0.006057814 -0.001555704   0.1644276 0.002927573
      Equation cpiUSA 0.058696937  0.085773731   0.2256604 0.619010979
      
      
      [[36]]
      [[36]]$Bdown
                              ECT        Const dolcan t -1 cpiUSA t -1
      Equation dolcan -0.03750343 -0.006852286   0.1250896 -0.01024718
      Equation cpiUSA  0.70847089  0.400683554   0.1907474  0.19818232
      
      [[36]]$Bmiddle
                              ECT        Const dolcan t -1 cpiUSA t -1
      Equation dolcan  0.04152116 -0.006533316   0.2580703  0.01122303
      Equation cpiUSA -0.35664491  0.184529474   0.4800999  0.27966286
      
      [[36]]$Bup
                                ECT        Const dolcan t -1 cpiUSA t -1
      Equation dolcan -0.0007475249 0.0007799121   0.1052364 0.002985671
      Equation cpiUSA -0.0967219539 0.1360188778  -0.7272839 0.680080621
      
      
      [[37]]
      [[37]]$Bdown
                              ECT dolcan t -1 cpiUSA t -1
      Equation dolcan -0.02202667    0.117653 -0.01446972
      Equation cpiUSA -0.55070230    1.016250  0.32487254
      
      [[37]]$Bup
                             ECT dolcan t -1 cpiUSA t -1
      Equation dolcan 0.00279744   0.1680402 0.001035012
      Equation cpiUSA 0.18971894  -0.2859286 0.750443210
      
      
      [[38]]
      [[38]]$Bdown
                              ECT dolcan t -1 cpiUSA t -1
      Equation dolcan -0.02522709   0.1674051 -0.01801245
      Equation cpiUSA -0.45272804   2.2623818  0.24296512
      
      [[38]]$Bmiddle
                              ECT dolcan t -1  cpiUSA t -1
      Equation dolcan -0.01545396  0.07333102 -0.009306737
      Equation cpiUSA -0.92219950 -0.54141665  0.301735719
      
      [[38]]$Bup
                              ECT dolcan t -1  cpiUSA t -1
      Equation dolcan 0.003049438   0.1689193 0.0007218864
      Equation cpiUSA 0.194655833  -0.2687060 0.7443087565
      
      
      [[39]]
      [[39]]$Bdown
                             ECT         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan -0.0599990 -4.729002e-05   0.2434441  -0.0116853
      Equation cpiUSA -0.5717297  4.019271e-06   1.4443120   0.1073026
      
      [[39]]$Bup
                              ECT        Trend dolcan t -1  cpiUSA t -1
      Equation dolcan 0.002438058 9.075573e-08   0.1271129 0.0007793217
      Equation cpiUSA 0.188165815 5.901968e-04  -0.1180760 0.5950325178
      
      
      [[40]]
      [[40]]$Bdown
                              ECT         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan -0.06592635 -4.893727e-05   0.1448262  -0.0165236
      Equation cpiUSA -0.74852793 -3.444909e-04   2.2377120   0.2191881
      
      [[40]]$Bmiddle
                              ECT         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan 0.002513236 -3.076313e-07   0.1907503 0.002161459
      Equation cpiUSA 0.080259228  7.917784e-04  -0.4850491 0.383556172
      
      [[40]]$Bup
                              ECT         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan 0.004302129 -0.0000163670  0.08396771 0.001005123
      Equation cpiUSA 0.129046087  0.0007655733 -0.67724374 0.707856779
      
      
      [[41]]
      [[41]]$Bdown
                              ECT      Const         Trend dolcan t -1  cpiUSA t -1
      Equation dolcan -0.01518849 0.01105132 -4.790127e-05   0.1838534 -0.004546173
      Equation cpiUSA -0.25113994 0.33469736 -8.430391e-04   0.5327942  0.337989478
                      dolcan t -2  cpiUSA t -2
      Equation dolcan  0.02939966 -0.002047247
      Equation cpiUSA -1.22933741 -0.162023742
      
      [[41]]$Bup
                              ECT       Const         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan -0.01284187  0.00991402 -3.535269e-05    0.125544 0.005828206
      Equation cpiUSA  0.18220945 -0.04146424  1.147708e-03   -1.067956 0.418093241
                      dolcan t -2 cpiUSA t -2
      Equation dolcan  -0.1979610 -0.01062449
      Equation cpiUSA   0.9412034  0.30705642
      
      
      [[42]]
      [[42]]$Bdown
                              ECT      Const         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan -0.01739221 0.01123809 -4.874173e-05   0.1570717 -0.01032538
      Equation cpiUSA -0.18905435 0.28407694 -5.662286e-04   0.4259146  0.28726509
                      dolcan t -2  cpiUSA t -2
      Equation dolcan  0.05235467  0.001209455
      Equation cpiUSA -1.18744412 -0.119469224
      
      [[42]]$Bmiddle
                             ECT      Const         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan -0.3473239 0.07524714 -5.319414e-05   0.1602594  0.03364773
      Equation cpiUSA -6.1188098 2.70914734 -9.002670e-03  -1.9196425  0.43896490
                      dolcan t -2 cpiUSA t -2
      Equation dolcan   0.1326966 -0.03037366
      Equation cpiUSA  -4.9984845 -0.24567650
      
      [[42]]$Bup
                             ECT        Const        Trend dolcan t -1  cpiUSA t -1
      Equation dolcan -0.0119483  0.009739977 2.161084e-05   0.1143009 -0.002826679
      Equation cpiUSA  0.3550181 -0.113831430 1.882903e-03  -1.3483856  0.368199325
                      dolcan t -2 cpiUSA t -2
      Equation dolcan  -0.2251498  -0.0106057
      Equation cpiUSA   0.4192614   0.2529322
      
      
      [[43]]
      [[43]]$Bdown
                               ECT       Const dolcan t -1  cpiUSA t -1 dolcan t -2
      Equation dolcan -0.002677087 0.001901141   0.1821789 -0.003349176  0.02289789
      Equation cpiUSA -0.030945320 0.173658592   0.5033226  0.359056047 -1.34376545
                        cpiUSA t -2
      Equation dolcan -0.0007533507
      Equation cpiUSA -0.1392518022
      
      [[43]]$Bup
                               ECT      Const dolcan t -1 cpiUSA t -1 dolcan t -2
      Equation dolcan -0.001946964 0.00412689   0.1173725 0.005439081  -0.2105055
      Equation cpiUSA -0.171488350 0.14641214  -0.8026721 0.430726006   1.3484553
                      cpiUSA t -2
      Equation dolcan -0.01099955
      Equation cpiUSA  0.31923260
      
      
      [[44]]
      [[44]]$Bdown
                             ECT       Const dolcan t -1 cpiUSA t -1 dolcan t -2
      Equation dolcan -0.0461568 -0.01011165   0.1582452 -0.01119553  0.01842045
      Equation cpiUSA  0.9601284  0.51899776  -0.3722760  0.20683381 -1.60201032
                       cpiUSA t -2
      Equation dolcan  0.003495803
      Equation cpiUSA -0.161519009
      
      [[44]]$Bmiddle
                             ECT        Const dolcan t -1 cpiUSA t -1 dolcan t -2
      Equation dolcan  0.0204315 -0.003717539   0.1867325  0.01718006  0.01977927
      Equation cpiUSA -0.1547078  0.142118177   1.4258148  0.58946004 -1.18049183
                       cpiUSA t -2
      Equation dolcan -0.009962532
      Equation cpiUSA -0.135171254
      
      [[44]]$Bup
                              ECT      Const dolcan t -1  cpiUSA t -1 dolcan t -2
      Equation dolcan -0.01764468 0.01267829   0.1202763 -0.002199864  -0.2165529
      Equation cpiUSA -0.14129491 0.14217690  -0.8277689  0.422812224   1.1682873
                      cpiUSA t -2
      Equation dolcan -0.01009317
      Equation cpiUSA  0.29758780
      
      
      [[45]]
      [[45]]$Bdown
                              ECT dolcan t -1 cpiUSA t -1  dolcan t -2  cpiUSA t -2
      Equation dolcan -0.01791467   0.1391735 -0.01328401  0.002174167  0.001436122
      Equation cpiUSA -0.56534507   0.6575878  0.31447878 -0.623941905 -0.088049070
      
      [[45]]$Bup
                              ECT  dolcan t -1 cpiUSA t -1 dolcan t -2 cpiUSA t -2
      Equation dolcan 0.006491931  0.160034945 0.007915646 -0.06608634 -0.01065185
      Equation cpiUSA 0.154651435 -0.008908433 0.671157310 -0.91581701  0.13182613
      
      
      [[46]]
      [[46]]$Bdown
                              ECT dolcan t -1  cpiUSA t -1 dolcan t -2  cpiUSA t -2
      Equation dolcan -0.01682907   0.2487917 -0.009049603  -0.1622036 -0.000472127
      Equation cpiUSA -0.63057527   0.7048539  0.113211738   0.7159188 -0.144705191
      
      [[46]]$Bmiddle
                              ECT dolcan t -1 cpiUSA t -1 dolcan t -2 cpiUSA t -2
      Equation dolcan -0.01504519  0.03394955 -0.02606875   0.2262093  0.01536797
      Equation cpiUSA -0.70956313 -0.33245260  0.43942156  -1.4005139 -0.03286632
      
      [[46]]$Bup
                              ECT dolcan t -1 cpiUSA t -1 dolcan t -2 cpiUSA t -2
      Equation dolcan 0.005762436   0.1725675  0.01049994 -0.08084231 -0.01247178
      Equation cpiUSA 0.152098278  -0.1468867  0.63963012 -1.13241543  0.16457327
      
      
      [[47]]
      [[47]]$Bdown
                               ECT        Trend dolcan t -1  cpiUSA t -1 dolcan t -2
      Equation dolcan -0.001594506 5.687498e-06   0.1856852 -0.001987102  0.02346587
      Equation cpiUSA  0.160563744 7.799357e-04   0.5882712  0.415492814 -1.40904656
                        cpiUSA t -2
      Equation dolcan  0.0005909138
      Equation cpiUSA -0.0821251317
      
      [[47]]$Bup
                              ECT        Trend dolcan t -1 cpiUSA t -1 dolcan t -2
      Equation dolcan 0.006034676 7.266438e-06   0.1190289 0.006094024  -0.2116016
      Equation cpiUSA 0.103260490 9.694584e-04  -1.0407069 0.416981490   0.9982535
                       cpiUSA t -2
      Equation dolcan -0.009724057
      Equation cpiUSA  0.303290461
      
      
      [[48]]
      [[48]]$Bdown
                               ECT        Trend dolcan t -1  cpiUSA t -1 dolcan t -2
      Equation dolcan -0.002995578 6.111664e-06   0.1616106 -0.007636568  0.04661039
      Equation cpiUSA  0.174864230 8.203576e-04   0.5406502  0.355232986 -1.33264807
                       cpiUSA t -2
      Equation dolcan  0.004104931
      Equation cpiUSA -0.046277248
      
      [[48]]$Bmiddle
                              ECT         Trend dolcan t -1 cpiUSA t -1 dolcan t -2
      Equation dolcan -0.05993876  6.599656e-05   0.1751623  0.04194837   0.2039924
      Equation cpiUSA  4.22801254 -4.711408e-03  -1.3830880  0.73781541  -2.4316002
                      cpiUSA t -2
      Equation dolcan -0.03413151
      Equation cpiUSA -0.38097151
      
      [[48]]$Bup
                              ECT        Trend dolcan t -1  cpiUSA t -1 dolcan t -2
      Equation dolcan 0.006794794 6.522883e-05   0.1066869 -0.002756671  -0.2406763
      Equation cpiUSA 0.135966975 1.373138e-03  -1.2593996  0.367381138   0.6007203
                      cpiUSA t -2
      Equation dolcan -0.01007585
      Equation cpiUSA  0.24673984
      
      

---

    Code
      sapply(mods, tsDyn:::coefMat.nlVar)
    Output
      [[1]]
                       Intercept         Trend  dolcan -1   cpiUSA -1
      Equation dolcan 0.01044232 -2.172806e-05  0.9871079 0.000145502
      Equation cpiUSA 0.37899558 -7.098272e-04 -0.1395487 1.002279482
      
      [[2]]
                       Intercept  dolcan -1    cpiUSA -1
      Equation dolcan 0.01278017  0.9870607 6.443813e-05
      Equation cpiUSA 0.45537005 -0.1410917 9.996312e-01
      
      [[3]]
                      dolcan -1    cpiUSA -1
      Equation dolcan 1.0004113 1.200969e-05
      Equation cpiUSA 0.3346056 9.977632e-01
      
      [[4]]
                              Trend  dolcan -1    cpiUSA -1
      Equation dolcan -0.0000766771 0.99195961 0.0003319243
      Equation cpiUSA -0.0027041580 0.03653954 1.0090455285
      
      [[5]]
                       Intercept         Trend dolcan -1    cpiUSA -1  dolcan -2
      Equation dolcan 0.01360877 -2.594234e-05 1.1503909 -0.003742394 -0.1660825
      Equation cpiUSA 0.28320638  3.472133e-04 0.5331984  1.457512609 -0.6049136
                         cpiUSA -2
      Equation dolcan  0.003911853
      Equation cpiUSA -0.459101557
      
      [[6]]
                       Intercept dolcan -1    cpiUSA -1  dolcan -2    cpiUSA -2
      Equation dolcan 0.01636855 1.1506724 -0.003709029 -0.1664027  0.003781728
      Equation cpiUSA 0.24626928 0.5294298  1.457066046 -0.6006284 -0.457359952
      
      [[7]]
                      dolcan -1   cpiUSA -1  dolcan -2    cpiUSA -2
      Equation dolcan 1.1609666 -0.00194731 -0.1602002  0.001956871
      Equation cpiUSA 0.6843077  1.48357158 -0.5073094 -0.484815406
      
      [[8]]
                              Trend dolcan -1    cpiUSA -1  dolcan -2    cpiUSA -2
      Equation dolcan -9.611456e-05 1.1534930 -0.003171405 -0.1628884  0.003578896
      Equation cpiUSA -1.113111e-03 0.5977555  1.469395219 -0.5384425 -0.466030590
      
      [[9]]
                              ECT  Intercept         Trend dolcan -1    cpiUSA -1
      Equation dolcan -0.01506974 0.01098652 -4.515916e-05 0.1655667 -0.003976831
      Equation cpiUSA -0.05033679 0.19306592 -3.133701e-04 0.5871833  0.456867918
      
      [[10]]
                               ECT   Intercept dolcan -1   cpiUSA -1
      Equation dolcan -0.001106223 0.002123837 0.1595354 -0.00300748
      Equation cpiUSA  0.046559351 0.131565673 0.5453308  0.46359448
      
      [[11]]
                                ECT dolcan -1  cpiUSA -1
      Equation dolcan -0.0006768995 0.1681196 0.00242903
      Equation cpiUSA  0.0731547016 1.0770919 0.80037086
      
      [[12]]
                              ECT        Trend dolcan -1    cpiUSA -1
      Equation dolcan 0.001372394 7.603074e-06 0.1610732 -0.001227277
      Equation cpiUSA 0.238600556 6.138198e-04 0.5082187  0.505185783
      
      [[13]]
                              ECT  Intercept         Trend dolcan -1    cpiUSA -1
      Equation dolcan -0.01524527 0.01156071 -4.607701e-05 0.1678529 -0.002542136
      Equation cpiUSA -0.03983784 0.18485483 -2.787031e-04 0.6918440  0.453852695
                        dolcan -2   cpiUSA -2
      Equation dolcan -0.02196011 -0.00308873
      Equation cpiUSA -0.63011337  0.01100919
      
      [[14]]
                               ECT   Intercept dolcan -1    cpiUSA -1   dolcan -2
      Equation dolcan -0.001031649 0.002457958 0.1635127 -0.001892005 -0.03033589
      Equation cpiUSA  0.046135183 0.129795583 0.6655917  0.457785100 -0.68077540
                         cpiUSA -2
      Equation dolcan -0.002298449
      Equation cpiUSA  0.015789314
      
      [[15]]
                                ECT dolcan -1  cpiUSA -1   dolcan -2   cpiUSA -2
      Equation dolcan -0.0007525866  0.172287 0.00152943 -0.02518247 0.001222138
      Equation cpiUSA  0.0608714297  1.128933 0.63845831 -0.40864241 0.201698403
      
      [[16]]
                              ECT        Trend dolcan -1     cpiUSA -1   dolcan -2
      Equation dolcan 0.001657334 8.424941e-06 0.1655519 -0.0006685009 -0.03020354
      Equation cpiUSA 0.230433428 5.927787e-04 0.6550524  0.4838119692 -0.76192515
                         cpiUSA -2
      Equation dolcan -0.001017307
      Equation cpiUSA  0.044131076
      
      [[17]]
                        Intercept         Trend dolcan -1    cpiUSA -1   Intercept
      Equation dolcan  0.08210166  0.0007054584 0.9822979 -0.002408619 -0.01354042
      Equation cpiUSA -0.55588003 -0.0017977750 0.7181117  1.003209833  2.20206986
                              Trend  dolcan -1    cpiUSA -1
      Equation dolcan -0.0001653546  0.9914561 0.0007168335
      Equation cpiUSA  0.0073398658 -0.7468369 0.9702240258
      
      [[18]]
                        Intercept         Trend dolcan -1    cpiUSA -1 Intercept
      Equation dolcan  0.08210166  0.0007054584 0.9822979 -0.002408619 -0.021278
      Equation cpiUSA -0.55588003 -0.0017977750 0.7181117  1.003209833  2.041256
                              Trend  dolcan -1    cpiUSA -1  Intercept        Trend
      Equation dolcan -0.0002016937  0.9952333 0.0008405793  0.3579121  0.001426282
      Equation cpiUSA  0.0064145304 -0.6685873 0.9732236367 -1.2735570 -0.003504812
                      dolcan -1    cpiUSA -1
      Equation dolcan 0.8674328 -0.005605241
      Equation cpiUSA 0.1599709  1.021491583
      
      [[19]]
                         Intercept dolcan -1    cpiUSA -1  Intercept  dolcan -1
      Equation dolcan  0.007691528 0.9939684 1.357332e-05 0.03045117  0.9689214
      Equation cpiUSA -1.012509544 1.3415788 9.969201e-01 0.46501765 -0.2412925
                        cpiUSA -1
      Equation dolcan 0.000145664
      Equation cpiUSA 1.000711673
      
      [[20]]
                        Intercept dolcan -1    cpiUSA -1  Intercept dolcan -1
      Equation dolcan -0.02978663 1.0310970 8.444228e-06 0.06718762  0.943391
      Equation cpiUSA -0.27074633 0.5937531 9.970502e-01 1.91845343 -1.055018
                         cpiUSA -1  Intercept  dolcan -1   cpiUSA -1
      Equation dolcan 1.089975e-05 0.03045117  0.9689214 0.000145664
      Equation cpiUSA 9.957302e-01 0.46501765 -0.2412925 1.000711673
      
      [[21]]
                      dolcan -1     cpiUSA -1 dolcan -1    cpiUSA -1
      Equation dolcan 1.0018378 -5.982604e-06 0.9949031 9.159023e-05
      Equation cpiUSA 0.3056463  9.994945e-01 0.1554722 9.998859e-01
      
      [[22]]
                      dolcan -1    cpiUSA -1 dolcan -1    cpiUSA -1 dolcan -1
      Equation dolcan 0.9994742 8.089218e-05 0.9998115 1.872682e-05 0.9949031
      Equation cpiUSA 0.3063167 9.977088e-01 0.5559962 9.959537e-01 0.1554722
                         cpiUSA -1
      Equation dolcan 9.159023e-05
      Equation cpiUSA 9.998859e-01
      
      [[23]]
                             Trend dolcan -1    cpiUSA -1         Trend    dolcan -1
      Equation dolcan 2.085515e-05 1.0038460 -0.000084888 -0.0002075003  0.970279328
      Equation cpiUSA 3.944688e-03 0.6854827  0.984569738 -0.0013705132 -0.007164849
                        cpiUSA -1
      Equation dolcan 0.000997695
      Equation cpiUSA 1.005870624
      
      [[24]]
                             Trend dolcan -1    cpiUSA -1         Trend dolcan -1
      Equation dolcan 0.0002877063 1.0280167 -0.001015842 -0.0001727056 0.9846775
      Equation cpiUSA 0.0010306707 0.4085665  0.993779851  0.0016630729 0.7017294
                         cpiUSA -1         Trend    dolcan -1   cpiUSA -1
      Equation dolcan 0.0006501501 -0.0002075003  0.970279328 0.000997695
      Equation cpiUSA 0.9898733656 -0.0013705132 -0.007164849 1.005870624
      
      [[25]]
                        Intercept         Trend   dolcan -1   cpiUSA -1  dolcan -2
      Equation dolcan  0.09231288  0.0006969576  1.11640771 0.005716144 -0.1476973
      Equation cpiUSA -0.38461556 -0.0006441612 -0.05789139 1.209808395  0.6664581
                         cpiUSA -2    Intercept         Trend dolcan -1    cpiUSA -1
      Equation dolcan -0.008053514 0.0006464748 -0.0001118549  1.141001 -0.003541846
      Equation cpiUSA -0.210683537 1.6398792639  0.0059730476  0.572203  1.334289684
                       dolcan -2   cpiUSA -2
      Equation dolcan -0.1549285  0.00404909
      Equation cpiUSA -1.1084546 -0.35828351
      
      [[26]]
                        Intercept         Trend   dolcan -1   cpiUSA -1  dolcan -2
      Equation dolcan  0.09231288  0.0006969576  1.11640771 0.005716144 -0.1476973
      Equation cpiUSA -0.38461556 -0.0006441612 -0.05789139 1.209808395  0.6664581
                         cpiUSA -2   Intercept        Trend dolcan -1    cpiUSA -1
      Equation dolcan -0.008053514 -0.01537662 -0.000157434 1.1806500 0.0009972436
      Equation cpiUSA -0.210683537  1.50018865  0.005426318 0.5086428 1.3201624875
                       dolcan -2     cpiUSA -2 Intercept       Trend dolcan -1
      Equation dolcan -0.1865112 -0.0003233887 0.3044483 0.001169578 0.9527351
      Equation cpiUSA -0.9522052 -0.3424836706 1.2815143 0.009182135 2.0104789
                        cpiUSA -1   dolcan -2   cpiUSA -2
      Equation dolcan -0.02974125 -0.07704594  0.02531077
      Equation cpiUSA  1.30743327 -2.39121154 -0.33933921
      
      [[27]]
                        Intercept dolcan -1   cpiUSA -1   dolcan -2   cpiUSA -2
      Equation dolcan  0.01395885 1.1438210 0.002560297 -0.15743623 -0.00252461
      Equation cpiUSA -0.64727526 0.8354151 1.394781489  0.01498696 -0.39677145
                       Intercept dolcan -1   cpiUSA -1  dolcan -2   cpiUSA -2
      Equation dolcan 0.04017077 1.1335259 -0.01023868 -0.1713070  0.01039812
      Equation cpiUSA 0.28300510 0.4492155  1.31096745 -0.5821929 -0.31059283
      
      [[28]]
                        Intercept dolcan -1   cpiUSA -1   dolcan -2   cpiUSA -2
      Equation dolcan  0.01395885 1.1438210 0.002560297 -0.15743623 -0.00252461
      Equation cpiUSA -0.64727526 0.8354151 1.394781489  0.01498696 -0.39677145
                        Intercept dolcan -1    cpiUSA -1  dolcan -2    cpiUSA -2
      Equation dolcan  0.02333395 1.2135735 -0.004115204 -0.2373632  0.004240793
      Equation cpiUSA -0.05032674 0.5533174  1.269024506 -0.3582838 -0.269630556
                        Intercept dolcan -1   cpiUSA -1   dolcan -2   cpiUSA -2
      Equation dolcan  0.05798003 0.9631114 -0.02791145 -0.02695616  0.02829938
      Equation cpiUSA -0.65346050 2.0919411  1.32179865 -1.99796630 -0.31587620
      
      [[29]]
                      dolcan -1     cpiUSA -1  dolcan -2     cpiUSA -2 dolcan -1
      Equation dolcan 1.1452641 -0.0005511544 -0.1435657  0.0005498422 1.1609381
      Equation cpiUSA 0.6624965  1.4943161966 -0.5053663 -0.4946045183 0.5957397
                         cpiUSA -1  dolcan -2    cpiUSA -2
      Equation dolcan -0.007105683 -0.1655764  0.007202736
      Equation cpiUSA  1.314457492 -0.4870101 -0.314595092
      
      [[30]]
                        dolcan -1   cpiUSA -1  dolcan -2    cpiUSA -2 dolcan -1
      Equation dolcan  1.19639319 0.003697079 -0.1986477 -0.003599167 1.1265025
      Equation cpiUSA -0.03686374 1.224304032  0.2957565 -0.226656499 0.9516258
                        cpiUSA -1  dolcan -2    cpiUSA -2 dolcan -1    cpiUSA -1
      Equation dolcan 0.003469589 -0.1289219 -0.003432864 1.1589176 -0.008044606
      Equation cpiUSA 1.378481433 -0.6090764 -0.380906894 0.6281012  1.326424793
                       dolcan -2    cpiUSA -2
      Equation dolcan -0.1627678  0.008132882
      Equation cpiUSA -0.5220339 -0.326551571
      
      [[31]]
                              Trend dolcan -1     cpiUSA -1  dolcan -2     cpiUSA -2
      Equation dolcan -1.622245e-05 1.1450492 -0.0004092175 -0.1449554  0.0004693694
      Equation cpiUSA  3.096566e-03 0.7035062  1.4672230481 -0.2400976 -0.4792437363
                              Trend dolcan -1    cpiUSA -1  dolcan -2   cpiUSA -2
      Equation dolcan -0.0002216198 1.1270540 -0.009180805 -0.1574311  0.01024224
      Equation cpiUSA -0.0005662404 0.5091657  1.309155537 -0.4661987 -0.30682914
      
      [[32]]
                              Trend dolcan -1   cpiUSA -1  dolcan -2    cpiUSA -2
      Equation dolcan -6.292985e-05  1.155769 0.001521385 -0.1604489 -0.001290011
      Equation cpiUSA  3.007284e-03  1.053381 1.507577719 -0.6132165 -0.519312189
                              Trend dolcan -1    cpiUSA -1   dolcan -2   cpiUSA -2
      Equation dolcan -0.0001870489 1.1844185 -0.005661468 -0.21109514  0.00657114
      Equation cpiUSA  0.0011972003 0.3427475  1.157796683 -0.04378136 -0.16328716
                              Trend dolcan -1   cpiUSA -1  dolcan -2   cpiUSA -2
      Equation dolcan -0.0001633523 0.9858677 -0.02757334 -0.0233897  0.02860914
      Equation cpiUSA  0.0035714290 2.1499439  1.31655866 -2.1653563 -0.32545538
      
      [[33]]
                              ECT       Const         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan -0.05165772 0.005288736 -3.885903e-05   0.1675227  -0.0121845
      Equation cpiUSA  0.73785216 0.623866396 -1.062156e-03  -0.1019883   0.1395116
                               ECT       Const         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan -0.001023032 0.003613178 -3.313512e-05   0.1312604 0.002326573
      Equation cpiUSA  0.015485020 0.101293958 -3.625345e-05   0.2371179 0.620902189
      
      [[34]]
                              ECT       Const         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan -0.05165772 0.005288736 -3.885903e-05   0.1675227  -0.0121845
      Equation cpiUSA  0.73785216 0.623866396 -1.062156e-03  -0.1019883   0.1395116
                               ECT      Const         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan -0.002269518 0.01440822 -0.0001021736  0.08504187 0.003521756
      Equation cpiUSA -0.337573130 0.36598485 -0.0011856393  1.18867823 0.481646459
                              ECT       Const        Trend dolcan t -1  cpiUSA t -1
      Equation dolcan -0.02900588  0.01680242 3.615931e-06   0.1850201 -0.007791395
      Equation cpiUSA  0.33504178 -0.12199902 1.994943e-03  -1.5209576  0.614156675
      
      [[35]]
                              ECT        Const dolcan t -1 cpiUSA t -1         ECT
      Equation dolcan -0.03750343 -0.006852286   0.1250896 -0.01024718 0.006057814
      Equation cpiUSA  0.70847089  0.400683554   0.1907474  0.19818232 0.058696937
                             Const dolcan t -1 cpiUSA t -1
      Equation dolcan -0.001555704   0.1644276 0.002927573
      Equation cpiUSA  0.085773731   0.2256604 0.619010979
      
      [[36]]
                              ECT        Const dolcan t -1 cpiUSA t -1         ECT
      Equation dolcan -0.03750343 -0.006852286   0.1250896 -0.01024718  0.04152116
      Equation cpiUSA  0.70847089  0.400683554   0.1907474  0.19818232 -0.35664491
                             Const dolcan t -1 cpiUSA t -1           ECT        Const
      Equation dolcan -0.006533316   0.2580703  0.01122303 -0.0007475249 0.0007799121
      Equation cpiUSA  0.184529474   0.4800999  0.27966286 -0.0967219539 0.1360188778
                      dolcan t -1 cpiUSA t -1
      Equation dolcan   0.1052364 0.002985671
      Equation cpiUSA  -0.7272839 0.680080621
      
      [[37]]
                              ECT dolcan t -1 cpiUSA t -1        ECT dolcan t -1
      Equation dolcan -0.02202667    0.117653 -0.01446972 0.00279744   0.1680402
      Equation cpiUSA -0.55070230    1.016250  0.32487254 0.18971894  -0.2859286
                      cpiUSA t -1
      Equation dolcan 0.001035012
      Equation cpiUSA 0.750443210
      
      [[38]]
                              ECT dolcan t -1 cpiUSA t -1         ECT dolcan t -1
      Equation dolcan -0.02522709   0.1674051 -0.01801245 -0.01545396  0.07333102
      Equation cpiUSA -0.45272804   2.2623818  0.24296512 -0.92219950 -0.54141665
                       cpiUSA t -1         ECT dolcan t -1  cpiUSA t -1
      Equation dolcan -0.009306737 0.003049438   0.1689193 0.0007218864
      Equation cpiUSA  0.301735719 0.194655833  -0.2687060 0.7443087565
      
      [[39]]
                             ECT         Trend dolcan t -1 cpiUSA t -1         ECT
      Equation dolcan -0.0599990 -4.729002e-05   0.2434441  -0.0116853 0.002438058
      Equation cpiUSA -0.5717297  4.019271e-06   1.4443120   0.1073026 0.188165815
                             Trend dolcan t -1  cpiUSA t -1
      Equation dolcan 9.075573e-08   0.1271129 0.0007793217
      Equation cpiUSA 5.901968e-04  -0.1180760 0.5950325178
      
      [[40]]
                              ECT         Trend dolcan t -1 cpiUSA t -1         ECT
      Equation dolcan -0.06592635 -4.893727e-05   0.1448262  -0.0165236 0.002513236
      Equation cpiUSA -0.74852793 -3.444909e-04   2.2377120   0.2191881 0.080259228
                              Trend dolcan t -1 cpiUSA t -1         ECT         Trend
      Equation dolcan -3.076313e-07   0.1907503 0.002161459 0.004302129 -0.0000163670
      Equation cpiUSA  7.917784e-04  -0.4850491 0.383556172 0.129046087  0.0007655733
                      dolcan t -1 cpiUSA t -1
      Equation dolcan  0.08396771 0.001005123
      Equation cpiUSA -0.67724374 0.707856779
      
      [[41]]
                              ECT      Const         Trend dolcan t -1  cpiUSA t -1
      Equation dolcan -0.01518849 0.01105132 -4.790127e-05   0.1838534 -0.004546173
      Equation cpiUSA -0.25113994 0.33469736 -8.430391e-04   0.5327942  0.337989478
                      dolcan t -2  cpiUSA t -2         ECT       Const         Trend
      Equation dolcan  0.02939966 -0.002047247 -0.01284187  0.00991402 -3.535269e-05
      Equation cpiUSA -1.22933741 -0.162023742  0.18220945 -0.04146424  1.147708e-03
                      dolcan t -1 cpiUSA t -1 dolcan t -2 cpiUSA t -2
      Equation dolcan    0.125544 0.005828206  -0.1979610 -0.01062449
      Equation cpiUSA   -1.067956 0.418093241   0.9412034  0.30705642
      
      [[42]]
                              ECT      Const         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan -0.01739221 0.01123809 -4.874173e-05   0.1570717 -0.01032538
      Equation cpiUSA -0.18905435 0.28407694 -5.662286e-04   0.4259146  0.28726509
                      dolcan t -2  cpiUSA t -2        ECT      Const         Trend
      Equation dolcan  0.05235467  0.001209455 -0.3473239 0.07524714 -5.319414e-05
      Equation cpiUSA -1.18744412 -0.119469224 -6.1188098 2.70914734 -9.002670e-03
                      dolcan t -1 cpiUSA t -1 dolcan t -2 cpiUSA t -2        ECT
      Equation dolcan   0.1602594  0.03364773   0.1326966 -0.03037366 -0.0119483
      Equation cpiUSA  -1.9196425  0.43896490  -4.9984845 -0.24567650  0.3550181
                             Const        Trend dolcan t -1  cpiUSA t -1 dolcan t -2
      Equation dolcan  0.009739977 2.161084e-05   0.1143009 -0.002826679  -0.2251498
      Equation cpiUSA -0.113831430 1.882903e-03  -1.3483856  0.368199325   0.4192614
                      cpiUSA t -2
      Equation dolcan  -0.0106057
      Equation cpiUSA   0.2529322
      
      [[43]]
                               ECT       Const dolcan t -1  cpiUSA t -1 dolcan t -2
      Equation dolcan -0.002677087 0.001901141   0.1821789 -0.003349176  0.02289789
      Equation cpiUSA -0.030945320 0.173658592   0.5033226  0.359056047 -1.34376545
                        cpiUSA t -2          ECT      Const dolcan t -1 cpiUSA t -1
      Equation dolcan -0.0007533507 -0.001946964 0.00412689   0.1173725 0.005439081
      Equation cpiUSA -0.1392518022 -0.171488350 0.14641214  -0.8026721 0.430726006
                      dolcan t -2 cpiUSA t -2
      Equation dolcan  -0.2105055 -0.01099955
      Equation cpiUSA   1.3484553  0.31923260
      
      [[44]]
                             ECT       Const dolcan t -1 cpiUSA t -1 dolcan t -2
      Equation dolcan -0.0461568 -0.01011165   0.1582452 -0.01119553  0.01842045
      Equation cpiUSA  0.9601284  0.51899776  -0.3722760  0.20683381 -1.60201032
                       cpiUSA t -2        ECT        Const dolcan t -1 cpiUSA t -1
      Equation dolcan  0.003495803  0.0204315 -0.003717539   0.1867325  0.01718006
      Equation cpiUSA -0.161519009 -0.1547078  0.142118177   1.4258148  0.58946004
                      dolcan t -2  cpiUSA t -2         ECT      Const dolcan t -1
      Equation dolcan  0.01977927 -0.009962532 -0.01764468 0.01267829   0.1202763
      Equation cpiUSA -1.18049183 -0.135171254 -0.14129491 0.14217690  -0.8277689
                       cpiUSA t -1 dolcan t -2 cpiUSA t -2
      Equation dolcan -0.002199864  -0.2165529 -0.01009317
      Equation cpiUSA  0.422812224   1.1682873  0.29758780
      
      [[45]]
                              ECT dolcan t -1 cpiUSA t -1  dolcan t -2  cpiUSA t -2
      Equation dolcan -0.01791467   0.1391735 -0.01328401  0.002174167  0.001436122
      Equation cpiUSA -0.56534507   0.6575878  0.31447878 -0.623941905 -0.088049070
                              ECT  dolcan t -1 cpiUSA t -1 dolcan t -2 cpiUSA t -2
      Equation dolcan 0.006491931  0.160034945 0.007915646 -0.06608634 -0.01065185
      Equation cpiUSA 0.154651435 -0.008908433 0.671157310 -0.91581701  0.13182613
      
      [[46]]
                              ECT dolcan t -1  cpiUSA t -1 dolcan t -2  cpiUSA t -2
      Equation dolcan -0.01682907   0.2487917 -0.009049603  -0.1622036 -0.000472127
      Equation cpiUSA -0.63057527   0.7048539  0.113211738   0.7159188 -0.144705191
                              ECT dolcan t -1 cpiUSA t -1 dolcan t -2 cpiUSA t -2
      Equation dolcan -0.01504519  0.03394955 -0.02606875   0.2262093  0.01536797
      Equation cpiUSA -0.70956313 -0.33245260  0.43942156  -1.4005139 -0.03286632
                              ECT dolcan t -1 cpiUSA t -1 dolcan t -2 cpiUSA t -2
      Equation dolcan 0.005762436   0.1725675  0.01049994 -0.08084231 -0.01247178
      Equation cpiUSA 0.152098278  -0.1468867  0.63963012 -1.13241543  0.16457327
      
      [[47]]
                               ECT        Trend dolcan t -1  cpiUSA t -1 dolcan t -2
      Equation dolcan -0.001594506 5.687498e-06   0.1856852 -0.001987102  0.02346587
      Equation cpiUSA  0.160563744 7.799357e-04   0.5882712  0.415492814 -1.40904656
                        cpiUSA t -2         ECT        Trend dolcan t -1 cpiUSA t -1
      Equation dolcan  0.0005909138 0.006034676 7.266438e-06   0.1190289 0.006094024
      Equation cpiUSA -0.0821251317 0.103260490 9.694584e-04  -1.0407069 0.416981490
                      dolcan t -2  cpiUSA t -2
      Equation dolcan  -0.2116016 -0.009724057
      Equation cpiUSA   0.9982535  0.303290461
      
      [[48]]
                               ECT        Trend dolcan t -1  cpiUSA t -1 dolcan t -2
      Equation dolcan -0.002995578 6.111664e-06   0.1616106 -0.007636568  0.04661039
      Equation cpiUSA  0.174864230 8.203576e-04   0.5406502  0.355232986 -1.33264807
                       cpiUSA t -2         ECT         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan  0.004104931 -0.05993876  6.599656e-05   0.1751623  0.04194837
      Equation cpiUSA -0.046277248  4.22801254 -4.711408e-03  -1.3830880  0.73781541
                      dolcan t -2 cpiUSA t -2         ECT        Trend dolcan t -1
      Equation dolcan   0.2039924 -0.03413151 0.006794794 6.522883e-05   0.1066869
      Equation cpiUSA  -2.4316002 -0.38097151 0.135966975 1.373138e-03  -1.2593996
                       cpiUSA t -1 dolcan t -2 cpiUSA t -2
      Equation dolcan -0.002756671  -0.2406763 -0.01007585
      Equation cpiUSA  0.367381138   0.6007203  0.24673984
      

---

    Code
      sapply(mods, tsDyn:::coefVec.nlVar)
    Output
      [[1]]
                                         [,1]
      Intercept Equation dolcan  1.044232e-02
      Trend Equation dolcan     -2.172806e-05
      dolcan -1 Equation dolcan  9.871079e-01
      cpiUSA -1 Equation dolcan  1.455020e-04
      Intercept Equation cpiUSA  3.789956e-01
      Trend Equation cpiUSA     -7.098272e-04
      dolcan -1 Equation cpiUSA -1.395487e-01
      cpiUSA -1 Equation cpiUSA  1.002279e+00
      
      [[2]]
                                         [,1]
      Intercept Equation dolcan  1.278017e-02
      dolcan -1 Equation dolcan  9.870607e-01
      cpiUSA -1 Equation dolcan  6.443813e-05
      Intercept Equation cpiUSA  4.553701e-01
      dolcan -1 Equation cpiUSA -1.410917e-01
      cpiUSA -1 Equation cpiUSA  9.996312e-01
      
      [[3]]
                                        [,1]
      dolcan -1 Equation dolcan 1.000411e+00
      cpiUSA -1 Equation dolcan 1.200969e-05
      dolcan -1 Equation cpiUSA 3.346056e-01
      cpiUSA -1 Equation cpiUSA 9.977632e-01
      
      [[4]]
                                         [,1]
      Trend Equation dolcan     -0.0000766771
      dolcan -1 Equation dolcan  0.9919596075
      cpiUSA -1 Equation dolcan  0.0003319243
      Trend Equation cpiUSA     -0.0027041580
      dolcan -1 Equation cpiUSA  0.0365395423
      cpiUSA -1 Equation cpiUSA  1.0090455285
      
      [[5]]
                                         [,1]
      Intercept Equation dolcan  1.360877e-02
      Trend Equation dolcan     -2.594234e-05
      dolcan -1 Equation dolcan  1.150391e+00
      cpiUSA -1 Equation dolcan -3.742394e-03
      dolcan -2 Equation dolcan -1.660825e-01
      cpiUSA -2 Equation dolcan  3.911853e-03
      Intercept Equation cpiUSA  2.832064e-01
      Trend Equation cpiUSA      3.472133e-04
      dolcan -1 Equation cpiUSA  5.331984e-01
      cpiUSA -1 Equation cpiUSA  1.457513e+00
      dolcan -2 Equation cpiUSA -6.049136e-01
      cpiUSA -2 Equation cpiUSA -4.591016e-01
      
      [[6]]
                                        [,1]
      Intercept Equation dolcan  0.016368553
      dolcan -1 Equation dolcan  1.150672447
      cpiUSA -1 Equation dolcan -0.003709029
      dolcan -2 Equation dolcan -0.166402698
      cpiUSA -2 Equation dolcan  0.003781728
      Intercept Equation cpiUSA  0.246269281
      dolcan -1 Equation cpiUSA  0.529429760
      cpiUSA -1 Equation cpiUSA  1.457066046
      dolcan -2 Equation cpiUSA -0.600628381
      cpiUSA -2 Equation cpiUSA -0.457359952
      
      [[7]]
                                        [,1]
      dolcan -1 Equation dolcan  1.160966574
      cpiUSA -1 Equation dolcan -0.001947310
      dolcan -2 Equation dolcan -0.160200153
      cpiUSA -2 Equation dolcan  0.001956871
      dolcan -1 Equation cpiUSA  0.684307680
      cpiUSA -1 Equation cpiUSA  1.483571582
      dolcan -2 Equation cpiUSA -0.507309431
      cpiUSA -2 Equation cpiUSA -0.484815406
      
      [[8]]
                                         [,1]
      Trend Equation dolcan     -9.611456e-05
      dolcan -1 Equation dolcan  1.153493e+00
      cpiUSA -1 Equation dolcan -3.171405e-03
      dolcan -2 Equation dolcan -1.628884e-01
      cpiUSA -2 Equation dolcan  3.578896e-03
      Trend Equation cpiUSA     -1.113111e-03
      dolcan -1 Equation cpiUSA  5.977555e-01
      cpiUSA -1 Equation cpiUSA  1.469395e+00
      dolcan -2 Equation cpiUSA -5.384425e-01
      cpiUSA -2 Equation cpiUSA -4.660306e-01
      
      [[9]]
                                         [,1]
      ECT Equation dolcan       -1.506974e-02
      Intercept Equation dolcan  1.098652e-02
      Trend Equation dolcan     -4.515916e-05
      dolcan -1 Equation dolcan  1.655667e-01
      cpiUSA -1 Equation dolcan -3.976831e-03
      ECT Equation cpiUSA       -5.033679e-02
      Intercept Equation cpiUSA  1.930659e-01
      Trend Equation cpiUSA     -3.133701e-04
      dolcan -1 Equation cpiUSA  5.871833e-01
      cpiUSA -1 Equation cpiUSA  4.568679e-01
      
      [[10]]
                                        [,1]
      ECT Equation dolcan       -0.001106223
      Intercept Equation dolcan  0.002123837
      dolcan -1 Equation dolcan  0.159535446
      cpiUSA -1 Equation dolcan -0.003007480
      ECT Equation cpiUSA        0.046559351
      Intercept Equation cpiUSA  0.131565673
      dolcan -1 Equation cpiUSA  0.545330785
      cpiUSA -1 Equation cpiUSA  0.463594476
      
      [[11]]
                                         [,1]
      ECT Equation dolcan       -0.0006768995
      dolcan -1 Equation dolcan  0.1681195532
      cpiUSA -1 Equation dolcan  0.0024290305
      ECT Equation cpiUSA        0.0731547016
      dolcan -1 Equation cpiUSA  1.0770918513
      cpiUSA -1 Equation cpiUSA  0.8003708588
      
      [[12]]
                                         [,1]
      ECT Equation dolcan        1.372394e-03
      Trend Equation dolcan      7.603074e-06
      dolcan -1 Equation dolcan  1.610732e-01
      cpiUSA -1 Equation dolcan -1.227277e-03
      ECT Equation cpiUSA        2.386006e-01
      Trend Equation cpiUSA      6.138198e-04
      dolcan -1 Equation cpiUSA  5.082187e-01
      cpiUSA -1 Equation cpiUSA  5.051858e-01
      
      [[13]]
                                         [,1]
      ECT Equation dolcan       -1.524527e-02
      Intercept Equation dolcan  1.156071e-02
      Trend Equation dolcan     -4.607701e-05
      dolcan -1 Equation dolcan  1.678529e-01
      cpiUSA -1 Equation dolcan -2.542136e-03
      dolcan -2 Equation dolcan -2.196011e-02
      cpiUSA -2 Equation dolcan -3.088730e-03
      ECT Equation cpiUSA       -3.983784e-02
      Intercept Equation cpiUSA  1.848548e-01
      Trend Equation cpiUSA     -2.787031e-04
      dolcan -1 Equation cpiUSA  6.918440e-01
      cpiUSA -1 Equation cpiUSA  4.538527e-01
      dolcan -2 Equation cpiUSA -6.301134e-01
      cpiUSA -2 Equation cpiUSA  1.100919e-02
      
      [[14]]
                                        [,1]
      ECT Equation dolcan       -0.001031649
      Intercept Equation dolcan  0.002457958
      dolcan -1 Equation dolcan  0.163512671
      cpiUSA -1 Equation dolcan -0.001892005
      dolcan -2 Equation dolcan -0.030335891
      cpiUSA -2 Equation dolcan -0.002298449
      ECT Equation cpiUSA        0.046135183
      Intercept Equation cpiUSA  0.129795583
      dolcan -1 Equation cpiUSA  0.665591703
      cpiUSA -1 Equation cpiUSA  0.457785100
      dolcan -2 Equation cpiUSA -0.680775403
      cpiUSA -2 Equation cpiUSA  0.015789314
      
      [[15]]
                                         [,1]
      ECT Equation dolcan       -0.0007525866
      dolcan -1 Equation dolcan  0.1722870283
      cpiUSA -1 Equation dolcan  0.0015294295
      dolcan -2 Equation dolcan -0.0251824688
      cpiUSA -2 Equation dolcan  0.0012221382
      ECT Equation cpiUSA        0.0608714297
      dolcan -1 Equation cpiUSA  1.1289328140
      cpiUSA -1 Equation cpiUSA  0.6384583142
      dolcan -2 Equation cpiUSA -0.4086424096
      cpiUSA -2 Equation cpiUSA  0.2016984025
      
      [[16]]
                                         [,1]
      ECT Equation dolcan        1.657334e-03
      Trend Equation dolcan      8.424941e-06
      dolcan -1 Equation dolcan  1.655519e-01
      cpiUSA -1 Equation dolcan -6.685009e-04
      dolcan -2 Equation dolcan -3.020354e-02
      cpiUSA -2 Equation dolcan -1.017307e-03
      ECT Equation cpiUSA        2.304334e-01
      Trend Equation cpiUSA      5.927787e-04
      dolcan -1 Equation cpiUSA  6.550524e-01
      cpiUSA -1 Equation cpiUSA  4.838120e-01
      dolcan -2 Equation cpiUSA -7.619252e-01
      cpiUSA -2 Equation cpiUSA  4.413108e-02
      
      [[17]]
                                         [,1]
      Intercept Equation dolcan  0.0821016594
      Trend Equation dolcan      0.0007054584
      dolcan -1 Equation dolcan  0.9822978662
      cpiUSA -1 Equation dolcan -0.0024086192
      Intercept Equation dolcan -0.0135404152
      Trend Equation dolcan     -0.0001653546
      dolcan -1 Equation dolcan  0.9914561139
      cpiUSA -1 Equation dolcan  0.0007168335
      Intercept Equation cpiUSA -0.5558800278
      Trend Equation cpiUSA     -0.0017977750
      dolcan -1 Equation cpiUSA  0.7181117320
      cpiUSA -1 Equation cpiUSA  1.0032098335
      Intercept Equation cpiUSA  2.2020698586
      Trend Equation cpiUSA      0.0073398658
      dolcan -1 Equation cpiUSA -0.7468369449
      cpiUSA -1 Equation cpiUSA  0.9702240258
      
      [[18]]
                                         [,1]
      Intercept Equation dolcan  0.0821016594
      Trend Equation dolcan      0.0007054584
      dolcan -1 Equation dolcan  0.9822978662
      cpiUSA -1 Equation dolcan -0.0024086192
      Intercept Equation dolcan -0.0212780040
      Trend Equation dolcan     -0.0002016937
      dolcan -1 Equation dolcan  0.9952333026
      cpiUSA -1 Equation dolcan  0.0008405793
      Intercept Equation dolcan  0.3579121229
      Trend Equation dolcan      0.0014262824
      dolcan -1 Equation dolcan  0.8674328244
      cpiUSA -1 Equation dolcan -0.0056052408
      Intercept Equation cpiUSA -0.5558800278
      Trend Equation cpiUSA     -0.0017977750
      dolcan -1 Equation cpiUSA  0.7181117320
      cpiUSA -1 Equation cpiUSA  1.0032098335
      Intercept Equation cpiUSA  2.0412562528
      Trend Equation cpiUSA      0.0064145304
      dolcan -1 Equation cpiUSA -0.6685873182
      cpiUSA -1 Equation cpiUSA  0.9732236367
      Intercept Equation cpiUSA -1.2735570315
      Trend Equation cpiUSA     -0.0035048123
      dolcan -1 Equation cpiUSA  0.1599708690
      cpiUSA -1 Equation cpiUSA  1.0214915833
      
      [[19]]
                                         [,1]
      Intercept Equation dolcan  7.691528e-03
      dolcan -1 Equation dolcan  9.939684e-01
      cpiUSA -1 Equation dolcan  1.357332e-05
      Intercept Equation dolcan  3.045117e-02
      dolcan -1 Equation dolcan  9.689214e-01
      cpiUSA -1 Equation dolcan  1.456640e-04
      Intercept Equation cpiUSA -1.012510e+00
      dolcan -1 Equation cpiUSA  1.341579e+00
      cpiUSA -1 Equation cpiUSA  9.969201e-01
      Intercept Equation cpiUSA  4.650176e-01
      dolcan -1 Equation cpiUSA -2.412925e-01
      cpiUSA -1 Equation cpiUSA  1.000712e+00
      
      [[20]]
                                         [,1]
      Intercept Equation dolcan -2.978663e-02
      dolcan -1 Equation dolcan  1.031097e+00
      cpiUSA -1 Equation dolcan  8.444228e-06
      Intercept Equation dolcan  6.718762e-02
      dolcan -1 Equation dolcan  9.433910e-01
      cpiUSA -1 Equation dolcan  1.089975e-05
      Intercept Equation dolcan  3.045117e-02
      dolcan -1 Equation dolcan  9.689214e-01
      cpiUSA -1 Equation dolcan  1.456640e-04
      Intercept Equation cpiUSA -2.707463e-01
      dolcan -1 Equation cpiUSA  5.937531e-01
      cpiUSA -1 Equation cpiUSA  9.970502e-01
      Intercept Equation cpiUSA  1.918453e+00
      dolcan -1 Equation cpiUSA -1.055018e+00
      cpiUSA -1 Equation cpiUSA  9.957302e-01
      Intercept Equation cpiUSA  4.650176e-01
      dolcan -1 Equation cpiUSA -2.412925e-01
      cpiUSA -1 Equation cpiUSA  1.000712e+00
      
      [[21]]
                                         [,1]
      dolcan -1 Equation dolcan  1.001838e+00
      cpiUSA -1 Equation dolcan -5.982604e-06
      dolcan -1 Equation dolcan  9.949031e-01
      cpiUSA -1 Equation dolcan  9.159023e-05
      dolcan -1 Equation cpiUSA  3.056463e-01
      cpiUSA -1 Equation cpiUSA  9.994945e-01
      dolcan -1 Equation cpiUSA  1.554722e-01
      cpiUSA -1 Equation cpiUSA  9.998859e-01
      
      [[22]]
                                        [,1]
      dolcan -1 Equation dolcan 9.994742e-01
      cpiUSA -1 Equation dolcan 8.089218e-05
      dolcan -1 Equation dolcan 9.998115e-01
      cpiUSA -1 Equation dolcan 1.872682e-05
      dolcan -1 Equation dolcan 9.949031e-01
      cpiUSA -1 Equation dolcan 9.159023e-05
      dolcan -1 Equation cpiUSA 3.063167e-01
      cpiUSA -1 Equation cpiUSA 9.977088e-01
      dolcan -1 Equation cpiUSA 5.559962e-01
      cpiUSA -1 Equation cpiUSA 9.959537e-01
      dolcan -1 Equation cpiUSA 1.554722e-01
      cpiUSA -1 Equation cpiUSA 9.998859e-01
      
      [[23]]
                                         [,1]
      Trend Equation dolcan      2.085515e-05
      dolcan -1 Equation dolcan  1.003846e+00
      cpiUSA -1 Equation dolcan -8.488800e-05
      Trend Equation dolcan     -2.075003e-04
      dolcan -1 Equation dolcan  9.702793e-01
      cpiUSA -1 Equation dolcan  9.976950e-04
      Trend Equation cpiUSA      3.944688e-03
      dolcan -1 Equation cpiUSA  6.854827e-01
      cpiUSA -1 Equation cpiUSA  9.845697e-01
      Trend Equation cpiUSA     -1.370513e-03
      dolcan -1 Equation cpiUSA -7.164849e-03
      cpiUSA -1 Equation cpiUSA  1.005871e+00
      
      [[24]]
                                         [,1]
      Trend Equation dolcan      0.0002877063
      dolcan -1 Equation dolcan  1.0280166792
      cpiUSA -1 Equation dolcan -0.0010158416
      Trend Equation dolcan     -0.0001727056
      dolcan -1 Equation dolcan  0.9846775202
      cpiUSA -1 Equation dolcan  0.0006501501
      Trend Equation dolcan     -0.0002075003
      dolcan -1 Equation dolcan  0.9702793281
      cpiUSA -1 Equation dolcan  0.0009976950
      Trend Equation cpiUSA      0.0010306707
      dolcan -1 Equation cpiUSA  0.4085665265
      cpiUSA -1 Equation cpiUSA  0.9937798507
      Trend Equation cpiUSA      0.0016630729
      dolcan -1 Equation cpiUSA  0.7017293518
      cpiUSA -1 Equation cpiUSA  0.9898733656
      Trend Equation cpiUSA     -0.0013705132
      dolcan -1 Equation cpiUSA -0.0071648488
      cpiUSA -1 Equation cpiUSA  1.0058706240
      
      [[25]]
                                         [,1]
      Intercept Equation dolcan  0.0923128826
      Trend Equation dolcan      0.0006969576
      dolcan -1 Equation dolcan  1.1164077100
      cpiUSA -1 Equation dolcan  0.0057161444
      dolcan -2 Equation dolcan -0.1476973419
      cpiUSA -2 Equation dolcan -0.0080535138
      Intercept Equation dolcan  0.0006464748
      Trend Equation dolcan     -0.0001118549
      dolcan -1 Equation dolcan  1.1410011259
      cpiUSA -1 Equation dolcan -0.0035418464
      dolcan -2 Equation dolcan -0.1549284601
      cpiUSA -2 Equation dolcan  0.0040490903
      Intercept Equation cpiUSA -0.3846155562
      Trend Equation cpiUSA     -0.0006441612
      dolcan -1 Equation cpiUSA -0.0578913909
      cpiUSA -1 Equation cpiUSA  1.2098083954
      dolcan -2 Equation cpiUSA  0.6664581104
      cpiUSA -2 Equation cpiUSA -0.2106835368
      Intercept Equation cpiUSA  1.6398792639
      Trend Equation cpiUSA      0.0059730476
      dolcan -1 Equation cpiUSA  0.5722029695
      cpiUSA -1 Equation cpiUSA  1.3342896838
      dolcan -2 Equation cpiUSA -1.1084545888
      cpiUSA -2 Equation cpiUSA -0.3582835140
      
      [[26]]
                                         [,1]
      Intercept Equation dolcan  0.0923128826
      Trend Equation dolcan      0.0006969576
      dolcan -1 Equation dolcan  1.1164077100
      cpiUSA -1 Equation dolcan  0.0057161444
      dolcan -2 Equation dolcan -0.1476973419
      cpiUSA -2 Equation dolcan -0.0080535138
      Intercept Equation dolcan -0.0153766233
      Trend Equation dolcan     -0.0001574340
      dolcan -1 Equation dolcan  1.1806499899
      cpiUSA -1 Equation dolcan  0.0009972436
      dolcan -2 Equation dolcan -0.1865111657
      cpiUSA -2 Equation dolcan -0.0003233887
      Intercept Equation dolcan  0.3044482547
      Trend Equation dolcan      0.0011695783
      dolcan -1 Equation dolcan  0.9527351077
      cpiUSA -1 Equation dolcan -0.0297412520
      dolcan -2 Equation dolcan -0.0770459401
      cpiUSA -2 Equation dolcan  0.0253107720
      Intercept Equation cpiUSA -0.3846155562
      Trend Equation cpiUSA     -0.0006441612
      dolcan -1 Equation cpiUSA -0.0578913909
      cpiUSA -1 Equation cpiUSA  1.2098083954
      dolcan -2 Equation cpiUSA  0.6664581104
      cpiUSA -2 Equation cpiUSA -0.2106835368
      Intercept Equation cpiUSA  1.5001886521
      Trend Equation cpiUSA      0.0054263183
      dolcan -1 Equation cpiUSA  0.5086428443
      cpiUSA -1 Equation cpiUSA  1.3201624875
      dolcan -2 Equation cpiUSA -0.9522051631
      cpiUSA -2 Equation cpiUSA -0.3424836706
      Intercept Equation cpiUSA  1.2815142982
      Trend Equation cpiUSA      0.0091821349
      dolcan -1 Equation cpiUSA  2.0104788680
      cpiUSA -1 Equation cpiUSA  1.3074332679
      dolcan -2 Equation cpiUSA -2.3912115414
      cpiUSA -2 Equation cpiUSA -0.3393392148
      
      [[27]]
                                        [,1]
      Intercept Equation dolcan  0.013958849
      dolcan -1 Equation dolcan  1.143821006
      cpiUSA -1 Equation dolcan  0.002560297
      dolcan -2 Equation dolcan -0.157436226
      cpiUSA -2 Equation dolcan -0.002524610
      Intercept Equation dolcan  0.040170775
      dolcan -1 Equation dolcan  1.133525867
      cpiUSA -1 Equation dolcan -0.010238678
      dolcan -2 Equation dolcan -0.171307037
      cpiUSA -2 Equation dolcan  0.010398124
      Intercept Equation cpiUSA -0.647275259
      dolcan -1 Equation cpiUSA  0.835415106
      cpiUSA -1 Equation cpiUSA  1.394781489
      dolcan -2 Equation cpiUSA  0.014986956
      cpiUSA -2 Equation cpiUSA -0.396771445
      Intercept Equation cpiUSA  0.283005105
      dolcan -1 Equation cpiUSA  0.449215476
      cpiUSA -1 Equation cpiUSA  1.310967448
      dolcan -2 Equation cpiUSA -0.582192922
      cpiUSA -2 Equation cpiUSA -0.310592829
      
      [[28]]
                                        [,1]
      Intercept Equation dolcan  0.013958849
      dolcan -1 Equation dolcan  1.143821006
      cpiUSA -1 Equation dolcan  0.002560297
      dolcan -2 Equation dolcan -0.157436226
      cpiUSA -2 Equation dolcan -0.002524610
      Intercept Equation dolcan  0.023333951
      dolcan -1 Equation dolcan  1.213573494
      cpiUSA -1 Equation dolcan -0.004115204
      dolcan -2 Equation dolcan -0.237363161
      cpiUSA -2 Equation dolcan  0.004240793
      Intercept Equation dolcan  0.057980029
      dolcan -1 Equation dolcan  0.963111394
      cpiUSA -1 Equation dolcan -0.027911455
      dolcan -2 Equation dolcan -0.026956162
      cpiUSA -2 Equation dolcan  0.028299384
      Intercept Equation cpiUSA -0.647275259
      dolcan -1 Equation cpiUSA  0.835415106
      cpiUSA -1 Equation cpiUSA  1.394781489
      dolcan -2 Equation cpiUSA  0.014986956
      cpiUSA -2 Equation cpiUSA -0.396771445
      Intercept Equation cpiUSA -0.050326740
      dolcan -1 Equation cpiUSA  0.553317420
      cpiUSA -1 Equation cpiUSA  1.269024506
      dolcan -2 Equation cpiUSA -0.358283786
      cpiUSA -2 Equation cpiUSA -0.269630556
      Intercept Equation cpiUSA -0.653460500
      dolcan -1 Equation cpiUSA  2.091941103
      cpiUSA -1 Equation cpiUSA  1.321798654
      dolcan -2 Equation cpiUSA -1.997966295
      cpiUSA -2 Equation cpiUSA -0.315876197
      
      [[29]]
                                         [,1]
      dolcan -1 Equation dolcan  1.1452640784
      cpiUSA -1 Equation dolcan -0.0005511544
      dolcan -2 Equation dolcan -0.1435656825
      cpiUSA -2 Equation dolcan  0.0005498422
      dolcan -1 Equation dolcan  1.1609380813
      cpiUSA -1 Equation dolcan -0.0071056830
      dolcan -2 Equation dolcan -0.1655764214
      cpiUSA -2 Equation dolcan  0.0072027360
      dolcan -1 Equation cpiUSA  0.6624965086
      cpiUSA -1 Equation cpiUSA  1.4943161966
      dolcan -2 Equation cpiUSA -0.5053663429
      cpiUSA -2 Equation cpiUSA -0.4946045183
      dolcan -1 Equation cpiUSA  0.5957396608
      cpiUSA -1 Equation cpiUSA  1.3144574925
      dolcan -2 Equation cpiUSA -0.4870101219
      cpiUSA -2 Equation cpiUSA -0.3145950922
      
      [[30]]
                                        [,1]
      dolcan -1 Equation dolcan  1.196393187
      cpiUSA -1 Equation dolcan  0.003697079
      dolcan -2 Equation dolcan -0.198647695
      cpiUSA -2 Equation dolcan -0.003599167
      dolcan -1 Equation dolcan  1.126502489
      cpiUSA -1 Equation dolcan  0.003469589
      dolcan -2 Equation dolcan -0.128921882
      cpiUSA -2 Equation dolcan -0.003432864
      dolcan -1 Equation dolcan  1.158917552
      cpiUSA -1 Equation dolcan -0.008044606
      dolcan -2 Equation dolcan -0.162767842
      cpiUSA -2 Equation dolcan  0.008132882
      dolcan -1 Equation cpiUSA -0.036863739
      cpiUSA -1 Equation cpiUSA  1.224304032
      dolcan -2 Equation cpiUSA  0.295756499
      cpiUSA -2 Equation cpiUSA -0.226656499
      dolcan -1 Equation cpiUSA  0.951625838
      cpiUSA -1 Equation cpiUSA  1.378481433
      dolcan -2 Equation cpiUSA -0.609076406
      cpiUSA -2 Equation cpiUSA -0.380906894
      dolcan -1 Equation cpiUSA  0.628101159
      cpiUSA -1 Equation cpiUSA  1.326424793
      dolcan -2 Equation cpiUSA -0.522033868
      cpiUSA -2 Equation cpiUSA -0.326551571
      
      [[31]]
                                         [,1]
      Trend Equation dolcan     -1.622245e-05
      dolcan -1 Equation dolcan  1.145049e+00
      cpiUSA -1 Equation dolcan -4.092175e-04
      dolcan -2 Equation dolcan -1.449554e-01
      cpiUSA -2 Equation dolcan  4.693694e-04
      Trend Equation dolcan     -2.216198e-04
      dolcan -1 Equation dolcan  1.127054e+00
      cpiUSA -1 Equation dolcan -9.180805e-03
      dolcan -2 Equation dolcan -1.574311e-01
      cpiUSA -2 Equation dolcan  1.024224e-02
      Trend Equation cpiUSA      3.096566e-03
      dolcan -1 Equation cpiUSA  7.035062e-01
      cpiUSA -1 Equation cpiUSA  1.467223e+00
      dolcan -2 Equation cpiUSA -2.400976e-01
      cpiUSA -2 Equation cpiUSA -4.792437e-01
      Trend Equation cpiUSA     -5.662404e-04
      dolcan -1 Equation cpiUSA  5.091657e-01
      cpiUSA -1 Equation cpiUSA  1.309156e+00
      dolcan -2 Equation cpiUSA -4.661987e-01
      cpiUSA -2 Equation cpiUSA -3.068291e-01
      
      [[32]]
                                         [,1]
      Trend Equation dolcan     -6.292985e-05
      dolcan -1 Equation dolcan  1.155769e+00
      cpiUSA -1 Equation dolcan  1.521385e-03
      dolcan -2 Equation dolcan -1.604489e-01
      cpiUSA -2 Equation dolcan -1.290011e-03
      Trend Equation dolcan     -1.870489e-04
      dolcan -1 Equation dolcan  1.184419e+00
      cpiUSA -1 Equation dolcan -5.661468e-03
      dolcan -2 Equation dolcan -2.110951e-01
      cpiUSA -2 Equation dolcan  6.571140e-03
      Trend Equation dolcan     -1.633523e-04
      dolcan -1 Equation dolcan  9.858677e-01
      cpiUSA -1 Equation dolcan -2.757334e-02
      dolcan -2 Equation dolcan -2.338970e-02
      cpiUSA -2 Equation dolcan  2.860914e-02
      Trend Equation cpiUSA      3.007284e-03
      dolcan -1 Equation cpiUSA  1.053381e+00
      cpiUSA -1 Equation cpiUSA  1.507578e+00
      dolcan -2 Equation cpiUSA -6.132165e-01
      cpiUSA -2 Equation cpiUSA -5.193122e-01
      Trend Equation cpiUSA      1.197200e-03
      dolcan -1 Equation cpiUSA  3.427475e-01
      cpiUSA -1 Equation cpiUSA  1.157797e+00
      dolcan -2 Equation cpiUSA -4.378136e-02
      cpiUSA -2 Equation cpiUSA -1.632872e-01
      Trend Equation cpiUSA      3.571429e-03
      dolcan -1 Equation cpiUSA  2.149944e+00
      cpiUSA -1 Equation cpiUSA  1.316559e+00
      dolcan -2 Equation cpiUSA -2.165356e+00
      cpiUSA -2 Equation cpiUSA -3.254554e-01
      
      [[33]]
                                           [,1]
      ECT Equation dolcan         -5.165772e-02
      Const Equation dolcan        5.288736e-03
      Trend Equation dolcan       -3.885903e-05
      dolcan t -1 Equation dolcan  1.675227e-01
      cpiUSA t -1 Equation dolcan -1.218450e-02
      ECT Equation dolcan         -1.023032e-03
      Const Equation dolcan        3.613178e-03
      Trend Equation dolcan       -3.313512e-05
      dolcan t -1 Equation dolcan  1.312604e-01
      cpiUSA t -1 Equation dolcan  2.326573e-03
      ECT Equation cpiUSA          7.378522e-01
      Const Equation cpiUSA        6.238664e-01
      Trend Equation cpiUSA       -1.062156e-03
      dolcan t -1 Equation cpiUSA -1.019883e-01
      cpiUSA t -1 Equation cpiUSA  1.395116e-01
      ECT Equation cpiUSA          1.548502e-02
      Const Equation cpiUSA        1.012940e-01
      Trend Equation cpiUSA       -3.625345e-05
      dolcan t -1 Equation cpiUSA  2.371179e-01
      cpiUSA t -1 Equation cpiUSA  6.209022e-01
      
      [[34]]
                                           [,1]
      ECT Equation dolcan         -5.165772e-02
      Const Equation dolcan        5.288736e-03
      Trend Equation dolcan       -3.885903e-05
      dolcan t -1 Equation dolcan  1.675227e-01
      cpiUSA t -1 Equation dolcan -1.218450e-02
      ECT Equation dolcan         -2.269518e-03
      Const Equation dolcan        1.440822e-02
      Trend Equation dolcan       -1.021736e-04
      dolcan t -1 Equation dolcan  8.504187e-02
      cpiUSA t -1 Equation dolcan  3.521756e-03
      ECT Equation dolcan         -2.900588e-02
      Const Equation dolcan        1.680242e-02
      Trend Equation dolcan        3.615931e-06
      dolcan t -1 Equation dolcan  1.850201e-01
      cpiUSA t -1 Equation dolcan -7.791395e-03
      ECT Equation cpiUSA          7.378522e-01
      Const Equation cpiUSA        6.238664e-01
      Trend Equation cpiUSA       -1.062156e-03
      dolcan t -1 Equation cpiUSA -1.019883e-01
      cpiUSA t -1 Equation cpiUSA  1.395116e-01
      ECT Equation cpiUSA         -3.375731e-01
      Const Equation cpiUSA        3.659849e-01
      Trend Equation cpiUSA       -1.185639e-03
      dolcan t -1 Equation cpiUSA  1.188678e+00
      cpiUSA t -1 Equation cpiUSA  4.816465e-01
      ECT Equation cpiUSA          3.350418e-01
      Const Equation cpiUSA       -1.219990e-01
      Trend Equation cpiUSA        1.994943e-03
      dolcan t -1 Equation cpiUSA -1.520958e+00
      cpiUSA t -1 Equation cpiUSA  6.141567e-01
      
      [[35]]
                                          [,1]
      ECT Equation dolcan         -0.037503431
      Const Equation dolcan       -0.006852286
      dolcan t -1 Equation dolcan  0.125089575
      cpiUSA t -1 Equation dolcan -0.010247183
      ECT Equation dolcan          0.006057814
      Const Equation dolcan       -0.001555704
      dolcan t -1 Equation dolcan  0.164427612
      cpiUSA t -1 Equation dolcan  0.002927573
      ECT Equation cpiUSA          0.708470889
      Const Equation cpiUSA        0.400683554
      dolcan t -1 Equation cpiUSA  0.190747386
      cpiUSA t -1 Equation cpiUSA  0.198182315
      ECT Equation cpiUSA          0.058696937
      Const Equation cpiUSA        0.085773731
      dolcan t -1 Equation cpiUSA  0.225660438
      cpiUSA t -1 Equation cpiUSA  0.619010979
      
      [[36]]
                                           [,1]
      ECT Equation dolcan         -0.0375034305
      Const Equation dolcan       -0.0068522856
      dolcan t -1 Equation dolcan  0.1250895753
      cpiUSA t -1 Equation dolcan -0.0102471831
      ECT Equation dolcan          0.0415211583
      Const Equation dolcan       -0.0065333156
      dolcan t -1 Equation dolcan  0.2580702693
      cpiUSA t -1 Equation dolcan  0.0112230307
      ECT Equation dolcan         -0.0007475249
      Const Equation dolcan        0.0007799121
      dolcan t -1 Equation dolcan  0.1052364407
      cpiUSA t -1 Equation dolcan  0.0029856709
      ECT Equation cpiUSA          0.7084708892
      Const Equation cpiUSA        0.4006835544
      dolcan t -1 Equation cpiUSA  0.1907473865
      cpiUSA t -1 Equation cpiUSA  0.1981823150
      ECT Equation cpiUSA         -0.3566449083
      Const Equation cpiUSA        0.1845294740
      dolcan t -1 Equation cpiUSA  0.4800998537
      cpiUSA t -1 Equation cpiUSA  0.2796628566
      ECT Equation cpiUSA         -0.0967219539
      Const Equation cpiUSA        0.1360188778
      dolcan t -1 Equation cpiUSA -0.7272838694
      cpiUSA t -1 Equation cpiUSA  0.6800806207
      
      [[37]]
                                          [,1]
      ECT Equation dolcan         -0.022026668
      dolcan t -1 Equation dolcan  0.117652983
      cpiUSA t -1 Equation dolcan -0.014469722
      ECT Equation dolcan          0.002797440
      dolcan t -1 Equation dolcan  0.168040225
      cpiUSA t -1 Equation dolcan  0.001035012
      ECT Equation cpiUSA         -0.550702304
      dolcan t -1 Equation cpiUSA  1.016249503
      cpiUSA t -1 Equation cpiUSA  0.324872542
      ECT Equation cpiUSA          0.189718939
      dolcan t -1 Equation cpiUSA -0.285928640
      cpiUSA t -1 Equation cpiUSA  0.750443210
      
      [[38]]
                                           [,1]
      ECT Equation dolcan         -0.0252270875
      dolcan t -1 Equation dolcan  0.1674051146
      cpiUSA t -1 Equation dolcan -0.0180124459
      ECT Equation dolcan         -0.0154539600
      dolcan t -1 Equation dolcan  0.0733310150
      cpiUSA t -1 Equation dolcan -0.0093067368
      ECT Equation dolcan          0.0030494378
      dolcan t -1 Equation dolcan  0.1689193339
      cpiUSA t -1 Equation dolcan  0.0007218864
      ECT Equation cpiUSA         -0.4527280357
      dolcan t -1 Equation cpiUSA  2.2623817652
      cpiUSA t -1 Equation cpiUSA  0.2429651158
      ECT Equation cpiUSA         -0.9221995022
      dolcan t -1 Equation cpiUSA -0.5414166520
      cpiUSA t -1 Equation cpiUSA  0.3017357194
      ECT Equation cpiUSA          0.1946558329
      dolcan t -1 Equation cpiUSA -0.2687059931
      cpiUSA t -1 Equation cpiUSA  0.7443087565
      
      [[39]]
                                           [,1]
      ECT Equation dolcan         -5.999900e-02
      Trend Equation dolcan       -4.729002e-05
      dolcan t -1 Equation dolcan  2.434441e-01
      cpiUSA t -1 Equation dolcan -1.168530e-02
      ECT Equation dolcan          2.438058e-03
      Trend Equation dolcan        9.075573e-08
      dolcan t -1 Equation dolcan  1.271129e-01
      cpiUSA t -1 Equation dolcan  7.793217e-04
      ECT Equation cpiUSA         -5.717297e-01
      Trend Equation cpiUSA        4.019271e-06
      dolcan t -1 Equation cpiUSA  1.444312e+00
      cpiUSA t -1 Equation cpiUSA  1.073026e-01
      ECT Equation cpiUSA          1.881658e-01
      Trend Equation cpiUSA        5.901968e-04
      dolcan t -1 Equation cpiUSA -1.180760e-01
      cpiUSA t -1 Equation cpiUSA  5.950325e-01
      
      [[40]]
                                           [,1]
      ECT Equation dolcan         -6.592635e-02
      Trend Equation dolcan       -4.893727e-05
      dolcan t -1 Equation dolcan  1.448262e-01
      cpiUSA t -1 Equation dolcan -1.652360e-02
      ECT Equation dolcan          2.513236e-03
      Trend Equation dolcan       -3.076313e-07
      dolcan t -1 Equation dolcan  1.907503e-01
      cpiUSA t -1 Equation dolcan  2.161459e-03
      ECT Equation dolcan          4.302129e-03
      Trend Equation dolcan       -1.636700e-05
      dolcan t -1 Equation dolcan  8.396771e-02
      cpiUSA t -1 Equation dolcan  1.005123e-03
      ECT Equation cpiUSA         -7.485279e-01
      Trend Equation cpiUSA       -3.444909e-04
      dolcan t -1 Equation cpiUSA  2.237712e+00
      cpiUSA t -1 Equation cpiUSA  2.191881e-01
      ECT Equation cpiUSA          8.025923e-02
      Trend Equation cpiUSA        7.917784e-04
      dolcan t -1 Equation cpiUSA -4.850491e-01
      cpiUSA t -1 Equation cpiUSA  3.835562e-01
      ECT Equation cpiUSA          1.290461e-01
      Trend Equation cpiUSA        7.655733e-04
      dolcan t -1 Equation cpiUSA -6.772437e-01
      cpiUSA t -1 Equation cpiUSA  7.078568e-01
      
      [[41]]
                                           [,1]
      ECT Equation dolcan         -1.518849e-02
      Const Equation dolcan        1.105132e-02
      Trend Equation dolcan       -4.790127e-05
      dolcan t -1 Equation dolcan  1.838534e-01
      cpiUSA t -1 Equation dolcan -4.546173e-03
      dolcan t -2 Equation dolcan  2.939966e-02
      cpiUSA t -2 Equation dolcan -2.047247e-03
      ECT Equation dolcan         -1.284187e-02
      Const Equation dolcan        9.914020e-03
      Trend Equation dolcan       -3.535269e-05
      dolcan t -1 Equation dolcan  1.255440e-01
      cpiUSA t -1 Equation dolcan  5.828206e-03
      dolcan t -2 Equation dolcan -1.979610e-01
      cpiUSA t -2 Equation dolcan -1.062449e-02
      ECT Equation cpiUSA         -2.511399e-01
      Const Equation cpiUSA        3.346974e-01
      Trend Equation cpiUSA       -8.430391e-04
      dolcan t -1 Equation cpiUSA  5.327942e-01
      cpiUSA t -1 Equation cpiUSA  3.379895e-01
      dolcan t -2 Equation cpiUSA -1.229337e+00
      cpiUSA t -2 Equation cpiUSA -1.620237e-01
      ECT Equation cpiUSA          1.822094e-01
      Const Equation cpiUSA       -4.146424e-02
      Trend Equation cpiUSA        1.147708e-03
      dolcan t -1 Equation cpiUSA -1.067956e+00
      cpiUSA t -1 Equation cpiUSA  4.180932e-01
      dolcan t -2 Equation cpiUSA  9.412034e-01
      cpiUSA t -2 Equation cpiUSA  3.070564e-01
      
      [[42]]
                                           [,1]
      ECT Equation dolcan         -1.739221e-02
      Const Equation dolcan        1.123809e-02
      Trend Equation dolcan       -4.874173e-05
      dolcan t -1 Equation dolcan  1.570717e-01
      cpiUSA t -1 Equation dolcan -1.032538e-02
      dolcan t -2 Equation dolcan  5.235467e-02
      cpiUSA t -2 Equation dolcan  1.209455e-03
      ECT Equation dolcan         -3.473239e-01
      Const Equation dolcan        7.524714e-02
      Trend Equation dolcan       -5.319414e-05
      dolcan t -1 Equation dolcan  1.602594e-01
      cpiUSA t -1 Equation dolcan  3.364773e-02
      dolcan t -2 Equation dolcan  1.326966e-01
      cpiUSA t -2 Equation dolcan -3.037366e-02
      ECT Equation dolcan         -1.194830e-02
      Const Equation dolcan        9.739977e-03
      Trend Equation dolcan        2.161084e-05
      dolcan t -1 Equation dolcan  1.143009e-01
      cpiUSA t -1 Equation dolcan -2.826679e-03
      dolcan t -2 Equation dolcan -2.251498e-01
      cpiUSA t -2 Equation dolcan -1.060570e-02
      ECT Equation cpiUSA         -1.890543e-01
      Const Equation cpiUSA        2.840769e-01
      Trend Equation cpiUSA       -5.662286e-04
      dolcan t -1 Equation cpiUSA  4.259146e-01
      cpiUSA t -1 Equation cpiUSA  2.872651e-01
      dolcan t -2 Equation cpiUSA -1.187444e+00
      cpiUSA t -2 Equation cpiUSA -1.194692e-01
      ECT Equation cpiUSA         -6.118810e+00
      Const Equation cpiUSA        2.709147e+00
      Trend Equation cpiUSA       -9.002670e-03
      dolcan t -1 Equation cpiUSA -1.919642e+00
      cpiUSA t -1 Equation cpiUSA  4.389649e-01
      dolcan t -2 Equation cpiUSA -4.998484e+00
      cpiUSA t -2 Equation cpiUSA -2.456765e-01
      ECT Equation cpiUSA          3.550181e-01
      Const Equation cpiUSA       -1.138314e-01
      Trend Equation cpiUSA        1.882903e-03
      dolcan t -1 Equation cpiUSA -1.348386e+00
      cpiUSA t -1 Equation cpiUSA  3.681993e-01
      dolcan t -2 Equation cpiUSA  4.192614e-01
      cpiUSA t -2 Equation cpiUSA  2.529322e-01
      
      [[43]]
                                           [,1]
      ECT Equation dolcan         -0.0026770867
      Const Equation dolcan        0.0019011406
      dolcan t -1 Equation dolcan  0.1821788744
      cpiUSA t -1 Equation dolcan -0.0033491756
      dolcan t -2 Equation dolcan  0.0228978890
      cpiUSA t -2 Equation dolcan -0.0007533507
      ECT Equation dolcan         -0.0019469643
      Const Equation dolcan        0.0041268900
      dolcan t -1 Equation dolcan  0.1173724718
      cpiUSA t -1 Equation dolcan  0.0054390809
      dolcan t -2 Equation dolcan -0.2105055304
      cpiUSA t -2 Equation dolcan -0.0109995510
      ECT Equation cpiUSA         -0.0309453200
      Const Equation cpiUSA        0.1736585917
      dolcan t -1 Equation cpiUSA  0.5033226493
      cpiUSA t -1 Equation cpiUSA  0.3590560466
      dolcan t -2 Equation cpiUSA -1.3437654458
      cpiUSA t -2 Equation cpiUSA -0.1392518022
      ECT Equation cpiUSA         -0.1714883501
      Const Equation cpiUSA        0.1464121394
      dolcan t -1 Equation cpiUSA -0.8026720620
      cpiUSA t -1 Equation cpiUSA  0.4307260063
      dolcan t -2 Equation cpiUSA  1.3484553395
      cpiUSA t -2 Equation cpiUSA  0.3192326008
      
      [[44]]
                                          [,1]
      ECT Equation dolcan         -0.046156798
      Const Equation dolcan       -0.010111647
      dolcan t -1 Equation dolcan  0.158245190
      cpiUSA t -1 Equation dolcan -0.011195528
      dolcan t -2 Equation dolcan  0.018420448
      cpiUSA t -2 Equation dolcan  0.003495803
      ECT Equation dolcan          0.020431500
      Const Equation dolcan       -0.003717539
      dolcan t -1 Equation dolcan  0.186732460
      cpiUSA t -1 Equation dolcan  0.017180057
      dolcan t -2 Equation dolcan  0.019779272
      cpiUSA t -2 Equation dolcan -0.009962532
      ECT Equation dolcan         -0.017644682
      Const Equation dolcan        0.012678289
      dolcan t -1 Equation dolcan  0.120276258
      cpiUSA t -1 Equation dolcan -0.002199864
      dolcan t -2 Equation dolcan -0.216552894
      cpiUSA t -2 Equation dolcan -0.010093165
      ECT Equation cpiUSA          0.960128385
      Const Equation cpiUSA        0.518997756
      dolcan t -1 Equation cpiUSA -0.372275967
      cpiUSA t -1 Equation cpiUSA  0.206833808
      dolcan t -2 Equation cpiUSA -1.602010322
      cpiUSA t -2 Equation cpiUSA -0.161519009
      ECT Equation cpiUSA         -0.154707761
      Const Equation cpiUSA        0.142118177
      dolcan t -1 Equation cpiUSA  1.425814827
      cpiUSA t -1 Equation cpiUSA  0.589460044
      dolcan t -2 Equation cpiUSA -1.180491828
      cpiUSA t -2 Equation cpiUSA -0.135171254
      ECT Equation cpiUSA         -0.141294905
      Const Equation cpiUSA        0.142176901
      dolcan t -1 Equation cpiUSA -0.827768858
      cpiUSA t -1 Equation cpiUSA  0.422812224
      dolcan t -2 Equation cpiUSA  1.168287279
      cpiUSA t -2 Equation cpiUSA  0.297587797
      
      [[45]]
                                          [,1]
      ECT Equation dolcan         -0.017914669
      dolcan t -1 Equation dolcan  0.139173536
      cpiUSA t -1 Equation dolcan -0.013284006
      dolcan t -2 Equation dolcan  0.002174167
      cpiUSA t -2 Equation dolcan  0.001436122
      ECT Equation dolcan          0.006491931
      dolcan t -1 Equation dolcan  0.160034945
      cpiUSA t -1 Equation dolcan  0.007915646
      dolcan t -2 Equation dolcan -0.066086343
      cpiUSA t -2 Equation dolcan -0.010651845
      ECT Equation cpiUSA         -0.565345068
      dolcan t -1 Equation cpiUSA  0.657587752
      cpiUSA t -1 Equation cpiUSA  0.314478775
      dolcan t -2 Equation cpiUSA -0.623941905
      cpiUSA t -2 Equation cpiUSA -0.088049070
      ECT Equation cpiUSA          0.154651435
      dolcan t -1 Equation cpiUSA -0.008908433
      cpiUSA t -1 Equation cpiUSA  0.671157310
      dolcan t -2 Equation cpiUSA -0.915817008
      cpiUSA t -2 Equation cpiUSA  0.131826135
      
      [[46]]
                                          [,1]
      ECT Equation dolcan         -0.016829069
      dolcan t -1 Equation dolcan  0.248791698
      cpiUSA t -1 Equation dolcan -0.009049603
      dolcan t -2 Equation dolcan -0.162203592
      cpiUSA t -2 Equation dolcan -0.000472127
      ECT Equation dolcan         -0.015045190
      dolcan t -1 Equation dolcan  0.033949549
      cpiUSA t -1 Equation dolcan -0.026068751
      dolcan t -2 Equation dolcan  0.226209316
      cpiUSA t -2 Equation dolcan  0.015367971
      ECT Equation dolcan          0.005762436
      dolcan t -1 Equation dolcan  0.172567504
      cpiUSA t -1 Equation dolcan  0.010499939
      dolcan t -2 Equation dolcan -0.080842306
      cpiUSA t -2 Equation dolcan -0.012471780
      ECT Equation cpiUSA         -0.630575270
      dolcan t -1 Equation cpiUSA  0.704853888
      cpiUSA t -1 Equation cpiUSA  0.113211738
      dolcan t -2 Equation cpiUSA  0.715918781
      cpiUSA t -2 Equation cpiUSA -0.144705191
      ECT Equation cpiUSA         -0.709563131
      dolcan t -1 Equation cpiUSA -0.332452599
      cpiUSA t -1 Equation cpiUSA  0.439421564
      dolcan t -2 Equation cpiUSA -1.400513919
      cpiUSA t -2 Equation cpiUSA -0.032866320
      ECT Equation cpiUSA          0.152098278
      dolcan t -1 Equation cpiUSA -0.146886730
      cpiUSA t -1 Equation cpiUSA  0.639630121
      dolcan t -2 Equation cpiUSA -1.132415428
      cpiUSA t -2 Equation cpiUSA  0.164573269
      
      [[47]]
                                           [,1]
      ECT Equation dolcan         -1.594506e-03
      Trend Equation dolcan        5.687498e-06
      dolcan t -1 Equation dolcan  1.856852e-01
      cpiUSA t -1 Equation dolcan -1.987102e-03
      dolcan t -2 Equation dolcan  2.346587e-02
      cpiUSA t -2 Equation dolcan  5.909138e-04
      ECT Equation dolcan          6.034676e-03
      Trend Equation dolcan        7.266438e-06
      dolcan t -1 Equation dolcan  1.190289e-01
      cpiUSA t -1 Equation dolcan  6.094024e-03
      dolcan t -2 Equation dolcan -2.116016e-01
      cpiUSA t -2 Equation dolcan -9.724057e-03
      ECT Equation cpiUSA          1.605637e-01
      Trend Equation cpiUSA        7.799357e-04
      dolcan t -1 Equation cpiUSA  5.882712e-01
      cpiUSA t -1 Equation cpiUSA  4.154928e-01
      dolcan t -2 Equation cpiUSA -1.409047e+00
      cpiUSA t -2 Equation cpiUSA -8.212513e-02
      ECT Equation cpiUSA          1.032605e-01
      Trend Equation cpiUSA        9.694584e-04
      dolcan t -1 Equation cpiUSA -1.040707e+00
      cpiUSA t -1 Equation cpiUSA  4.169815e-01
      dolcan t -2 Equation cpiUSA  9.982535e-01
      cpiUSA t -2 Equation cpiUSA  3.032905e-01
      
      [[48]]
                                           [,1]
      ECT Equation dolcan         -2.995578e-03
      Trend Equation dolcan        6.111664e-06
      dolcan t -1 Equation dolcan  1.616106e-01
      cpiUSA t -1 Equation dolcan -7.636568e-03
      dolcan t -2 Equation dolcan  4.661039e-02
      cpiUSA t -2 Equation dolcan  4.104931e-03
      ECT Equation dolcan         -5.993876e-02
      Trend Equation dolcan        6.599656e-05
      dolcan t -1 Equation dolcan  1.751623e-01
      cpiUSA t -1 Equation dolcan  4.194837e-02
      dolcan t -2 Equation dolcan  2.039924e-01
      cpiUSA t -2 Equation dolcan -3.413151e-02
      ECT Equation dolcan          6.794794e-03
      Trend Equation dolcan        6.522883e-05
      dolcan t -1 Equation dolcan  1.066869e-01
      cpiUSA t -1 Equation dolcan -2.756671e-03
      dolcan t -2 Equation dolcan -2.406763e-01
      cpiUSA t -2 Equation dolcan -1.007585e-02
      ECT Equation cpiUSA          1.748642e-01
      Trend Equation cpiUSA        8.203576e-04
      dolcan t -1 Equation cpiUSA  5.406502e-01
      cpiUSA t -1 Equation cpiUSA  3.552330e-01
      dolcan t -2 Equation cpiUSA -1.332648e+00
      cpiUSA t -2 Equation cpiUSA -4.627725e-02
      ECT Equation cpiUSA          4.228013e+00
      Trend Equation cpiUSA       -4.711408e-03
      dolcan t -1 Equation cpiUSA -1.383088e+00
      cpiUSA t -1 Equation cpiUSA  7.378154e-01
      dolcan t -2 Equation cpiUSA -2.431600e+00
      cpiUSA t -2 Equation cpiUSA -3.809715e-01
      ECT Equation cpiUSA          1.359670e-01
      Trend Equation cpiUSA        1.373138e-03
      dolcan t -1 Equation cpiUSA -1.259400e+00
      cpiUSA t -1 Equation cpiUSA  3.673811e-01
      dolcan t -2 Equation cpiUSA  6.007203e-01
      cpiUSA t -2 Equation cpiUSA  2.467398e-01
      

---

    Code
      sapply(mods_nonLIn, coef, regime = "L")
    Output
      [[1]]
                        Intercept         Trend dolcan -1    cpiUSA -1
      Equation dolcan  0.08210166  0.0007054584 0.9822979 -0.002408619
      Equation cpiUSA -0.55588003 -0.0017977750 0.7181117  1.003209833
      
      [[2]]
                        Intercept         Trend dolcan -1    cpiUSA -1
      Equation dolcan  0.08210166  0.0007054584 0.9822979 -0.002408619
      Equation cpiUSA -0.55588003 -0.0017977750 0.7181117  1.003209833
      
      [[3]]
                         Intercept dolcan -1    cpiUSA -1
      Equation dolcan  0.007691528 0.9939684 1.357332e-05
      Equation cpiUSA -1.012509544 1.3415788 9.969201e-01
      
      [[4]]
                        Intercept dolcan -1    cpiUSA -1
      Equation dolcan -0.02978663 1.0310970 8.444228e-06
      Equation cpiUSA -0.27074633 0.5937531 9.970502e-01
      
      [[5]]
                      dolcan -1     cpiUSA -1
      Equation dolcan 1.0018378 -5.982604e-06
      Equation cpiUSA 0.3056463  9.994945e-01
      
      [[6]]
                      dolcan -1    cpiUSA -1
      Equation dolcan 0.9994742 8.089218e-05
      Equation cpiUSA 0.3063167 9.977088e-01
      
      [[7]]
                             Trend dolcan -1    cpiUSA -1
      Equation dolcan 2.085515e-05 1.0038460 -0.000084888
      Equation cpiUSA 3.944688e-03 0.6854827  0.984569738
      
      [[8]]
                             Trend dolcan -1    cpiUSA -1
      Equation dolcan 0.0002877063 1.0280167 -0.001015842
      Equation cpiUSA 0.0010306707 0.4085665  0.993779851
      
      [[9]]
                        Intercept         Trend   dolcan -1   cpiUSA -1  dolcan -2
      Equation dolcan  0.09231288  0.0006969576  1.11640771 0.005716144 -0.1476973
      Equation cpiUSA -0.38461556 -0.0006441612 -0.05789139 1.209808395  0.6664581
                         cpiUSA -2
      Equation dolcan -0.008053514
      Equation cpiUSA -0.210683537
      
      [[10]]
                        Intercept         Trend   dolcan -1   cpiUSA -1  dolcan -2
      Equation dolcan  0.09231288  0.0006969576  1.11640771 0.005716144 -0.1476973
      Equation cpiUSA -0.38461556 -0.0006441612 -0.05789139 1.209808395  0.6664581
                         cpiUSA -2
      Equation dolcan -0.008053514
      Equation cpiUSA -0.210683537
      
      [[11]]
                        Intercept dolcan -1   cpiUSA -1   dolcan -2   cpiUSA -2
      Equation dolcan  0.01395885 1.1438210 0.002560297 -0.15743623 -0.00252461
      Equation cpiUSA -0.64727526 0.8354151 1.394781489  0.01498696 -0.39677145
      
      [[12]]
                        Intercept dolcan -1   cpiUSA -1   dolcan -2   cpiUSA -2
      Equation dolcan  0.01395885 1.1438210 0.002560297 -0.15743623 -0.00252461
      Equation cpiUSA -0.64727526 0.8354151 1.394781489  0.01498696 -0.39677145
      
      [[13]]
                      dolcan -1     cpiUSA -1  dolcan -2     cpiUSA -2
      Equation dolcan 1.1452641 -0.0005511544 -0.1435657  0.0005498422
      Equation cpiUSA 0.6624965  1.4943161966 -0.5053663 -0.4946045183
      
      [[14]]
                        dolcan -1   cpiUSA -1  dolcan -2    cpiUSA -2
      Equation dolcan  1.19639319 0.003697079 -0.1986477 -0.003599167
      Equation cpiUSA -0.03686374 1.224304032  0.2957565 -0.226656499
      
      [[15]]
                              Trend dolcan -1     cpiUSA -1  dolcan -2     cpiUSA -2
      Equation dolcan -1.622245e-05 1.1450492 -0.0004092175 -0.1449554  0.0004693694
      Equation cpiUSA  3.096566e-03 0.7035062  1.4672230481 -0.2400976 -0.4792437363
      
      [[16]]
                              Trend dolcan -1   cpiUSA -1  dolcan -2    cpiUSA -2
      Equation dolcan -6.292985e-05  1.155769 0.001521385 -0.1604489 -0.001290011
      Equation cpiUSA  3.007284e-03  1.053381 1.507577719 -0.6132165 -0.519312189
      
      [[17]]
                              ECT       Const         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan -0.05165772 0.005288736 -3.885903e-05   0.1675227  -0.0121845
      Equation cpiUSA  0.73785216 0.623866396 -1.062156e-03  -0.1019883   0.1395116
      
      [[18]]
                              ECT       Const         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan -0.05165772 0.005288736 -3.885903e-05   0.1675227  -0.0121845
      Equation cpiUSA  0.73785216 0.623866396 -1.062156e-03  -0.1019883   0.1395116
      
      [[19]]
                              ECT        Const dolcan t -1 cpiUSA t -1
      Equation dolcan -0.03750343 -0.006852286   0.1250896 -0.01024718
      Equation cpiUSA  0.70847089  0.400683554   0.1907474  0.19818232
      
      [[20]]
                              ECT        Const dolcan t -1 cpiUSA t -1
      Equation dolcan -0.03750343 -0.006852286   0.1250896 -0.01024718
      Equation cpiUSA  0.70847089  0.400683554   0.1907474  0.19818232
      
      [[21]]
                              ECT dolcan t -1 cpiUSA t -1
      Equation dolcan -0.02202667    0.117653 -0.01446972
      Equation cpiUSA -0.55070230    1.016250  0.32487254
      
      [[22]]
                              ECT dolcan t -1 cpiUSA t -1
      Equation dolcan -0.02522709   0.1674051 -0.01801245
      Equation cpiUSA -0.45272804   2.2623818  0.24296512
      
      [[23]]
                             ECT         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan -0.0599990 -4.729002e-05   0.2434441  -0.0116853
      Equation cpiUSA -0.5717297  4.019271e-06   1.4443120   0.1073026
      
      [[24]]
                              ECT         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan -0.06592635 -4.893727e-05   0.1448262  -0.0165236
      Equation cpiUSA -0.74852793 -3.444909e-04   2.2377120   0.2191881
      
      [[25]]
                              ECT      Const         Trend dolcan t -1  cpiUSA t -1
      Equation dolcan -0.01518849 0.01105132 -4.790127e-05   0.1838534 -0.004546173
      Equation cpiUSA -0.25113994 0.33469736 -8.430391e-04   0.5327942  0.337989478
                      dolcan t -2  cpiUSA t -2
      Equation dolcan  0.02939966 -0.002047247
      Equation cpiUSA -1.22933741 -0.162023742
      
      [[26]]
                              ECT      Const         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan -0.01739221 0.01123809 -4.874173e-05   0.1570717 -0.01032538
      Equation cpiUSA -0.18905435 0.28407694 -5.662286e-04   0.4259146  0.28726509
                      dolcan t -2  cpiUSA t -2
      Equation dolcan  0.05235467  0.001209455
      Equation cpiUSA -1.18744412 -0.119469224
      
      [[27]]
                               ECT       Const dolcan t -1  cpiUSA t -1 dolcan t -2
      Equation dolcan -0.002677087 0.001901141   0.1821789 -0.003349176  0.02289789
      Equation cpiUSA -0.030945320 0.173658592   0.5033226  0.359056047 -1.34376545
                        cpiUSA t -2
      Equation dolcan -0.0007533507
      Equation cpiUSA -0.1392518022
      
      [[28]]
                             ECT       Const dolcan t -1 cpiUSA t -1 dolcan t -2
      Equation dolcan -0.0461568 -0.01011165   0.1582452 -0.01119553  0.01842045
      Equation cpiUSA  0.9601284  0.51899776  -0.3722760  0.20683381 -1.60201032
                       cpiUSA t -2
      Equation dolcan  0.003495803
      Equation cpiUSA -0.161519009
      
      [[29]]
                              ECT dolcan t -1 cpiUSA t -1  dolcan t -2  cpiUSA t -2
      Equation dolcan -0.01791467   0.1391735 -0.01328401  0.002174167  0.001436122
      Equation cpiUSA -0.56534507   0.6575878  0.31447878 -0.623941905 -0.088049070
      
      [[30]]
                              ECT dolcan t -1  cpiUSA t -1 dolcan t -2  cpiUSA t -2
      Equation dolcan -0.01682907   0.2487917 -0.009049603  -0.1622036 -0.000472127
      Equation cpiUSA -0.63057527   0.7048539  0.113211738   0.7159188 -0.144705191
      
      [[31]]
                               ECT        Trend dolcan t -1  cpiUSA t -1 dolcan t -2
      Equation dolcan -0.001594506 5.687498e-06   0.1856852 -0.001987102  0.02346587
      Equation cpiUSA  0.160563744 7.799357e-04   0.5882712  0.415492814 -1.40904656
                        cpiUSA t -2
      Equation dolcan  0.0005909138
      Equation cpiUSA -0.0821251317
      
      [[32]]
                               ECT        Trend dolcan t -1  cpiUSA t -1 dolcan t -2
      Equation dolcan -0.002995578 6.111664e-06   0.1616106 -0.007636568  0.04661039
      Equation cpiUSA  0.174864230 8.203576e-04   0.5406502  0.355232986 -1.33264807
                       cpiUSA t -2
      Equation dolcan  0.004104931
      Equation cpiUSA -0.046277248
      

---

    Code
      sapply(mods_nonLIn, coef, regime = "H")
    Output
      [[1]]
                        Intercept         Trend  dolcan -1    cpiUSA -1
      Equation dolcan -0.01354042 -0.0001653546  0.9914561 0.0007168335
      Equation cpiUSA  2.20206986  0.0073398658 -0.7468369 0.9702240258
      
      [[2]]
                       Intercept        Trend dolcan -1    cpiUSA -1
      Equation dolcan  0.3579121  0.001426282 0.8674328 -0.005605241
      Equation cpiUSA -1.2735570 -0.003504812 0.1599709  1.021491583
      
      [[3]]
                       Intercept  dolcan -1   cpiUSA -1
      Equation dolcan 0.03045117  0.9689214 0.000145664
      Equation cpiUSA 0.46501765 -0.2412925 1.000711673
      
      [[4]]
                       Intercept  dolcan -1   cpiUSA -1
      Equation dolcan 0.03045117  0.9689214 0.000145664
      Equation cpiUSA 0.46501765 -0.2412925 1.000711673
      
      [[5]]
                      dolcan -1    cpiUSA -1
      Equation dolcan 0.9949031 9.159023e-05
      Equation cpiUSA 0.1554722 9.998859e-01
      
      [[6]]
                      dolcan -1    cpiUSA -1
      Equation dolcan 0.9949031 9.159023e-05
      Equation cpiUSA 0.1554722 9.998859e-01
      
      [[7]]
                              Trend    dolcan -1   cpiUSA -1
      Equation dolcan -0.0002075003  0.970279328 0.000997695
      Equation cpiUSA -0.0013705132 -0.007164849 1.005870624
      
      [[8]]
                              Trend    dolcan -1   cpiUSA -1
      Equation dolcan -0.0002075003  0.970279328 0.000997695
      Equation cpiUSA -0.0013705132 -0.007164849 1.005870624
      
      [[9]]
                         Intercept         Trend dolcan -1    cpiUSA -1  dolcan -2
      Equation dolcan 0.0006464748 -0.0001118549  1.141001 -0.003541846 -0.1549285
      Equation cpiUSA 1.6398792639  0.0059730476  0.572203  1.334289684 -1.1084546
                        cpiUSA -2
      Equation dolcan  0.00404909
      Equation cpiUSA -0.35828351
      
      [[10]]
                      Intercept       Trend dolcan -1   cpiUSA -1   dolcan -2
      Equation dolcan 0.3044483 0.001169578 0.9527351 -0.02974125 -0.07704594
      Equation cpiUSA 1.2815143 0.009182135 2.0104789  1.30743327 -2.39121154
                        cpiUSA -2
      Equation dolcan  0.02531077
      Equation cpiUSA -0.33933921
      
      [[11]]
                       Intercept dolcan -1   cpiUSA -1  dolcan -2   cpiUSA -2
      Equation dolcan 0.04017077 1.1335259 -0.01023868 -0.1713070  0.01039812
      Equation cpiUSA 0.28300510 0.4492155  1.31096745 -0.5821929 -0.31059283
      
      [[12]]
                        Intercept dolcan -1   cpiUSA -1   dolcan -2   cpiUSA -2
      Equation dolcan  0.05798003 0.9631114 -0.02791145 -0.02695616  0.02829938
      Equation cpiUSA -0.65346050 2.0919411  1.32179865 -1.99796630 -0.31587620
      
      [[13]]
                      dolcan -1    cpiUSA -1  dolcan -2    cpiUSA -2
      Equation dolcan 1.1609381 -0.007105683 -0.1655764  0.007202736
      Equation cpiUSA 0.5957397  1.314457492 -0.4870101 -0.314595092
      
      [[14]]
                      dolcan -1    cpiUSA -1  dolcan -2    cpiUSA -2
      Equation dolcan 1.1589176 -0.008044606 -0.1627678  0.008132882
      Equation cpiUSA 0.6281012  1.326424793 -0.5220339 -0.326551571
      
      [[15]]
                              Trend dolcan -1    cpiUSA -1  dolcan -2   cpiUSA -2
      Equation dolcan -0.0002216198 1.1270540 -0.009180805 -0.1574311  0.01024224
      Equation cpiUSA -0.0005662404 0.5091657  1.309155537 -0.4661987 -0.30682914
      
      [[16]]
                              Trend dolcan -1   cpiUSA -1  dolcan -2   cpiUSA -2
      Equation dolcan -0.0001633523 0.9858677 -0.02757334 -0.0233897  0.02860914
      Equation cpiUSA  0.0035714290 2.1499439  1.31655866 -2.1653563 -0.32545538
      
      [[17]]
                               ECT       Const         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan -0.001023032 0.003613178 -3.313512e-05   0.1312604 0.002326573
      Equation cpiUSA  0.015485020 0.101293958 -3.625345e-05   0.2371179 0.620902189
      
      [[18]]
                              ECT       Const        Trend dolcan t -1  cpiUSA t -1
      Equation dolcan -0.02900588  0.01680242 3.615931e-06   0.1850201 -0.007791395
      Equation cpiUSA  0.33504178 -0.12199902 1.994943e-03  -1.5209576  0.614156675
      
      [[19]]
                              ECT        Const dolcan t -1 cpiUSA t -1
      Equation dolcan 0.006057814 -0.001555704   0.1644276 0.002927573
      Equation cpiUSA 0.058696937  0.085773731   0.2256604 0.619010979
      
      [[20]]
                                ECT        Const dolcan t -1 cpiUSA t -1
      Equation dolcan -0.0007475249 0.0007799121   0.1052364 0.002985671
      Equation cpiUSA -0.0967219539 0.1360188778  -0.7272839 0.680080621
      
      [[21]]
                             ECT dolcan t -1 cpiUSA t -1
      Equation dolcan 0.00279744   0.1680402 0.001035012
      Equation cpiUSA 0.18971894  -0.2859286 0.750443210
      
      [[22]]
                              ECT dolcan t -1  cpiUSA t -1
      Equation dolcan 0.003049438   0.1689193 0.0007218864
      Equation cpiUSA 0.194655833  -0.2687060 0.7443087565
      
      [[23]]
                              ECT        Trend dolcan t -1  cpiUSA t -1
      Equation dolcan 0.002438058 9.075573e-08   0.1271129 0.0007793217
      Equation cpiUSA 0.188165815 5.901968e-04  -0.1180760 0.5950325178
      
      [[24]]
                              ECT         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan 0.004302129 -0.0000163670  0.08396771 0.001005123
      Equation cpiUSA 0.129046087  0.0007655733 -0.67724374 0.707856779
      
      [[25]]
                              ECT       Const         Trend dolcan t -1 cpiUSA t -1
      Equation dolcan -0.01284187  0.00991402 -3.535269e-05    0.125544 0.005828206
      Equation cpiUSA  0.18220945 -0.04146424  1.147708e-03   -1.067956 0.418093241
                      dolcan t -2 cpiUSA t -2
      Equation dolcan  -0.1979610 -0.01062449
      Equation cpiUSA   0.9412034  0.30705642
      
      [[26]]
                             ECT        Const        Trend dolcan t -1  cpiUSA t -1
      Equation dolcan -0.0119483  0.009739977 2.161084e-05   0.1143009 -0.002826679
      Equation cpiUSA  0.3550181 -0.113831430 1.882903e-03  -1.3483856  0.368199325
                      dolcan t -2 cpiUSA t -2
      Equation dolcan  -0.2251498  -0.0106057
      Equation cpiUSA   0.4192614   0.2529322
      
      [[27]]
                               ECT      Const dolcan t -1 cpiUSA t -1 dolcan t -2
      Equation dolcan -0.001946964 0.00412689   0.1173725 0.005439081  -0.2105055
      Equation cpiUSA -0.171488350 0.14641214  -0.8026721 0.430726006   1.3484553
                      cpiUSA t -2
      Equation dolcan -0.01099955
      Equation cpiUSA  0.31923260
      
      [[28]]
                              ECT      Const dolcan t -1  cpiUSA t -1 dolcan t -2
      Equation dolcan -0.01764468 0.01267829   0.1202763 -0.002199864  -0.2165529
      Equation cpiUSA -0.14129491 0.14217690  -0.8277689  0.422812224   1.1682873
                      cpiUSA t -2
      Equation dolcan -0.01009317
      Equation cpiUSA  0.29758780
      
      [[29]]
                              ECT  dolcan t -1 cpiUSA t -1 dolcan t -2 cpiUSA t -2
      Equation dolcan 0.006491931  0.160034945 0.007915646 -0.06608634 -0.01065185
      Equation cpiUSA 0.154651435 -0.008908433 0.671157310 -0.91581701  0.13182613
      
      [[30]]
                              ECT dolcan t -1 cpiUSA t -1 dolcan t -2 cpiUSA t -2
      Equation dolcan 0.005762436   0.1725675  0.01049994 -0.08084231 -0.01247178
      Equation cpiUSA 0.152098278  -0.1468867  0.63963012 -1.13241543  0.16457327
      
      [[31]]
                              ECT        Trend dolcan t -1 cpiUSA t -1 dolcan t -2
      Equation dolcan 0.006034676 7.266438e-06   0.1190289 0.006094024  -0.2116016
      Equation cpiUSA 0.103260490 9.694584e-04  -1.0407069 0.416981490   0.9982535
                       cpiUSA t -2
      Equation dolcan -0.009724057
      Equation cpiUSA  0.303290461
      
      [[32]]
                              ECT        Trend dolcan t -1  cpiUSA t -1 dolcan t -2
      Equation dolcan 0.006794794 6.522883e-05   0.1066869 -0.002756671  -0.2406763
      Equation cpiUSA 0.135966975 1.373138e-03  -1.2593996  0.367381138   0.6007203
                      cpiUSA t -2
      Equation dolcan -0.01007585
      Equation cpiUSA  0.24673984
      

# utilities: get_lag, get_nVar, df.residual

    Code
      sapply(mods, tsDyn:::get_series)
    Output
           [,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]     [,8]    
      [1,] "dolcan" "dolcan" "dolcan" "dolcan" "dolcan" "dolcan" "dolcan" "dolcan"
      [2,] "cpiUSA" "cpiUSA" "cpiUSA" "cpiUSA" "cpiUSA" "cpiUSA" "cpiUSA" "cpiUSA"
           [,9]     [,10]    [,11]    [,12]    [,13]    [,14]    [,15]    [,16]   
      [1,] "dolcan" "dolcan" "dolcan" "dolcan" "dolcan" "dolcan" "dolcan" "dolcan"
      [2,] "cpiUSA" "cpiUSA" "cpiUSA" "cpiUSA" "cpiUSA" "cpiUSA" "cpiUSA" "cpiUSA"
           [,17]    [,18]    [,19]    [,20]    [,21]    [,22]    [,23]    [,24]   
      [1,] "dolcan" "dolcan" "dolcan" "dolcan" "dolcan" "dolcan" "dolcan" "dolcan"
      [2,] "cpiUSA" "cpiUSA" "cpiUSA" "cpiUSA" "cpiUSA" "cpiUSA" "cpiUSA" "cpiUSA"
           [,25]    [,26]    [,27]    [,28]    [,29]    [,30]    [,31]    [,32]   
      [1,] "dolcan" "dolcan" "dolcan" "dolcan" "dolcan" "dolcan" "dolcan" "dolcan"
      [2,] "cpiUSA" "cpiUSA" "cpiUSA" "cpiUSA" "cpiUSA" "cpiUSA" "cpiUSA" "cpiUSA"
           [,33]    [,34]    [,35]    [,36]    [,37]    [,38]    [,39]    [,40]   
      [1,] "dolcan" "dolcan" "dolcan" "dolcan" "dolcan" "dolcan" "dolcan" "dolcan"
      [2,] "cpiUSA" "cpiUSA" "cpiUSA" "cpiUSA" "cpiUSA" "cpiUSA" "cpiUSA" "cpiUSA"
           [,41]    [,42]    [,43]    [,44]    [,45]    [,46]    [,47]    [,48]   
      [1,] "dolcan" "dolcan" "dolcan" "dolcan" "dolcan" "dolcan" "dolcan" "dolcan"
      [2,] "cpiUSA" "cpiUSA" "cpiUSA" "cpiUSA" "cpiUSA" "cpiUSA" "cpiUSA" "cpiUSA"

---

    Code
      models_multivariate %>% mutate(lag_out = map_int(object, tsDyn:::get_lag.nlVar),
      nVar_out = map_int(object, tsDyn:::get_nVar.nlVar), df_out = map(object,
        df.residual)) %>% select(-starts_with("object")) %>% as.data.frame()
    Output
         lag include model nthresh lag_out nVar_out df_out
      1    1    both   VAR      NA       1        2    319
      2    1   const   VAR      NA       1        2    320
      3    1    none   VAR      NA       1        2    321
      4    1   trend   VAR      NA       1        2    320
      5    2    both   VAR      NA       2        2    316
      6    2   const   VAR      NA       2        2    317
      7    2    none   VAR      NA       2        2    318
      8    2   trend   VAR      NA       2        2    317
      9    1    both  VECM      NA       1        2    317
      10   1   const  VECM      NA       1        2    318
      11   1    none  VECM      NA       1        2    319
      12   1   trend  VECM      NA       1        2    318
      13   2    both  VECM      NA       2        2    314
      14   2   const  VECM      NA       2        2    315
      15   2    none  VECM      NA       2        2    316
      16   2   trend  VECM      NA       2        2    315
      17   1    both  TVAR       1       1        2    315
      18   1    both  TVAR       2       1        2    311
      19   1   const  TVAR       1       1        2    317
      20   1   const  TVAR       2       1        2    314
      21   1    none  TVAR       1       1        2    319
      22   1    none  TVAR       2       1        2    317
      23   1   trend  TVAR       1       1        2    317
      24   1   trend  TVAR       2       1        2    314
      25   2    both  TVAR       1       2        2    310
      26   2    both  TVAR       2       2        2    304
      27   2   const  TVAR       1       2        2    312
      28   2   const  TVAR       2       2        2    307
      29   2    none  TVAR       1       2        2    314
      30   2    none  TVAR       2       2        2    310
      31   2   trend  TVAR       1       2        2    312
      32   2   trend  TVAR       2       2        2    307
      33   1    both TVECM       1       1        2    317
      34   1    both TVECM       2       1        2    317
      35   1   const TVECM       1       1        2    318
      36   1   const TVECM       2       1        2    318
      37   1    none TVECM       1       1        2    319
      38   1    none TVECM       2       1        2    319
      39   1   trend TVECM       1       1        2    318
      40   1   trend TVECM       2       1        2    318
      41   2    both TVECM       1       2        2    314
      42   2    both TVECM       2       2        2    314
      43   2   const TVECM       1       2        2    315
      44   2   const TVECM       2       2        2    315
      45   2    none TVECM       1       2        2    316
      46   2    none TVECM       2       2        2    316
      47   2   trend TVECM       1       2        2    315
      48   2   trend TVECM       2       2        2    315

# AIC et al

    Code
      as.data.frame(uni_stats)
    Output
         lag include model nthresh  deviance       AIC       BIC     logLik
      1    1    both   VAR      NA 10.044863 -3918.837 -3888.616 1050.78416
      2    1   const   VAR      NA 10.056006 -3922.400 -3899.734 1050.56567
      3    1    none   VAR      NA 10.572927 -3907.237 -3892.126 1040.98396
      4    1   trend   VAR      NA 10.200193 -3916.979 -3894.313 1047.85541
      5    2    both   VAR      NA  7.949829 -3982.323 -3937.029 1089.36528
      6    2   const   VAR      NA  7.952476 -3986.131 -3948.385 1089.26900
      7    2    none   VAR      NA  8.095284 -3980.289 -3950.092 1084.34791
      8    2   trend   VAR      NA  8.034911 -3981.552 -3943.806 1086.97953
      9    1    both  VECM      NA  7.960313 -3983.832 -3942.312  -12.44047
      10   1   const  VECM      NA  7.984527 -3983.729 -3949.758  -12.42773
      11   1    none  VECM      NA  9.764497 -3919.627 -3893.205  -12.21623
      12   1   trend  VECM      NA  8.181358 -3974.330 -3940.359  -12.39854
      13   2    both  VECM      NA  7.937843 -3962.912 -3906.340  -12.43898
      14   2   const  VECM      NA  7.956742 -3962.974 -3913.945  -12.42671
      15   2    none  VECM      NA  9.344372 -3911.881 -3870.395  -12.25508
      16   2   trend  VECM      NA  8.131732 -3954.238 -3905.209  -12.39950
      17   1    both  TVAR       1  8.161472 -3974.104 -3909.884 1087.41751
      18   1    both  TVAR       2  7.794644 -3974.970 -3876.751 1096.85081
      19   1   const  TVAR       1  8.348522 -3971.026 -3921.917 1081.87871
      20   1   const  TVAR       2  7.915594 -3976.778 -3901.225 1091.75449
      21   1    none  TVAR       1  9.029764 -3950.762 -3916.763 1067.74675
      22   1    none  TVAR       2  8.123147 -3975.972 -3923.085 1085.35180
      23   1   trend  TVAR       1  8.842832 -3954.841 -3905.732 1073.78631
      24   1   trend  TVAR       2  8.069110 -3974.302 -3898.749 1090.51655
      25   2    both  TVAR       1  7.045296 -3999.979 -3905.615 1111.19308
      26   2    both  TVAR       2  6.794512 -3995.990 -3852.557 1122.19836
      27   2   const  TVAR       1  7.278526 -3997.421 -3918.155 1105.91411
      28   2   const  TVAR       2  6.972166 -3995.703 -3874.917 1116.05490
      29   2    none  TVAR       1  7.512746 -3990.483 -3926.315 1098.44488
      30   2    none  TVAR       2  7.141784 -3990.169 -3892.031 1107.28805
      31   2   trend  TVAR       1  7.415455 -3992.260 -3912.995 1103.33370
      32   2   trend  TVAR       2  7.025195 -3994.773 -3873.988 1115.59023
      33   1    both TVECM       1  7.298715 -4002.429 -3923.163 1108.41789
      34   1    both TVECM       2  7.007150 -4001.215 -3880.429 1118.81102
      35   1   const TVECM       1  7.411220 -4001.207 -3937.040 1103.80725
      36   1   const TVECM       2  7.128824 -4000.627 -3902.489 1112.51714
      37   1    none TVECM       1  8.235945 -3974.743 -3925.674 1086.57525
      38   1    none TVECM       2  7.768609 -3980.323 -3904.831 1096.36484
      39   1   trend TVECM       1  7.674577 -3986.448 -3922.280 1096.42744
      40   1   trend TVECM       2  7.313595 -3985.205 -3887.067 1104.80603
      41   2    both TVECM       1  7.026082 -3978.469 -3869.098 1107.27617
      42   2    both TVECM       2  6.526916 -3987.609 -3821.666 1126.84620
      43   2   const TVECM       1  7.174833 -3977.099 -3882.813 1102.59081
      44   2   const TVECM       2  6.958434 -3977.753 -3834.438 1115.91804
      45   2    none TVECM       1  8.129782 -3951.727 -3872.527 1085.90505
      46   2    none TVECM       2  7.778024 -3957.032 -3836.346 1099.55728
      47   2   trend TVECM       1  7.476702 -3962.656 -3868.370 1095.36946
      48   2   trend TVECM       2  7.116788 -3965.483 -3822.168 1109.78272

# More methods

    Code
      sapply(mods, function(x) dim(residuals(x, initVal = FALSE)))
    Output
           [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]
      [1,]  323  323  323  323  322  322  322  322  322   322   322   322   321   321
      [2,]    2    2    2    2    2    2    2    2    2     2     2     2     2     2
           [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24] [,25] [,26]
      [1,]   321   321   323   323   323   323   323   323   323   323   322   322
      [2,]     2     2     2     2     2     2     2     2     2     2     2     2
           [,27] [,28] [,29] [,30] [,31] [,32] [,33] [,34] [,35] [,36] [,37] [,38]
      [1,]   322   322   322   322   322   322   322   322   322   322   322   322
      [2,]     2     2     2     2     2     2     2     2     2     2     2     2
           [,39] [,40] [,41] [,42] [,43] [,44] [,45] [,46] [,47] [,48]
      [1,]   322   322   321   321   321   321   321   321   321   321
      [2,]     2     2     2     2     2     2     2     2     2     2

---

    Code
      sapply(mods, function(x) dim(residuals(x, initVal = TRUE)))
    Output
           [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]
      [1,]  324  324  324  324  324  324  324  324  324   324   324   324   324   324
      [2,]    2    2    2    2    2    2    2    2    2     2     2     2     2     2
           [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24] [,25] [,26]
      [1,]   324   324   324   324   324   324   324   324   324   324   324   324
      [2,]     2     2     2     2     2     2     2     2     2     2     2     2
           [,27] [,28] [,29] [,30] [,31] [,32] [,33] [,34] [,35] [,36] [,37] [,38]
      [1,]   324   324   324   324   324   324   324   324   324   324   324   324
      [2,]     2     2     2     2     2     2     2     2     2     2     2     2
           [,39] [,40] [,41] [,42] [,43] [,44] [,45] [,46] [,47] [,48]
      [1,]   324   324   324   324   324   324   324   324   324   324
      [2,]     2     2     2     2     2     2     2     2     2     2

---

    Code
      sapply(mods, function(x) head(residuals(x), 3))
    Output
                   [,1]          [,2]          [,3]          [,4]          [,5]
      [1,] -0.005308679 -0.0053528912 -0.0044470911 -4.875801e-03 -0.0002761132
      [2,] -0.000562297 -0.0006130099  0.0003521491 -9.193893e-05  0.0021370882
      [3,]  0.002334327  0.0022838251  0.0032484540  2.803963e-03 -0.0025009876
      [4,] -0.112614651 -0.1140590044 -0.0817844594 -9.690370e-02  0.0176156297
      [5,] -0.032854255 -0.0345109777 -0.0001214182 -1.578299e-02 -0.1018553909
      [6,] -0.112606385 -0.1142562095 -0.0798855362 -9.556132e-02 -0.0766399682
                    [,6]         [,7]          [,8]          [,9]         [,10]
      [1,] -0.0003212195  0.001032152  0.0003538675 -0.0001994337  0.0007452142
      [2,]  0.0020802604  0.003316316  0.0027547418  0.0022212265  0.0031165244
      [3,] -0.0025557374 -0.001209697 -0.0018476732 -0.0024170157 -0.0014856339
      [4,]  0.0182193343  0.038581166  0.0307258957  0.0202515098  0.0268066443
      [5,] -0.1010948050 -0.082498017 -0.0890016583 -0.0989631161 -0.0927504322
      [6,] -0.0759071945 -0.055655648 -0.0630441145 -0.0737534121 -0.0672903334
                   [,11]        [,12]        [,13]        [,14]        [,15]
      [1,]  0.0016321643  0.001168065  0.001814330  0.002760607  0.003482956
      [2,]  0.0035287153  0.003389858 -0.002362168 -0.001434433 -0.000770474
      [3,] -0.0006638402 -0.001084396 -0.004412278 -0.003414396 -0.002414532
      [4,]  0.0817506885  0.044282550 -0.100555090 -0.094831405 -0.056686869
      [5,] -0.0672163725 -0.078426778 -0.073791351 -0.068179817 -0.033118638
      [6,] -0.0163825460 -0.050335303 -0.043677558 -0.037641732  0.015157327
                  [,16]        [,17]        [,18]         [,19]         [,20]
      [1,]  0.003166844 -0.001424210 -0.001424210 -0.0057438958 -0.0052248731
      [2,] -0.001050438  0.003062471  0.003062471 -0.0009687918 -0.0003114186
      [3,] -0.002909126  0.005926812  0.005926812  0.0019341782  0.0025520947
      [4,] -0.078928511 -0.059742392 -0.059742392 -0.0519551103 -0.0500568651
      [5,] -0.052816895  0.024102528  0.024102528  0.0335939086  0.0327004754
      [6,] -0.019642286 -0.055756275 -0.055756275 -0.0470502604 -0.0471562181
                   [,21]         [,22]         [,23]        [,24]        [,25]
      [1,] -0.0053691455 -0.0054379213 -0.0051889773 -0.003561525  0.004102126
      [2,] -0.0005612088 -0.0006552363 -0.0003794736  0.001247440  0.005617448
      [3,]  0.0023383849  0.0022235012  0.0025183605  0.004103193  0.001805708
      [4,] -0.1012870876 -0.0519933722 -0.0672088253 -0.045271423  0.031527658
      [5,] -0.0200601424  0.0295753355  0.0143145229  0.036391430 -0.061917000
      [6,] -0.1002598556 -0.0501429777 -0.0662180249 -0.043409223 -0.054753418
                  [,26]         [,27]         [,28]         [,29]         [,30]
      [1,]  0.004102126 -0.0002059923 -0.0002059923  8.170417e-05  0.0006375665
      [2,]  0.005617448  0.0017416863  0.0017416863  2.335062e-03  0.0022721680
      [3,]  0.001805708 -0.0023652295 -0.0023652295 -2.031810e-03 -0.0019255622
      [4,]  0.031527658  0.0512300076  0.0512300076  2.957823e-02  0.0343563228
      [5,] -0.061917000 -0.0608379247 -0.0608379247 -9.251089e-02 -0.0620061819
      [6,] -0.054753418 -0.0420762634 -0.0420762634 -6.498183e-02 -0.0536045705
                   [,31]         [,32]         [,33]        [,34]         [,35]
      [1,] -4.617848e-05 -0.0003474697 -0.0018647997  0.002666191 -0.0004805491
      [2,]  2.195463e-03  0.0016847810  0.0002490016  0.005208250  0.0014165758
      [3,] -2.158039e-03 -0.0025396333 -0.0038113147  0.000175325 -0.0027311159
      [4,]  5.398867e-02  0.0629246852  0.0429958251  0.080479236  0.0365552931
      [5,] -6.586403e-02 -0.0618868869 -0.0877304605 -0.042313378 -0.0938493139
      [6,] -4.088708e-02 -0.0346586317 -0.0487259889 -0.001118014 -0.0550318552
                   [,36]         [,37]         [,38]         [,39]        [,40]
      [1,]  0.0005356284  1.384438e-05 -5.954496e-05  0.0001118291 -0.001076092
      [2,]  0.0026890999  2.033583e-03  1.981847e-03  0.0023472748  0.001370979
      [3,] -0.0012766080 -2.273531e-03 -2.353100e-03 -0.0018613160 -0.002679665
      [4,]  0.0529823134  2.386660e-02  2.242882e-02  0.0547875865  0.062880401
      [5,] -0.0782618171 -1.141770e-01 -1.151906e-01 -0.0722191478 -0.070838222
      [6,] -0.0316654940 -6.348857e-02 -6.504743e-02 -0.0350355286 -0.023128661
                   [,41]         [,42]         [,43]         [,44]         [,45]
      [1,]  0.0004232282  0.0023172126  0.0006064402  0.0021432667  4.087696e-05
      [2,] -0.0015282666 -0.0002227808 -0.0013067853 -0.0004197752 -2.526233e-03
      [3,] -0.0037580293 -0.0025569680 -0.0035697717 -0.0027031857 -4.964813e-03
      [4,] -0.0317054075 -0.0284696452 -0.0376532966 -0.0436251449 -1.012334e-01
      [5,] -0.0355441350 -0.0295319264 -0.0427344185 -0.0466955940 -6.379911e-02
      [6,]  0.0056861494  0.0070541105 -0.0004255431 -0.0056854926 -2.347663e-02
                   [,46]         [,47]         [,48]
      [1,]  5.409569e-06  0.0001197018  0.0020260673
      [2,] -2.174246e-03 -0.0018566902 -0.0005134883
      [3,] -4.643001e-03 -0.0039856278 -0.0027786570
      [4,] -9.824560e-02 -0.0304359439 -0.0250670208
      [5,] -6.451122e-02 -0.0341705412 -0.0261344179
      [6,] -2.180155e-02  0.0066380537  0.0096449964

---

    Code
      sapply(mods, function(x) head(residuals(x, initVal = TRUE), 3))
    Output
                   [,1]          [,2]          [,3]          [,4]          [,5]
      [1,]           NA            NA            NA            NA            NA
      [2,] -0.005308679 -0.0053528912 -0.0044470911 -4.875801e-03            NA
      [3,] -0.000562297 -0.0006130099  0.0003521491 -9.193893e-05 -0.0002761132
      [4,]           NA            NA            NA            NA            NA
      [5,] -0.112614651 -0.1140590044 -0.0817844594 -9.690370e-02            NA
      [6,] -0.032854255 -0.0345109777 -0.0001214182 -1.578299e-02  0.0176156297
                    [,6]        [,7]         [,8]          [,9]        [,10]
      [1,]            NA          NA           NA            NA           NA
      [2,]            NA          NA           NA            NA           NA
      [3,] -0.0003212195 0.001032152 0.0003538675 -0.0001994337 0.0007452142
      [4,]            NA          NA           NA            NA           NA
      [5,]            NA          NA           NA            NA           NA
      [6,]  0.0182193343 0.038581166 0.0307258957  0.0202515098 0.0268066443
                 [,11]       [,12] [,13] [,14] [,15] [,16]        [,17]        [,18]
      [1,]          NA          NA    NA    NA    NA    NA           NA           NA
      [2,]          NA          NA    NA    NA    NA    NA -0.001424210 -0.001424210
      [3,] 0.001632164 0.001168065    NA    NA    NA    NA  0.003062471  0.003062471
      [4,]          NA          NA    NA    NA    NA    NA           NA           NA
      [5,]          NA          NA    NA    NA    NA    NA -0.059742392 -0.059742392
      [6,] 0.081750688 0.044282550    NA    NA    NA    NA  0.024102528  0.024102528
                   [,19]         [,20]         [,21]         [,22]         [,23]
      [1,]            NA            NA            NA            NA            NA
      [2,] -0.0057438958 -0.0052248731 -0.0053691455 -0.0054379213 -0.0051889773
      [3,] -0.0009687918 -0.0003114186 -0.0005612088 -0.0006552363 -0.0003794736
      [4,]            NA            NA            NA            NA            NA
      [5,] -0.0519551103 -0.0500568651 -0.1012870876 -0.0519933722 -0.0672088253
      [6,]  0.0335939086  0.0327004754 -0.0200601424  0.0295753355  0.0143145229
                  [,24]       [,25]       [,26]         [,27]         [,28]
      [1,]           NA          NA          NA            NA            NA
      [2,] -0.003561525          NA          NA            NA            NA
      [3,]  0.001247440 0.004102126 0.004102126 -0.0002059923 -0.0002059923
      [4,]           NA          NA          NA            NA            NA
      [5,] -0.045271423          NA          NA            NA            NA
      [6,]  0.036391430 0.031527658 0.031527658  0.0512300076  0.0512300076
                  [,29]        [,30]         [,31]         [,32]       [,33]
      [1,]           NA           NA            NA            NA          NA
      [2,]           NA           NA            NA            NA          NA
      [3,] 8.170417e-05 0.0006375665 -4.617848e-05 -0.0003474697 -0.00186480
      [4,]           NA           NA            NA            NA          NA
      [5,]           NA           NA            NA            NA          NA
      [6,] 2.957823e-02 0.0343563228  5.398867e-02  0.0629246852  0.04299583
                 [,34]         [,35]        [,36]        [,37]         [,38]
      [1,]          NA            NA           NA           NA            NA
      [2,]          NA            NA           NA           NA            NA
      [3,] 0.002666191 -0.0004805491 0.0005356284 1.384438e-05 -5.954496e-05
      [4,]          NA            NA           NA           NA            NA
      [5,]          NA            NA           NA           NA            NA
      [6,] 0.080479236  0.0365552931 0.0529823134 2.386660e-02  2.242882e-02
                  [,39]        [,40] [,41] [,42] [,43] [,44] [,45] [,46] [,47] [,48]
      [1,]           NA           NA    NA    NA    NA    NA    NA    NA    NA    NA
      [2,]           NA           NA    NA    NA    NA    NA    NA    NA    NA    NA
      [3,] 0.0001118291 -0.001076092    NA    NA    NA    NA    NA    NA    NA    NA
      [4,]           NA           NA    NA    NA    NA    NA    NA    NA    NA    NA
      [5,]           NA           NA    NA    NA    NA    NA    NA    NA    NA    NA
      [6,] 0.0547875865  0.062880401    NA    NA    NA    NA    NA    NA    NA    NA

---

    Code
      sapply(mods, function(x) tail(residuals(x), 3))
    Output
                    [,1]          [,2]         [,3]          [,4]         [,5]
      [1,]  0.0000438739 -0.0002672998 -0.001430829  0.0004183927  0.004408037
      [2,] -0.0106570524 -0.0109737178 -0.012133437 -0.0102672949 -0.010300829
      [3,]  0.0042229912  0.0038897905  0.002869918  0.0047041355  0.005821738
      [4,]  0.0038329653 -0.0063326719 -0.047790370  0.0174258190 -0.142307519
      [5,] -0.1258433297 -0.1361883725 -0.177510303 -0.1116973985 -0.135430154
      [6,] -0.1967164630 -0.2076016940 -0.243940747 -0.1792537242 -0.140241858
                   [,6]         [,7]         [,8]         [,9]        [,10]
      [1,]  0.004002238  0.001911235  0.004720864  0.004635935  0.002210391
      [2,] -0.010676874 -0.012076582 -0.009809016 -0.010080144 -0.012103246
      [3,]  0.005446472  0.004390315  0.006440390  0.006049157  0.004200300
      [4,] -0.136876284 -0.168335980 -0.135797413 -0.134473459 -0.151304887
      [5,] -0.130397150 -0.151456121 -0.125195241 -0.127844025 -0.141882807
      [6,] -0.135219285 -0.151109458 -0.127367349 -0.132424257 -0.145253914
                  [,11]        [,12]        [,13]        [,14]        [,15]
      [1,]  0.001761170  0.001668408  0.004448133  0.002105365  0.002063668
      [2,] -0.010954132 -0.012102074 -0.009489065 -0.011935364 -0.011837849
      [3,]  0.006152866  0.004469782  0.006226286  0.004324666  0.006005050
      [4,] -0.179132795 -0.186621759 -0.131767469 -0.145938021 -0.148139879
      [5,] -0.070698483 -0.163375366 -0.143500593 -0.158297368 -0.153147969
      [6,] -0.024298031 -0.160178596 -0.133286039 -0.144788246 -0.056053417
                  [,16]         [,17]         [,18]         [,19]         [,20]
      [1,]  0.001663540  0.0007586967 -0.0018483821 -9.425697e-05 -9.425697e-05
      [2,] -0.012351173 -0.0099150435 -0.0127873328 -1.080785e-02 -1.080785e-02
      [3,]  0.004519487  0.0051129851 -0.0001734336  3.864951e-03  3.864951e-03
      [4,] -0.176292885  0.0257029014 -0.0056338404  1.298343e-02  1.298343e-02
      [5,] -0.189265482 -0.1053083511 -0.1365073302 -1.170383e-01 -1.170383e-01
      [6,] -0.160577611 -0.1881816356 -0.2028752259 -1.895493e-01 -1.895493e-01
                  [,21]        [,22]        [,23]        [,24]        [,25]
      [1,] -0.002062814 -0.002062814  0.001081961  0.001081961  0.004930888
      [2,] -0.012778583 -0.012778583 -0.009595217 -0.009595217 -0.009619263
      [3,]  0.002163017  0.002163017  0.005239294  0.005239294  0.006446378
      [4,] -0.017078256 -0.017078256  0.003692584  0.003692584 -0.083978311
      [5,] -0.147133175 -0.147133175 -0.126107445 -0.126107445 -0.118634088
      [6,] -0.215539373 -0.215539373 -0.195220957 -0.195220957 -0.139960013
                  [,26]        [,27]         [,28]        [,29]        [,30]
      [1,]  0.007958548  0.006733075  9.232744e-03  0.003069195  0.003310677
      [2,] -0.012341188 -0.010151020 -1.101619e-02 -0.012701341 -0.012709521
      [3,] -0.001934612  0.004933680 -8.505483e-05  0.003079866  0.002929368
      [4,] -0.079977878 -0.079994680 -6.997441e-02 -0.101192255 -0.104105640
      [5,] -0.143694020 -0.120277812 -1.332917e-01 -0.138633038 -0.138302578
      [6,] -0.154826151 -0.145053874 -1.403056e-01 -0.161405995 -0.159174350
                  [,31]         [,32]        [,33]        [,34]         [,35]
      [1,]  0.007200325  0.0093251297  0.004720619  0.004720619  0.0012557049
      [2,] -0.009366025 -0.0109431121 -0.012754825 -0.012754825 -0.0149060797
      [3,]  0.005824471  0.0002960907  0.001909557  0.001909557 -0.0003260228
      [4,] -0.090637183 -0.0742255268  0.012645928  0.012645928 -0.0467923780
      [5,] -0.130111278 -0.1378090983 -0.068035697 -0.068035697 -0.1145342558
      [6,] -0.154393508 -0.1454364519 -0.111601265 -0.111601265 -0.1486770065
                   [,36]         [,37]       [,38]        [,39]        [,40]
      [1,]  0.0012557049  3.222807e-03  0.00489317  0.006737791  0.006326841
      [2,] -0.0149060797 -1.412074e-02 -0.01441868 -0.011860549 -0.012197368
      [3,] -0.0003260228  7.483059e-06 -0.00025493  0.003592208  0.001505638
      [4,] -0.0467923780 -1.199745e-01 -0.02753768 -0.013840469 -0.003389256
      [5,] -0.1145342558 -1.643901e-01 -0.11858171 -0.135057082 -0.102111120
      [6,] -0.1486770065 -1.875231e-01 -0.13791887 -0.182140168 -0.127815820
                  [,41]        [,42]        [,43]         [,44]        [,45]
      [1,]  0.005074408  0.006004894  0.002423600  0.0016221496  0.003170901
      [2,] -0.008998245 -0.009858670 -0.011719170 -0.0159526623 -0.014469657
      [3,]  0.006110133  0.004848736  0.004143152 -0.0004264964  0.000512633
      [4,] -0.074738870 -0.076982358 -0.121391818 -0.0432468242 -0.119770936
      [5,] -0.088122464 -0.115138703 -0.136009434 -0.0886935978 -0.151409247
      [6,] -0.141934394 -0.163337681 -0.176552322 -0.1499705620 -0.197728345
                  [,46]        [,47]        [,48]
      [1,]  0.004289740  0.001890719  0.002797694
      [2,] -0.016653157 -0.012215001 -0.013164409
      [3,]  0.002182878  0.004289260  0.003071763
      [4,] -0.027440621 -0.171159229 -0.158054074
      [5,] -0.082521807 -0.185544263 -0.198701324
      [6,] -0.195612179 -0.197080887 -0.208256076

---

    Code
      sapply(mods, function(x) head(fitted(x), 3))
    Output
                 [,1]       [,2]       [,3]       [,4]       [,5]       [,6]
      [1,]  1.0009087  1.0009529  1.0000471  1.0004758  0.9969761  0.9970212
      [2,]  0.9972623  0.9973130  0.9963479  0.9967919  0.9985629  0.9986197
      [3,]  0.9983657  0.9984162  0.9974515  0.9978960  1.0033010  1.0033557
      [4,] 28.2826147 28.2840590 28.2517845 28.2669037 28.4223844 28.4217807
      [5,] 28.4728543 28.4745110 28.4401214 28.4557830 28.7318554 28.7310948
      [6,] 28.7426064 28.7442562 28.7098855 28.7255613 28.8866400 28.8859072
                 [,7]       [,8]        [,9]        [,10]         [,11]         [,12]
      [1,]  0.9956678  0.9963461 0.001299434 0.0003547858 -0.0005321643 -6.806549e-05
      [2,]  0.9973837  0.9979453 0.001778774 0.0008834756  0.0004712847  6.101423e-04
      [3,]  1.0020097  1.0026477 0.002517016 0.0015856339  0.0007638402  1.184396e-03
      [4,] 28.4014188 28.4092741 0.249748490 0.2431933557  0.1882493115  2.257175e-01
      [5,] 28.7124980 28.7190017 0.288963116 0.2827504322  0.2572163725  2.684268e-01
      [6,] 28.8656556 28.8730441 0.253753412 0.2472903334  0.1963825460  2.303353e-01
                 [,13]       [,14]        [,15]        [,16]      [,17]      [,18]
      [1,] 0.002185670 0.001239393 5.170443e-04 0.0008331564  0.9970242  0.9970242
      [2,] 0.002462168 0.001534433 8.704740e-04 0.0011504379  0.9936375  0.9936375
      [3,] 0.002012278 0.001014396 1.453227e-05 0.0005091261  0.9947732  0.9947732
      [4,] 0.290555090 0.284831405 2.466869e-01 0.2689285113 28.2297424 28.2297424
      [5,] 0.253791351 0.248179817 2.131186e-01 0.2328168951 28.4158975 28.4158975
      [6,] 0.243677558 0.237641732 1.848427e-01 0.2196422856 28.6857563 28.6857563
                [,19]      [,20]      [,21]      [,22]      [,23]      [,24]
      [1,]  1.0013439  1.0008249  1.0009691  1.0010379  1.0007890  0.9991615
      [2,]  0.9976688  0.9970114  0.9972612  0.9973552  0.9970795  0.9954526
      [3,]  0.9987658  0.9981479  0.9983616  0.9984765  0.9981816  0.9965968
      [4,] 28.2219551 28.2200569 28.2712871 28.2219934 28.2372088 28.2152714
      [5,] 28.4064061 28.4072995 28.4600601 28.4104247 28.4256855 28.4036086
      [6,] 28.6770503 28.6771562 28.7302599 28.6801430 28.6962180 28.6734092
                [,25]      [,26]      [,27]      [,28]      [,29]      [,30]
      [1,]  0.9925979  0.9925979  0.9969060  0.9969060  0.9966183  0.9960624
      [2,]  0.9950826  0.9950826  0.9989583  0.9989583  0.9983649  0.9984278
      [3,]  0.9989943  0.9989943  1.0031652  1.0031652  1.0028318  1.0027256
      [4,] 28.4084723 28.4084723 28.3887700 28.3887700 28.4104218 28.4056437
      [5,] 28.6919170 28.6919170 28.6908379 28.6908379 28.7225109 28.6920062
      [6,] 28.8647534 28.8647534 28.8520763 28.8520763 28.8749818 28.8636046
                [,31]      [,32]       [,33]         [,34]       [,35]        [,36]
      [1,]  0.9967462  0.9970475 0.002964800 -1.566191e-03 0.001580549 0.0005643716
      [2,]  0.9985045  0.9990152 0.003750998 -1.208250e-03 0.002583424 0.0013109001
      [3,]  1.0029580  1.0033396 0.003911315 -7.532505e-05 0.002831116 0.0013766080
      [4,] 28.3860113 28.3770753 0.227004175  1.895208e-01 0.233444707 0.2170176866
      [5,] 28.6958640 28.6918869 0.277730460  2.323134e-01 0.283849314 0.2682618171
      [6,] 28.8508871 28.8446586 0.228725989  1.811180e-01 0.235031855 0.2116654940
                 [,37]       [,38]        [,39]       [,40]       [,41]        [,42]
      [1,] 0.001086156 0.001159545 0.0009881709 0.002176092 0.003576772 0.0016827874
      [2,] 0.001966417 0.002018153 0.0016527252 0.002629021 0.001628267 0.0003227808
      [3,] 0.002373531 0.002453100 0.0019613160 0.002779665 0.001358029 0.0001569680
      [4,] 0.246133404 0.247571177 0.2152124135 0.207119599 0.221705408 0.2184696452
      [5,] 0.304177047 0.305190615 0.2622191478 0.260838222 0.215544135 0.2095319264
      [6,] 0.243488573 0.245047428 0.2150355286 0.203128661 0.194313851 0.1929458895
                 [,43]        [,44]       [,45]       [,46]       [,47]        [,48]
      [1,] 0.003393560 0.0018567333 0.003959123 0.003994590 0.003880298 0.0019739327
      [2,] 0.001406785 0.0005197752 0.002626233 0.002274246 0.001956690 0.0006134883
      [3,] 0.001169772 0.0003031857 0.002564813 0.002243001 0.001585628 0.0003786570
      [4,] 0.227653297 0.2336251449 0.291233410 0.288245598 0.220435944 0.2150670208
      [5,] 0.222734419 0.2266955940 0.243799106 0.244511219 0.214170541 0.2061344179
      [6,] 0.200425543 0.2056854926 0.223476630 0.221801549 0.193361946 0.1903550036

---

    Code
      sapply(mods, function(x) tail(fitted(x), 3))
    Output
                 [,1]       [,2]       [,3]       [,4]       [,5]       [,6]
      [1,]   1.477556   1.477867   1.479031   1.477182   1.473192   1.473598
      [2,]   1.478057   1.478374   1.479533   1.477667   1.477701   1.478077
      [3,]   1.467977   1.468310   1.469330   1.467496   1.466378   1.466754
      [4,] 110.376167 110.386333 110.427790 110.362574 110.522308 110.516876
      [5,] 110.575843 110.586188 110.627510 110.561697 110.585430 110.580397
      [6,] 110.646716 110.657602 110.693941 110.629254 110.590242 110.585219
                 [,7]       [,8]          [,9]         [,10]        [,11]
      [1,]   1.475689   1.472879 -0.0041359353 -0.0017103905 -0.001261170
      [2,]   1.479477   1.477209 -0.0001198563  0.0019032456  0.000754132
      [3,]   1.467810   1.465760 -0.0012491574  0.0005996997 -0.001352866
      [4,] 110.548336 110.515797  0.3344734591  0.3513048875  0.379132795
      [5,] 110.601456 110.575195  0.1978440251  0.2118828075  0.140698483
      [6,] 110.601109 110.577367  0.1324242575  0.1452539143  0.024298031
                  [,12]         [,13]         [,14]        [,15]         [,16]
      [1,] -0.001168408 -0.0039481332 -0.0016053652 -0.001563668 -0.0011635402
      [2,]  0.001902074 -0.0007109355  0.0017353637  0.001637849  0.0021511734
      [3,]  0.000330218 -0.0014262855  0.0004753343 -0.001205050  0.0002805132
      [4,]  0.386621759  0.3317674689  0.3459380207  0.348139879  0.3762928845
      [5,]  0.233375366  0.2135005935  0.2282973675  0.223147969  0.2592654816
      [6,]  0.160178596  0.1332860394  0.1447882464  0.056053417  0.1605776114
                [,17]      [,18]      [,19]      [,20]      [,21]      [,22]
      [1,]   1.476841   1.479448   1.477694   1.477694   1.479663   1.479663
      [2,]   1.477315   1.480187   1.478208   1.478208   1.480179   1.480179
      [3,]   1.467087   1.472373   1.468335   1.468335   1.470037   1.470037
      [4,] 110.354297 110.385634 110.367017 110.367017 110.397078 110.397078
      [5,] 110.555308 110.586507 110.567038 110.567038 110.597133 110.597133
      [6,] 110.638182 110.652875 110.639549 110.639549 110.665539 110.665539
                [,23]      [,24]      [,25]      [,26]      [,27]      [,28]
      [1,]   1.476518   1.476518   1.472669   1.469641   1.470867   1.468367
      [2,]   1.476995   1.476995   1.477019   1.479741   1.477551   1.478416
      [3,]   1.466961   1.466961   1.465754   1.474135   1.467266   1.472285
      [4,] 110.376307 110.376307 110.463978 110.459978 110.459995 110.449974
      [5,] 110.576107 110.576107 110.568634 110.593694 110.570278 110.583292
      [6,] 110.645221 110.645221 110.589960 110.604826 110.595054 110.590306
                [,29]      [,30]      [,31]      [,32]        [,33]        [,34]
      [1,]   1.474531   1.474289   1.470400   1.468275 -0.004220619 -0.004220619
      [2,]   1.480101   1.480110   1.476766   1.478343  0.002554825  0.002554825
      [3,]   1.469120   1.469271   1.466376   1.471904  0.002890443  0.002890443
      [4,] 110.481192 110.484106 110.470637 110.454226  0.187354072  0.187354072
      [5,] 110.588633 110.588303 110.580111 110.587809  0.138035697  0.138035697
      [6,] 110.611406 110.609174 110.604394 110.595436  0.111601265  0.111601265
                   [,35]         [,36]        [,37]        [,38]        [,39]
      [1,] -0.0007557049 -0.0007557049 -0.002722807 -0.004393170 -0.006237791
      [2,]  0.0047060797  0.0047060797  0.003920743  0.004218684  0.001660549
      [3,]  0.0051260228  0.0051260228  0.004792517  0.005054930  0.001207792
      [4,]  0.2467923780  0.2467923780  0.319974470  0.227537683  0.213840469
      [5,]  0.1845342558  0.1845342558  0.234390088  0.188581705  0.205057082
      [6,]  0.1486770065  0.1486770065  0.187523106  0.137918871  0.182140168
                  [,40]        [,41]         [,42]         [,43]        [,44]
      [1,] -0.005826841 -0.004574408 -5.504894e-03 -0.0019235999 -0.001122150
      [2,]  0.001997368 -0.001201755 -3.413303e-04  0.0015191704  0.005752662
      [3,]  0.003294362 -0.001310133 -4.873627e-05  0.0006568483  0.005226496
      [4,]  0.203389256  0.274738870  2.769824e-01  0.3213918182  0.243246824
      [5,]  0.172111120  0.158122464  1.851387e-01  0.2060094340  0.158693598
      [6,]  0.127815820  0.141934394  1.633377e-01  0.1765523223  0.149970562
                  [,45]        [,46]         [,47]        [,48]
      [1,] -0.002670901 -0.003789740 -0.0013907190 -0.002297694
      [2,]  0.004269657  0.006453157  0.0020150009  0.002964409
      [3,]  0.004287367  0.002617122  0.0005107403  0.001728237
      [4,]  0.319770936  0.227440621  0.3711592295  0.358054074
      [5,]  0.221409247  0.152521807  0.2555442634  0.268701324
      [6,]  0.197728345  0.195612179  0.1970808866  0.208256076

# mod_refit_check

    Code
      suppressMessages(suppressWarnings(sapply(mods, tsDyn:::mod_refit_check)))
    Output
       [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
      [16] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
      [31] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
      [46] TRUE TRUE TRUE

# Linear models

    Code
      sapply(mods_Linear, function(x) round(confint(x), 3))
    Output
      [[1]]
                        2.5 % 97.5 %
      dolcan:Intercept -0.014  0.035
      dolcan:Trend      0.379  0.379
      dolcan:dolcan-1  -0.018  0.018
      dolcan:cpiUSA-1  -0.001  0.000
      cpiUSA:Intercept  0.652  1.322
      cpiUSA:Trend     -0.142 -0.137
      cpiUSA:dolcan-1  -0.239  0.239
      cpiUSA:cpiUSA-1   0.993  1.011
      
      [[2]]
                        2.5 % 97.5 %
      dolcan:Intercept -0.003  0.029
      dolcan:dolcan-1   0.438  0.473
      dolcan:cpiUSA-1   0.987  0.987
      cpiUSA:Intercept -0.361  0.079
      cpiUSA:dolcan-1  -0.239  0.239
      cpiUSA:cpiUSA-1   0.998  1.001
      
      [[3]]
                       2.5 % 97.5 %
      dolcan:dolcan-1  0.996  1.005
      dolcan:cpiUSA-1  0.335  0.335
      cpiUSA:dolcan-1 -0.066  0.066
      cpiUSA:cpiUSA-1  0.997  0.999
      
      [[4]]
                       2.5 % 97.5 %
      dolcan:Trend     0.000  0.000
      dolcan:dolcan-1 -0.016  0.011
      dolcan:cpiUSA-1  0.991  0.992
      cpiUSA:Trend     0.035  0.038
      cpiUSA:dolcan-1 -0.183  0.183
      cpiUSA:cpiUSA-1  1.002  1.016
      
      [[5]]
                        2.5 % 97.5 %
      dolcan:Intercept -0.011  0.038
      dolcan:Trend      0.283  0.283
      dolcan:dolcan-1  -0.109  0.109
      dolcan:cpiUSA-1  -0.008  0.008
      dolcan:dolcan-2   1.041  1.259
      dolcan:cpiUSA-2   0.525  0.541
      cpiUSA:Intercept -0.306  0.299
      cpiUSA:Trend      1.455  1.460
      cpiUSA:dolcan-1  -1.511  1.179
      cpiUSA:cpiUSA-1  -0.704 -0.506
      cpiUSA:dolcan-2  -1.340  1.348
      cpiUSA:cpiUSA-2  -0.558 -0.360
      
      [[6]]
                        2.5 % 97.5 %
      dolcan:Intercept  0.000  0.033
      dolcan:dolcan-1   0.137  0.355
      dolcan:cpiUSA-1   1.143  1.159
      dolcan:dolcan-2   0.421  0.638
      dolcan:cpiUSA-2  -0.012  0.004
      cpiUSA:Intercept  1.254  1.660
      cpiUSA:dolcan-1  -1.510  1.177
      cpiUSA:cpiUSA-1  -0.699 -0.502
      cpiUSA:dolcan-2  -1.338  1.346
      cpiUSA:cpiUSA-2  -0.556 -0.359
      
      [[7]]
                       2.5 % 97.5 %
      dolcan:dolcan-1  1.052  1.270
      dolcan:cpiUSA-1  0.676  0.692
      dolcan:dolcan-2 -0.111  0.107
      dolcan:cpiUSA-2  1.476  1.491
      cpiUSA:dolcan-1 -1.507  1.187
      cpiUSA:cpiUSA-1 -0.604 -0.411
      cpiUSA:dolcan-2 -1.348  1.352
      cpiUSA:cpiUSA-2 -0.581 -0.388
      
      [[8]]
                       2.5 % 97.5 %
      dolcan:Trend     0.000  0.000
      dolcan:dolcan-1 -0.110  0.108
      dolcan:cpiUSA-1  1.146  1.161
      dolcan:dolcan-2  0.489  0.707
      dolcan:cpiUSA-2 -0.011  0.005
      cpiUSA:Trend     1.468  1.471
      cpiUSA:dolcan-1 -1.511  1.186
      cpiUSA:cpiUSA-1 -0.637 -0.440
      cpiUSA:dolcan-2 -1.344  1.351
      cpiUSA:cpiUSA-2 -0.565 -0.367
      
      [[9]]
                        2.5 % 97.5 %
      dolcan:ECT       -0.032  0.002
      dolcan:Intercept -0.061 -0.040
      dolcan:Trend      0.011  0.011
      dolcan:dolcan-1   0.084  0.302
      dolcan:cpiUSA-1  -0.008  0.008
      cpiUSA:ECT       -0.205  0.204
      cpiUSA:Intercept  0.038  0.293
      cpiUSA:Trend      0.587  0.588
      cpiUSA:dolcan-1  -1.346  1.338
      cpiUSA:cpiUSA-1   0.358  0.556
      
      [[10]]
                        2.5 % 97.5 %
      dolcan:ECT       -0.006  0.004
      dolcan:Intercept  0.044  0.049
      dolcan:dolcan-1  -0.107  0.111
      dolcan:cpiUSA-1   0.124  0.140
      cpiUSA:ECT        0.099  0.220
      cpiUSA:Intercept  0.515  0.576
      cpiUSA:dolcan-1  -1.342  1.336
      cpiUSA:cpiUSA-1   0.366  0.562
      
      [[11]]
                       2.5 % 97.5 %
      dolcan:ECT      -0.006  0.004
      dolcan:dolcan-1 -0.036  0.182
      dolcan:cpiUSA-1  0.163  0.173
      cpiUSA:ECT       1.010  1.144
      cpiUSA:dolcan-1 -1.471  1.476
      cpiUSA:cpiUSA-1  0.736  0.865
      
      [[12]]
                       2.5 % 97.5 %
      dolcan:ECT      -0.005  0.007
      dolcan:Trend     0.239  0.239
      dolcan:dolcan-1 -0.109  0.109
      dolcan:cpiUSA-1 -0.007  0.008
      cpiUSA:ECT       0.087  0.235
      cpiUSA:Trend     0.508  0.508
      cpiUSA:dolcan-1 -1.359  1.356
      cpiUSA:cpiUSA-1  0.410  0.600
      
      [[13]]
                        2.5 % 97.5 %
      dolcan:ECT       -0.032  0.002
      dolcan:Intercept -0.050 -0.029
      dolcan:Trend      0.012  0.012
      dolcan:dolcan-1   0.074  0.295
      dolcan:cpiUSA-1  -0.009  0.009
      dolcan:dolcan-2  -0.111  0.111
      dolcan:cpiUSA-2   0.159  0.177
      cpiUSA:ECT        0.485  0.899
      cpiUSA:Intercept -0.134  0.129
      cpiUSA:Trend      0.453  0.454
      cpiUSA:dolcan-1  -1.386  1.342
      cpiUSA:cpiUSA-1  -0.741 -0.519
      cpiUSA:dolcan-2  -1.373  1.367
      cpiUSA:cpiUSA-2  -0.100  0.122
      
      [[14]]
                        2.5 % 97.5 %
      dolcan:ECT       -0.006  0.004
      dolcan:Intercept  0.043  0.049
      dolcan:dolcan-1  -0.108  0.113
      dolcan:cpiUSA-1   0.121  0.139
      dolcan:dolcan-2   0.053  0.274
      dolcan:cpiUSA-2   0.657  0.675
      cpiUSA:ECT       -0.063  0.060
      cpiUSA:Intercept  0.423  0.492
      cpiUSA:dolcan-1  -1.393  1.332
      cpiUSA:cpiUSA-1  -0.792 -0.570
      cpiUSA:dolcan-2  -1.366  1.362
      cpiUSA:cpiUSA-2  -0.095  0.127
      
      [[15]]
                       2.5 % 97.5 %
      dolcan:ECT      -0.006  0.004
      dolcan:dolcan-1 -0.050  0.172
      dolcan:cpiUSA-1  0.164  0.180
      dolcan:dolcan-2  1.018  1.240
      dolcan:cpiUSA-2 -0.007  0.010
      cpiUSA:ECT       0.572  0.705
      cpiUSA:dolcan-1 -1.494  1.443
      cpiUSA:cpiUSA-1 -0.517 -0.300
      cpiUSA:dolcan-2 -1.473  1.476
      cpiUSA:cpiUSA-2  0.094  0.309
      
      [[16]]
                       2.5 % 97.5 %
      dolcan:ECT      -0.005  0.008
      dolcan:Trend     0.230  0.230
      dolcan:dolcan-1 -0.111  0.111
      dolcan:cpiUSA-1 -0.008  0.009
      dolcan:dolcan-2  0.054  0.277
      dolcan:cpiUSA-2  0.646  0.664
      cpiUSA:ECT      -0.079  0.078
      cpiUSA:Trend     0.484  0.484
      cpiUSA:dolcan-1 -1.409  1.348
      cpiUSA:cpiUSA-1 -0.872 -0.652
      cpiUSA:dolcan-2 -1.382  1.380
      cpiUSA:cpiUSA-2 -0.066  0.154
      

---

    Code
      map_dfr(mods_Linear, function(x) mutate(tidy(x), across(where(is.numeric),
      ~ round(., 3))))
    Output
          equation      term estimate std.error statistic p.value
      1     dolcan Intercept    0.010     0.012     0.838   0.403
      2     dolcan     Trend    0.000     0.000    -0.249   0.803
      3     dolcan dolcan -1    0.987     0.009   110.837   0.000
      4     dolcan cpiUSA -1    0.000     0.000     0.442   0.659
      5     cpiUSA Intercept    0.379     0.170     2.226   0.027
      6     cpiUSA     Trend   -0.001     0.001    -0.596   0.551
      7     cpiUSA dolcan -1   -0.140     0.122    -1.147   0.252
      8     cpiUSA cpiUSA -1    1.002     0.004   222.856   0.000
      9     dolcan Intercept    0.013     0.008     1.559   0.120
      10    dolcan dolcan -1    0.987     0.009   111.020   0.000
      11    dolcan cpiUSA -1    0.000     0.000     1.252   0.212
      12    cpiUSA Intercept    0.455     0.112     4.065   0.000
      13    cpiUSA dolcan -1   -0.141     0.121    -1.162   0.246
      14    cpiUSA cpiUSA -1    1.000     0.001  1421.261   0.000
      15    dolcan dolcan -1    1.000     0.002   419.096   0.000
      16    dolcan cpiUSA -1    0.000     0.000     0.307   0.759
      17    cpiUSA dolcan -1    0.335     0.033    10.043   0.000
      18    cpiUSA cpiUSA -1    0.998     0.001  1830.173   0.000
      19    dolcan     Trend    0.000     0.000    -1.336   0.183
      20    dolcan dolcan -1    0.992     0.007   146.709   0.000
      21    dolcan cpiUSA -1    0.000     0.000     1.368   0.172
      22    cpiUSA     Trend   -0.003     0.001    -3.427   0.001
      23    cpiUSA dolcan -1    0.037     0.093     0.393   0.694
      24    cpiUSA cpiUSA -1    1.009     0.003   302.527   0.000
      25    dolcan Intercept    0.014     0.012     1.092   0.275
      26    dolcan     Trend    0.000     0.000    -0.299   0.765
      27    dolcan dolcan -1    1.150     0.055    20.756   0.000
      28    dolcan cpiUSA -1   -0.004     0.004    -0.920   0.358
      29    dolcan dolcan -2   -0.166     0.055    -2.999   0.003
      30    dolcan cpiUSA -2    0.004     0.004     0.957   0.339
      31    cpiUSA Intercept    0.283     0.154     1.843   0.066
      32    cpiUSA     Trend    0.000     0.001     0.325   0.746
      33    cpiUSA dolcan -1    0.533     0.684     0.780   0.436
      34    cpiUSA cpiUSA -1    1.458     0.050    29.055   0.000
      35    cpiUSA dolcan -2   -0.605     0.683    -0.885   0.377
      36    cpiUSA cpiUSA -2   -0.459     0.050    -9.106   0.000
      37    dolcan Intercept    0.016     0.008     1.959   0.051
      38    dolcan dolcan -1    1.151     0.055    20.794   0.000
      39    dolcan cpiUSA -1   -0.004     0.004    -0.914   0.362
      40    dolcan dolcan -2   -0.166     0.055    -3.009   0.003
      41    dolcan cpiUSA -2    0.004     0.004     0.932   0.352
      42    cpiUSA Intercept    0.246     0.103     2.388   0.018
      43    cpiUSA dolcan -1    0.529     0.683     0.776   0.439
      44    cpiUSA cpiUSA -1    1.457     0.050    29.098   0.000
      45    cpiUSA dolcan -2   -0.601     0.682    -0.880   0.379
      46    cpiUSA cpiUSA -2   -0.457     0.050    -9.136   0.000
      47    dolcan dolcan -1    1.161     0.055    20.982   0.000
      48    dolcan cpiUSA -1   -0.002     0.004    -0.490   0.625
      49    dolcan dolcan -2   -0.160     0.055    -2.889   0.004
      50    dolcan cpiUSA -2    0.002     0.004     0.493   0.622
      51    cpiUSA dolcan -1    0.684     0.685     1.000   0.318
      52    cpiUSA cpiUSA -1    1.484     0.049    30.161   0.000
      53    cpiUSA dolcan -2   -0.507     0.686    -0.739   0.460
      54    cpiUSA cpiUSA -2   -0.485     0.049    -9.877   0.000
      55    dolcan     Trend    0.000     0.000    -1.649   0.100
      56    dolcan dolcan -1    1.153     0.055    20.833   0.000
      57    dolcan cpiUSA -1   -0.003     0.004    -0.786   0.432
      58    dolcan dolcan -2   -0.163     0.055    -2.944   0.003
      59    dolcan cpiUSA -2    0.004     0.004     0.878   0.381
      60    cpiUSA     Trend   -0.001     0.001    -1.543   0.124
      61    cpiUSA dolcan -1    0.598     0.685     0.872   0.384
      62    cpiUSA cpiUSA -1    1.469     0.050    29.426   0.000
      63    cpiUSA dolcan -2   -0.538     0.685    -0.786   0.432
      64    cpiUSA cpiUSA -2   -0.466     0.050    -9.234   0.000
      65    dolcan       ECT   -0.015     0.008    -1.787   0.075
      66    dolcan Intercept    0.011     0.005     2.088   0.038
      67    dolcan     Trend    0.000     0.000    -1.735   0.084
      68    dolcan dolcan -1    0.166     0.055     2.996   0.003
      69    dolcan cpiUSA -1   -0.004     0.004    -0.977   0.329
      70    cpiUSA       ECT   -0.050     0.104    -0.484   0.629
      71    cpiUSA Intercept    0.193     0.065     2.972   0.003
      72    cpiUSA     Trend    0.000     0.000    -0.975   0.330
      73    cpiUSA dolcan -1    0.587     0.682     0.861   0.390
      74    cpiUSA cpiUSA -1    0.457     0.050     9.091   0.000
      75    dolcan       ECT   -0.001     0.003    -0.439   0.661
      76    dolcan Intercept    0.002     0.001     1.678   0.094
      77    dolcan dolcan -1    0.160     0.055     2.884   0.004
      78    dolcan cpiUSA -1   -0.003     0.004    -0.743   0.458
      79    cpiUSA       ECT    0.047     0.031     1.503   0.134
      80    cpiUSA Intercept    0.132     0.016     8.446   0.000
      81    cpiUSA dolcan -1    0.545     0.681     0.801   0.424
      82    cpiUSA cpiUSA -1    0.464     0.050     9.314   0.000
      83    dolcan       ECT   -0.001     0.003    -0.269   0.788
      84    dolcan dolcan -1    0.168     0.055     3.043   0.003
      85    dolcan cpiUSA -1    0.002     0.002     1.000   0.318
      86    cpiUSA       ECT    0.073     0.034     2.149   0.032
      87    cpiUSA dolcan -1    1.077     0.749     1.438   0.151
      88    cpiUSA cpiUSA -1    0.800     0.033    24.314   0.000
      89    dolcan       ECT    0.001     0.003     0.453   0.651
      90    dolcan     Trend    0.000     0.000     1.212   0.227
      91    dolcan dolcan -1    0.161     0.056     2.902   0.004
      92    dolcan cpiUSA -1   -0.001     0.004    -0.317   0.752
      93    cpiUSA       ECT    0.239     0.038     6.342   0.000
      94    cpiUSA     Trend    0.001     0.000     7.869   0.000
      95    cpiUSA dolcan -1    0.508     0.690     0.737   0.462
      96    cpiUSA cpiUSA -1    0.505     0.048    10.495   0.000
      97    dolcan       ECT   -0.015     0.009    -1.788   0.075
      98    dolcan Intercept    0.012     0.005     2.140   0.033
      99    dolcan     Trend    0.000     0.000    -1.746   0.082
      100   dolcan dolcan -1    0.168     0.056     2.987   0.003
      101   dolcan cpiUSA -1   -0.003     0.005    -0.555   0.579
      102   dolcan dolcan -2   -0.022     0.056    -0.389   0.697
      103   dolcan cpiUSA -2   -0.003     0.005    -0.673   0.502
      104   cpiUSA       ECT   -0.040     0.105    -0.379   0.705
      105   cpiUSA Intercept    0.185     0.067     2.773   0.006
      106   cpiUSA     Trend    0.000     0.000    -0.856   0.393
      107   cpiUSA dolcan -1    0.692     0.693     0.998   0.319
      108   cpiUSA cpiUSA -1    0.454     0.057     8.031   0.000
      109   cpiUSA dolcan -2   -0.630     0.696    -0.905   0.366
      110   cpiUSA cpiUSA -2    0.011     0.057     0.194   0.846
      111   dolcan       ECT   -0.001     0.003    -0.406   0.685
      112   dolcan Intercept    0.002     0.001     1.731   0.084
      113   dolcan dolcan -1    0.164     0.056     2.903   0.004
      114   dolcan cpiUSA -1   -0.002     0.005    -0.413   0.680
      115   dolcan dolcan -2   -0.030     0.056    -0.538   0.591
      116   dolcan cpiUSA -2   -0.002     0.005    -0.501   0.616
      117   cpiUSA       ECT    0.046     0.031     1.476   0.141
      118   cpiUSA Intercept    0.130     0.017     7.435   0.000
      119   cpiUSA dolcan -1    0.666     0.692     0.961   0.337
      120   cpiUSA cpiUSA -1    0.458     0.056     8.131   0.000
      121   cpiUSA dolcan -2   -0.681     0.693    -0.982   0.327
      122   cpiUSA cpiUSA -2    0.016     0.056     0.280   0.780
      123   dolcan       ECT   -0.001     0.003    -0.296   0.768
      124   dolcan dolcan -1    0.172     0.056     3.062   0.002
      125   dolcan cpiUSA -1    0.002     0.004     0.369   0.712
      126   dolcan dolcan -2   -0.025     0.057    -0.446   0.656
      127   dolcan cpiUSA -2    0.001     0.004     0.297   0.767
      128   cpiUSA       ECT    0.061     0.034     1.803   0.072
      129   cpiUSA dolcan -1    1.129     0.746     1.512   0.131
      130   cpiUSA cpiUSA -1    0.638     0.055    11.614   0.000
      131   cpiUSA dolcan -2   -0.409     0.749    -0.545   0.586
      132   cpiUSA cpiUSA -2    0.202     0.055     3.689   0.000
      133   dolcan       ECT    0.002     0.003     0.513   0.608
      134   dolcan     Trend    0.000     0.000     1.211   0.227
      135   dolcan dolcan -1    0.166     0.057     2.930   0.004
      136   dolcan cpiUSA -1   -0.001     0.005    -0.148   0.883
      137   dolcan dolcan -2   -0.030     0.057    -0.534   0.594
      138   dolcan cpiUSA -2   -0.001     0.005    -0.225   0.822
      139   cpiUSA       ECT    0.230     0.040     5.757   0.000
      140   cpiUSA     Trend    0.001     0.000     6.875   0.000
      141   cpiUSA dolcan -1    0.655     0.701     0.935   0.350
      142   cpiUSA cpiUSA -1    0.484     0.056     8.631   0.000
      143   cpiUSA dolcan -2   -0.762     0.702    -1.086   0.278
      144   cpiUSA cpiUSA -2    0.044     0.056     0.788   0.431

# Non linear functions

    Code
      sapply(mods_nonLIn, function(x) head(regime(x), 3))
    Output
           [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]
      [1,]   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA    NA    NA    NA    NA
      [2,]    1    1    1    1    1    1    1    1   NA    NA    NA    NA    NA    NA
      [3,]    1    1    1    1    1    1    1    1    1     1     1     1     1     1
           [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24] [,25] [,26]
      [1,]    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
      [2,]    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
      [3,]     1     1     2     3     2     3     2     3     2     3    NA    NA
           [,27] [,28] [,29] [,30] [,31] [,32]
      [1,]    NA    NA    NA    NA    NA    NA
      [2,]    NA    NA    NA    NA    NA    NA
      [3,]    NA    NA    NA    NA    NA    NA

---

    Code
      sapply(mods_nonLIn, function(x) tail(regime(x), 3))
    Output
           [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]
      [1,]    2    3    2    3    2    3    2    3    2     3     2     3     2     3
      [2,]    2    3    2    3    2    3    2    3    2     3     2     3     2     3
      [3,]    2    3    2    3    2    3    2    3    2     3     2     3     2     3
           [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24] [,25] [,26]
      [1,]     2     3     1     1     1     1     1     1     1     1     1     1
      [2,]     2     3     1     1     1     1     1     1     1     1     1     1
      [3,]     2     3     1     1     1     1     1     1     1     1     1     1
           [,27] [,28] [,29] [,30] [,31] [,32]
      [1,]     1     1     1     1     1     1
      [2,]     1     1     1     1     1     1
      [3,]     1     1     1     1     1     1

---

    Code
      sapply(mods_nonLIn, function(x) head(regime(x, initVal = FALSE), 3))
    Output
           [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]
      [1,]    1    1    1    1    1    1    1    1    1     1     1     1     1     1
      [2,]    1    1    1    1    1    1    1    1    1     1     1     1     1     1
      [3,]    1    1    1    1    1    1    1    1    1     1     1     1     1     1
           [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24] [,25] [,26]
      [1,]     1     1     2     3     2     3     2     3     2     3     2     3
      [2,]     1     1     2     3     2     3     2     3     2     3     2     3
      [3,]     1     1     2     3     2     3     2     3     2     3     2     3
           [,27] [,28] [,29] [,30] [,31] [,32]
      [1,]     2     3     2     3     2     3
      [2,]     2     3     2     3     2     3
      [3,]     2     3     2     3     2     3

---

    Code
      sapply(mods_nonLIn, function(x) tail(regime(x, initVal = FALSE), 3))
    Output
           [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]
      [1,]    2    3    2    3    2    3    2    3    2     3     2     3     2     3
      [2,]    2    3    2    3    2    3    2    3    2     3     2     3     2     3
      [3,]    2    3    2    3    2    3    2    3    2     3     2     3     2     3
           [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24] [,25] [,26]
      [1,]     2     3     1     1     1     1     1     1     1     1     1     1
      [2,]     2     3     1     1     1     1     1     1     1     1     1     1
      [3,]     2     3     1     1     1     1     1     1     1     1     1     1
           [,27] [,28] [,29] [,30] [,31] [,32]
      [1,]     1     1     1     1     1     1
      [2,]     1     1     1     1     1     1
      [3,]     1     1     1     1     1     1

---

    Code
      sapply(mods_nonLIn, function(x) head(regime(x, time = FALSE), 3))
    Output
           [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]
      [1,]   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA    NA    NA    NA    NA
      [2,]    1    1    1    1    1    1    1    1   NA    NA    NA    NA    NA    NA
      [3,]    1    1    1    1    1    1    1    1    1     1     1     1     1     1
           [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24] [,25] [,26]
      [1,]    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
      [2,]    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
      [3,]     1     1     2     3     2     3     2     3     2     3    NA    NA
           [,27] [,28] [,29] [,30] [,31] [,32]
      [1,]    NA    NA    NA    NA    NA    NA
      [2,]    NA    NA    NA    NA    NA    NA
      [3,]    NA    NA    NA    NA    NA    NA

---

    Code
      sapply(mods_nonLIn, function(x) head(regime(x, time = FALSE, initVal = FALSE),
      3))
    Output
           [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]
      [1,]    1    1    1    1    1    1    1    1    1     1     1     1     1     1
      [2,]    1    1    1    1    1    1    1    1    1     1     1     1     1     1
      [3,]    1    1    1    1    1    1    1    1    1     1     1     1     1     1
           [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24] [,25] [,26]
      [1,]     1     1     2     3     2     3     2     3     2     3     2     3
      [2,]     1     1     2     3     2     3     2     3     2     3     2     3
      [3,]     1     1     2     3     2     3     2     3     2     3     2     3
           [,27] [,28] [,29] [,30] [,31] [,32]
      [1,]     2     3     2     3     2     3
      [2,]     2     3     2     3     2     3
      [3,]     2     3     2     3     2     3

# toLatex

    Code
      sapply(mods, toLatex)
    Output
      [[1]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \begin{smatrix}     %const
      0.0104 \\ 0.3790
      \end{smatrix}+\begin{smatrix}     %trend
      -2.2\text{e-05} \\ -0.0007
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.9871 & 0.0001 \\
      -0.1395 & 1.0023 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      \end{equation}
      
      [[2]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \begin{smatrix}     %const
      0.0128 \\ 0.4554
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.9871 & 6.4\text{e-05} \\
      -0.1411 & 0.9996 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      \end{equation}
      
      [[3]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      +\begin{smatrix}      %Lag1
      1.0004 & 1.2\text{e-05} \\
      0.3346 & 0.9978 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      \end{equation}
      
      [[4]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \begin{smatrix}     %trend
      -7.7\text{e-05} \\ -0.0027
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      0.9920 & 0.0003 \\
      0.0365 & 1.0090 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      \end{equation}
      
      [[5]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \begin{smatrix}     %const
      0.0136 \\ 0.2832
      \end{smatrix}+\begin{smatrix}     %trend
      -2.6\text{e-05} \\ 0.0003
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      1.1504 & -0.0037 \\
      0.5332 & 1.4575 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1661 & 0.0039 \\
      -0.6049 & -0.4591 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      \end{equation}
      
      [[6]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \begin{smatrix}     %const
      0.0164 \\ 0.2463
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      1.1507 & -0.0037 \\
      0.5294 & 1.4571 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1664 & 0.0038 \\
      -0.6006 & -0.4574 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      \end{equation}
      
      [[7]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      +\begin{smatrix}      %Lag1
      1.1610 & -0.0019 \\
      0.6843 & 1.4836 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1602 & 0.0020 \\
      -0.5073 & -0.4848 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      \end{equation}
      
      [[8]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \begin{smatrix}     %trend
      -9.6\text{e-05} \\ -0.0011
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      1.1535 & -0.0032 \\
      0.5978 & 1.4694 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1629 & 0.0036 \\
      -0.5384 & -0.4660 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      \end{equation}
      
      [[9]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      +\begin{smatrix}  %ECT
      -0.0151 \\ -0.0503
      \end{smatrix}ECT_{-1}
      \begin{smatrix}     %const
      0.0110 \\ 0.1931
      \end{smatrix}+\begin{smatrix}     %trend
      -4.5\text{e-05} \\ -0.0003
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.1656 & -0.0040 \\
      0.5872 & 0.4569 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      \end{equation}
      
      [[10]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      +\begin{smatrix}  %ECT
      -0.0011 \\ 0.0466
      \end{smatrix}ECT_{-1}
      \begin{smatrix}     %const
      0.0021 \\ 0.1316
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.1595 & -0.0030 \\
      0.5453 & 0.4636 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      \end{equation}
      
      [[11]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      +\begin{smatrix}  %ECT
      -0.0007 \\ 0.0732
      \end{smatrix}ECT_{-1}
      +\begin{smatrix}      %Lag1
      0.1681 & 0.0024 \\
      1.0771 & 0.8004 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      \end{equation}
      
      [[12]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      +\begin{smatrix}  %ECT
      0.0014 \\ 0.2386
      \end{smatrix}ECT_{-1}
      \begin{smatrix}     %trend
      7.6\text{e-06} \\ 0.0006
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      0.1611 & -0.0012 \\
      0.5082 & 0.5052 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      \end{equation}
      
      [[13]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      +\begin{smatrix}  %ECT
      -0.0152 \\ -0.0398
      \end{smatrix}ECT_{-1}
      \begin{smatrix}     %const
      0.0116 \\ 0.1849
      \end{smatrix}+\begin{smatrix}     %trend
      -4.6\text{e-05} \\ -0.0003
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.1679 & -0.0025 \\
      0.6918 & 0.4539 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.0220 & -0.0031 \\
      -0.6301 & 0.0110 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      \end{equation}
      
      [[14]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      +\begin{smatrix}  %ECT
      -0.0010 \\ 0.0461
      \end{smatrix}ECT_{-1}
      \begin{smatrix}     %const
      0.0025 \\ 0.1298
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.1635 & -0.0019 \\
      0.6656 & 0.4578 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.0303 & -0.0023 \\
      -0.6808 & 0.0158 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      \end{equation}
      
      [[15]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      +\begin{smatrix}  %ECT
      -0.0008 \\ 0.0609
      \end{smatrix}ECT_{-1}
      +\begin{smatrix}      %Lag1
      0.1723 & 0.0015 \\
      1.1289 & 0.6385 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.0252 & 0.0012 \\
      -0.4086 & 0.2017 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      \end{equation}
      
      [[16]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      +\begin{smatrix}  %ECT
      0.0017 \\ 0.2304
      \end{smatrix}ECT_{-1}
      \begin{smatrix}     %trend
      8.4\text{e-06} \\ 0.0006
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      0.1656 & -0.0007 \\
      0.6551 & 0.4838 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.0302 & -0.0010 \\
      -0.7619 & 0.0441 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      \end{equation}
      
      [[17]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix}     %const
      0.0821 \\ -0.5559
      \end{smatrix}+\begin{smatrix}     %trend
      0.0007 \\ -0.0018
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.9823 & -0.0024 \\
      0.7181 & 1.0032 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}< 1.1439 \\
      \begin{smatrix}     %const
      -0.0135 \\ 2.2021
      \end{smatrix}+\begin{smatrix}     %trend
      -0.0002 \\ 0.0073
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.9915 & 0.0007 \\
      -0.7468 & 0.9702 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}> 1.1439 \\
      \end{array}
      \right.
      \end{equation}
      
      [[18]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix}     %const
      0.0821 \\ -0.5559
      \end{smatrix}+\begin{smatrix}     %trend
      0.0007 \\ -0.0018
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.9823 & -0.0024 \\
      0.7181 & 1.0032 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}< 1.1439 \\
      \begin{smatrix}     %const
      -0.0213 \\ 2.0413
      \end{smatrix}+\begin{smatrix}     %trend
      -0.0002 \\ 0.0064
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.9952 & 0.0008 \\
      -0.6686 & 0.9732 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if } 1.1439 < \text{Th} < 1.3893 \\
      \begin{smatrix}     %const
      0.3579 \\ -1.2736
      \end{smatrix}+\begin{smatrix}     %trend
      0.0014 \\ -0.0035
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.8674 & -0.0056 \\
      0.1600 & 1.0215 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}> 1.3893 \\
      \end{array}
      \right.
      \end{equation}
      
      [[19]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix}     %const
      0.0077 \\ -1.0125
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.9940 & 1.4\text{e-05} \\
      1.3416 & 0.9969 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}< 1.2252 \\
      \begin{smatrix}     %const
      0.0305 \\ 0.4650
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.9689 & 0.0001 \\
      -0.2413 & 1.0007 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}> 1.2252 \\
      \end{array}
      \right.
      \end{equation}
      
      [[20]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix}     %const
      -0.0298 \\ -0.2707
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      1.0311 & 8.4\text{e-06} \\
      0.5938 & 0.9971 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}< 1.1439 \\
      \begin{smatrix}     %const
      0.0672 \\ 1.9185
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.9434 & 1.1\text{e-05} \\
      -1.0550 & 0.9957 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if } 1.1439 < \text{Th} < 1.2252 \\
      \begin{smatrix}     %const
      0.0305 \\ 0.4650
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.9689 & 0.0001 \\
      -0.2413 & 1.0007 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}> 1.2252 \\
      \end{array}
      \right.
      \end{equation}
      
      [[21]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      +\begin{smatrix}      %Lag1
      1.0018 & -6\text{e-06} \\
      0.3056 & 0.9995 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}< 1.2252 \\
      +\begin{smatrix}      %Lag1
      0.9949 & 9.2\text{e-05} \\
      0.1555 & 0.9999 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}> 1.2252 \\
      \end{array}
      \right.
      \end{equation}
      
      [[22]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      +\begin{smatrix}      %Lag1
      0.9995 & 8.1\text{e-05} \\
      0.3063 & 0.9977 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}< 1.1439 \\
      +\begin{smatrix}      %Lag1
      0.9998 & 1.9\text{e-05} \\
      0.5560 & 0.9960 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if } 1.1439 < \text{Th} < 1.2252 \\
      +\begin{smatrix}      %Lag1
      0.9949 & 9.2\text{e-05} \\
      0.1555 & 0.9999 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}> 1.2252 \\
      \end{array}
      \right.
      \end{equation}
      
      [[23]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix}     %trend
       2.1\text{e-05} \\ 0.0039
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      1.0038 & -8.5\text{e-05} \\
      0.6855 & 0.9846 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}< 1.2252 \\
      \begin{smatrix}     %trend
      -0.0002 \\ -0.0014
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      0.9703 & 0.0010 \\
      -0.0072 & 1.0059 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}> 1.2252 \\
      \end{array}
      \right.
      \end{equation}
      
      [[24]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix}     %trend
      0.0003 \\ 0.0010
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      1.0280 & -0.0010 \\
      0.4086 & 0.9938 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}< 1.1439 \\
      \begin{smatrix}     %trend
      -0.0002 \\ 0.0017
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      0.9847 & 0.0007 \\
      0.7017 & 0.9899 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if } 1.1439 < \text{Th} < 1.2252 \\
      \begin{smatrix}     %trend
      -0.0002 \\ -0.0014
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      0.9703 & 0.0010 \\
      -0.0072 & 1.0059 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}> 1.2252 \\
      \end{array}
      \right.
      \end{equation}
      
      [[25]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix}     %const
      0.0923 \\ -0.3846
      \end{smatrix}+\begin{smatrix}     %trend
      0.0007 \\ -0.0006
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      1.1164 & 0.0057 \\
      -0.0579 & 1.2098 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1477 & -0.0081 \\
      0.6665 & -0.2107 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}< 1.1439 \\
      \begin{smatrix}     %const
      0.0006 \\ 1.6399
      \end{smatrix}+\begin{smatrix}     %trend
      -0.0001 \\ 0.0060
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      1.1410 & -0.0035 \\
      0.5722 & 1.3343 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1549 & 0.0040 \\
      -1.1085 & -0.3583 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}> 1.1439 \\
      \end{array}
      \right.
      \end{equation}
      
      [[26]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix}     %const
      0.0923 \\ -0.3846
      \end{smatrix}+\begin{smatrix}     %trend
      0.0007 \\ -0.0006
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      1.1164 & 0.0057 \\
      -0.0579 & 1.2098 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1477 & -0.0081 \\
      0.6665 & -0.2107 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}< 1.1439 \\
      \begin{smatrix}     %const
      -0.0154 \\ 1.5002
      \end{smatrix}+\begin{smatrix}     %trend
      -0.0002 \\ 0.0054
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      1.1806 & 0.0010 \\
      0.5086 & 1.3202 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1865 & -0.0003 \\
      -0.9522 & -0.3425 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if } 1.1439 < \text{Th} < 1.3893 \\
      \begin{smatrix}     %const
      0.3044 \\ 1.2815
      \end{smatrix}+\begin{smatrix}     %trend
      0.0012 \\ 0.0092
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.9527 & -0.0297 \\
      2.0105 & 1.3074 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.0770 & 0.0253 \\
      -2.3912 & -0.3393 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}> 1.3893 \\
      \end{array}
      \right.
      \end{equation}
      
      [[27]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix}     %const
      0.0140 \\ -0.6473
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      1.1438 & 0.0026 \\
      0.8354 & 1.3948 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1574 & -0.0025 \\
      0.0150 & -0.3968 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}< 1.2252 \\
      \begin{smatrix}     %const
      0.0402 \\ 0.2830
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      1.1335 & -0.0102 \\
      0.4492 & 1.3110 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1713 & 0.0104 \\
      -0.5822 & -0.3106 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}> 1.2252 \\
      \end{array}
      \right.
      \end{equation}
      
      [[28]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix}     %const
      0.0140 \\ -0.6473
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      1.1438 & 0.0026 \\
      0.8354 & 1.3948 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1574 & -0.0025 \\
      0.0150 & -0.3968 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}< 1.2252 \\
      \begin{smatrix}     %const
      0.0233 \\ -0.0503
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      1.2136 & -0.0041 \\
      0.5533 & 1.2690 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.2374 & 0.0042 \\
      -0.3583 & -0.2696 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if } 1.2252 < \text{Th} < 1.3893 \\
      \begin{smatrix}     %const
      0.0580 \\ -0.6535
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.9631 & -0.0279 \\
      2.0919 & 1.3218 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.0270 & 0.0283 \\
      -1.9980 & -0.3159 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}> 1.3893 \\
      \end{array}
      \right.
      \end{equation}
      
      [[29]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      +\begin{smatrix}      %Lag1
      1.1453 & -0.0006 \\
      0.6625 & 1.4943 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1436 & 0.0005 \\
      -0.5054 & -0.4946 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}< 1.2263 \\
      +\begin{smatrix}      %Lag1
      1.1609 & -0.0071 \\
      0.5957 & 1.3145 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1656 & 0.0072 \\
      -0.4870 & -0.3146 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}> 1.2263 \\
      \end{array}
      \right.
      \end{equation}
      
      [[30]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      +\begin{smatrix}      %Lag1
      1.1964 & 0.0037 \\
      -0.0369 & 1.2243 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1986 & -0.0036 \\
      0.2958 & -0.2267 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}< 1.1439 \\
      +\begin{smatrix}      %Lag1
      1.1265 & 0.0035 \\
      0.9516 & 1.3785 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1289 & -0.0034 \\
      -0.6091 & -0.3809 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if } 1.1439 < \text{Th} < 1.2252 \\
      +\begin{smatrix}      %Lag1
      1.1589 & -0.0080 \\
      0.6281 & 1.3264 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1628 & 0.0081 \\
      -0.5220 & -0.3266 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}> 1.2252 \\
      \end{array}
      \right.
      \end{equation}
      
      [[31]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix}     %trend
      -1.6\text{e-05} \\ 0.0031
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      1.1450 & -0.0004 \\
      0.7035 & 1.4672 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1450 & 0.0005 \\
      -0.2401 & -0.4792 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}< 1.2263 \\
      \begin{smatrix}     %trend
      -0.0002 \\ -0.0006
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      1.1271 & -0.0092 \\
      0.5092 & 1.3092 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1574 & 0.0102 \\
      -0.4662 & -0.3068 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}> 1.2263 \\
      \end{array}
      \right.
      \end{equation}
      
      [[32]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix}     %trend
      -6.3\text{e-05} \\ 0.0030
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      1.1558 & 0.0015 \\
      1.0534 & 1.5076 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1604 & -0.0013 \\
      -0.6132 & -0.5193 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}< 1.2338 \\
      \begin{smatrix}     %trend
      -0.0002 \\ 0.0012
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      1.1844 & -0.0057 \\
      0.3427 & 1.1578 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.2111 & 0.0066 \\
      -0.0438 & -0.1633 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if } 1.2338 < \text{Th} < 1.3893 \\
      \begin{smatrix}     %trend
      -0.0002 \\ 0.0036
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      0.9859 & -0.0276 \\
      2.1499 & 1.3166 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.0234 & 0.0286 \\
      -2.1654 & -0.3255 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}> 1.3893 \\
      \end{array}
      \right.
      \end{equation}
      
      [[33]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix} %ECT
      -0.0517 \\ 0.7379
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      0.0053 \\ 0.6239
      \end{smatrix}+\begin{smatrix}     %trend
      -3.9\text{e-05} \\ -0.0011
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.1675 & -0.0122 \\
      -0.1020 & 0.1395 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}< -0.1176 \\
      \begin{smatrix} %ECT
      -0.0010 \\ 0.0155
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      0.0036 \\ 0.1013
      \end{smatrix}+\begin{smatrix}     %trend
      -3.3\text{e-05} \\ -3.6\text{e-05}
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.1313 & 0.0023 \\
      0.2371 & 0.6209 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}> -0.1176 \\
      \end{array}
      \right.
      \end{equation}
      
      [[34]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix} %ECT
      -0.0517 \\ 0.7379
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      0.0053 \\ 0.6239
      \end{smatrix}+\begin{smatrix}     %trend
      -3.9\text{e-05} \\ -0.0011
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.1675 & -0.0122 \\
      -0.1020 & 0.1395 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}< -0.1176 \\
      \begin{smatrix} %ECT
      -0.0023 \\ -0.3376
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      0.0144 \\ 0.3660
      \end{smatrix}+\begin{smatrix}     %trend
      -0.0001 \\ -0.0012
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.0850 & 0.0035 \\
      1.1887 & 0.4816 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if } -0.1176 < \text{Th} < 0.3236 \\
      \begin{smatrix} %ECT
      -0.0290 \\ 0.3350
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      0.0168 \\ -0.1220
      \end{smatrix}+\begin{smatrix}     %trend
      3.6\text{e-06} \\ 0.0020
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.1850 & -0.0078 \\
      -1.5210 & 0.6142 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}> 0.3236 \\
      \end{array}
      \right.
      \end{equation}
      
      [[35]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix} %ECT
      -0.0375 \\ 0.7085
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      -0.0069 \\ 0.4007
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.1251 & -0.0102 \\
      0.1907 & 0.1982 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}< -0.1176 \\
      \begin{smatrix} %ECT
      0.0061 \\ 0.0587
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      -0.0016 \\ 0.0858
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.1644 & 0.0029 \\
      0.2257 & 0.6190 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}> -0.1176 \\
      \end{array}
      \right.
      \end{equation}
      
      [[36]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix} %ECT
      -0.0375 \\ 0.7085
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      -0.0069 \\ 0.4007
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.1251 & -0.0102 \\
      0.1907 & 0.1982 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}< -0.1205 \\
      \begin{smatrix} %ECT
      0.0415 \\ -0.3566
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      -0.0065 \\ 0.1845
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.2581 & 0.0112 \\
      0.4801 & 0.2797 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if } -0.1205 < \text{Th} < 0.1832 \\
      \begin{smatrix} %ECT
      -0.0007 \\ -0.0967
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      0.0008 \\ 0.1360
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.1052 & 0.0030 \\
      -0.7273 & 0.6801 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}> 0.1832 \\
      \end{array}
      \right.
      \end{equation}
      
      [[37]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix} %ECT
      -0.0220 \\ -0.5507
      \end{smatrix}ECT_{-1}+
      +\begin{smatrix}      %Lag1
      0.1177 & -0.0145 \\
      1.0162 & 0.3249 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}< -0.1176 \\
      \begin{smatrix} %ECT
      0.0028 \\ 0.1897
      \end{smatrix}ECT_{-1}+
      +\begin{smatrix}      %Lag1
      0.1680 & 0.0010 \\
      -0.2859 & 0.7504 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}> -0.1176 \\
      \end{array}
      \right.
      \end{equation}
      
      [[38]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix} %ECT
      -0.0252 \\ -0.4527
      \end{smatrix}ECT_{-1}+
      +\begin{smatrix}      %Lag1
      0.1674 & -0.0180 \\
      2.2624 & 0.2430 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}< -0.2599 \\
      \begin{smatrix} %ECT
      -0.0155 \\ -0.9222
      \end{smatrix}ECT_{-1}+
      +\begin{smatrix}      %Lag1
      0.0733 & -0.0093 \\
      -0.5414 & 0.3017 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if } -0.2599 < \text{Th} < -0.1051 \\
      \begin{smatrix} %ECT
      0.0030 \\ 0.1947
      \end{smatrix}ECT_{-1}+
      +\begin{smatrix}      %Lag1
      0.1689 & 0.0007 \\
      -0.2687 & 0.7443 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}> -0.1051 \\
      \end{array}
      \right.
      \end{equation}
      
      [[39]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix} %ECT
      -0.0600 \\ -0.5717
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %trend
      -4.7\text{e-05} \\  4.0\text{e-06}
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      0.2434 & -0.0117 \\
      1.4443 & 0.1073 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}< -0.2639 \\
      \begin{smatrix} %ECT
      0.0024 \\ 0.1882
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %trend
      9.1\text{e-08} \\ 0.0006
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      0.1271 & 0.0008 \\
      -0.1181 & 0.5950 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}> -0.2639 \\
      \end{array}
      \right.
      \end{equation}
      
      [[40]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix} %ECT
      -0.0659 \\ -0.7485
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %trend
      -4.9\text{e-05} \\ -0.0003
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      0.1448 & -0.0165 \\
      2.2377 & 0.2192 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}< -0.2796 \\
      \begin{smatrix} %ECT
      0.0025 \\ 0.0803
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %trend
      -3.1\text{e-07} \\ 0.0008
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      0.1908 & 0.0022 \\
      -0.4850 & 0.3836 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if } -0.2796 < \text{Th} < 0.2373 \\
      \begin{smatrix} %ECT
      0.0043 \\ 0.1290
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %trend
      -1.6\text{e-05} \\ 0.0008
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      0.0840 & 0.0010 \\
      -0.6772 & 0.7079 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}> 0.2373 \\
      \end{array}
      \right.
      \end{equation}
      
      [[41]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix} %ECT
      -0.0152 \\ -0.2511
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      0.0111 \\ 0.3347
      \end{smatrix}+\begin{smatrix}     %trend
      -4.8\text{e-05} \\ -0.0008
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.1839 & -0.0045 \\
      0.5328 & 0.3380 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      0.0294 & -0.0020 \\
      -1.2293 & -0.1620 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}< 0.2123 \\
      \begin{smatrix} %ECT
      -0.0128 \\ 0.1822
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      0.0099 \\ -0.0415
      \end{smatrix}+\begin{smatrix}     %trend
      -3.5\text{e-05} \\ 0.0011
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.1255 & 0.0058 \\
      -1.0680 & 0.4181 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1980 & -0.0106 \\
      0.9412 & 0.3071 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}> 0.2123 \\
      \end{array}
      \right.
      \end{equation}
      
      [[42]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix} %ECT
      -0.0174 \\ -0.1891
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      0.0112 \\ 0.2841
      \end{smatrix}+\begin{smatrix}     %trend
      -4.9\text{e-05} \\ -0.0006
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.1571 & -0.0103 \\
      0.4259 & 0.2873 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      0.0524 & 0.0012 \\
      -1.1874 & -0.1195 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}< 0.1839 \\
      \begin{smatrix} %ECT
      -0.3473 \\ -6.1188
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      0.0752 \\ 2.7091
      \end{smatrix}+\begin{smatrix}     %trend
      -5.3\text{e-05} \\ -0.0090
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.1603 & 0.0336 \\
      -1.9196 & 0.4390 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      0.1327 & -0.0304 \\
      -4.9985 & -0.2457 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if } 0.1839 < \text{Th} < 0.2161 \\
      \begin{smatrix} %ECT
      -0.0119 \\ 0.3550
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      0.0097 \\ -0.1138
      \end{smatrix}+\begin{smatrix}     %trend
      2.2\text{e-05} \\ 0.0019
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.1143 & -0.0028 \\
      -1.3484 & 0.3682 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.2251 & -0.0106 \\
      0.4193 & 0.2529 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}> 0.2161 \\
      \end{array}
      \right.
      \end{equation}
      
      [[43]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix} %ECT
      -0.0027 \\ -0.0309
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      0.0019 \\ 0.1737
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.1822 & -0.0033 \\
      0.5033 & 0.3591 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      0.0229 & -0.0008 \\
      -1.3438 & -0.1393 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}< 0.2123 \\
      \begin{smatrix} %ECT
      -0.0019 \\ -0.1715
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      0.0041 \\ 0.1464
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.1174 & 0.0054 \\
      -0.8027 & 0.4307 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.2105 & -0.0110 \\
      1.3485 & 0.3192 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}> 0.2123 \\
      \end{array}
      \right.
      \end{equation}
      
      [[44]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix} %ECT
      -0.0462 \\ 0.9601
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      -0.0101 \\ 0.5190
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.1582 & -0.0112 \\
      -0.3723 & 0.2068 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      0.0184 & 0.0035 \\
      -1.6020 & -0.1615 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}< -0.2136 \\
      \begin{smatrix} %ECT
      0.0204 \\ -0.1547
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      -0.0037 \\ 0.1421
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.1867 & 0.0172 \\
      1.4258 & 0.5895 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      0.0198 & -0.0100 \\
      -1.1805 & -0.1352 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if } -0.2136 < \text{Th} < 0.2161 \\
      \begin{smatrix} %ECT
      -0.0176 \\ -0.1413
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      0.0127 \\ 0.1422
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.1203 & -0.0022 \\
      -0.8278 & 0.4228 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.2166 & -0.0101 \\
      1.1683 & 0.2976 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}> 0.2161 \\
      \end{array}
      \right.
      \end{equation}
      
      [[45]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix} %ECT
      -0.0179 \\ -0.5653
      \end{smatrix}ECT_{-1}+
      +\begin{smatrix}      %Lag1
      0.1392 & -0.0133 \\
      0.6576 & 0.3145 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      0.0022 & 0.0014 \\
      -0.6239 & -0.0880 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}< -0.1898 \\
      \begin{smatrix} %ECT
      0.0065 \\ 0.1547
      \end{smatrix}ECT_{-1}+
      +\begin{smatrix}      %Lag1
      0.1600 & 0.0079 \\
      -0.0089 & 0.6712 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.0661 & -0.0107 \\
      -0.9158 & 0.1318 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}> -0.1898 \\
      \end{array}
      \right.
      \end{equation}
      
      [[46]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix} %ECT
      -0.0168 \\ -0.6306
      \end{smatrix}ECT_{-1}+
      +\begin{smatrix}      %Lag1
      0.2488 & -0.0090 \\
      0.7049 & 0.1132 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1622 & -0.0005 \\
      0.7159 & -0.1447 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}< -0.2834 \\
      \begin{smatrix} %ECT
      -0.0150 \\ -0.7096
      \end{smatrix}ECT_{-1}+
      +\begin{smatrix}      %Lag1
      0.0339 & -0.0261 \\
      -0.3325 & 0.4394 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      0.2262 & 0.0154 \\
      -1.4005 & -0.0329 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if } -0.2834 < \text{Th} < -0.2136 \\
      \begin{smatrix} %ECT
      0.0058 \\ 0.1521
      \end{smatrix}ECT_{-1}+
      +\begin{smatrix}      %Lag1
      0.1726 & 0.0105 \\
      -0.1469 & 0.6396 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.0808 & -0.0125 \\
      -1.1324 & 0.1646 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}> -0.2136 \\
      \end{array}
      \right.
      \end{equation}
      
      [[47]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix} %ECT
      -0.0016 \\ 0.1606
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %trend
      5.7\text{e-06} \\ 0.0008
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      0.1857 & -0.0020 \\
      0.5883 & 0.4155 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      0.0235 & 0.0006 \\
      -1.4090 & -0.0821 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}< 0.2123 \\
      \begin{smatrix} %ECT
      0.0060 \\ 0.1033
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %trend
      7.3\text{e-06} \\ 0.0010
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      0.1190 & 0.0061 \\
      -1.0407 & 0.4170 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.2116 & -0.0097 \\
      0.9983 & 0.3033 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}> 0.2123 \\
      \end{array}
      \right.
      \end{equation}
      
      [[48]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix} %ECT
      -0.0030 \\ 0.1749
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %trend
      6.1\text{e-06} \\ 0.0008
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      0.1616 & -0.0076 \\
      0.5407 & 0.3552 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      0.0466 & 0.0041 \\
      -1.3326 & -0.0463 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}< 0.1839 \\
      \begin{smatrix} %ECT
      -0.0599 \\ 4.2280
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %trend
      6.6\text{e-05} \\ -0.0047
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      0.1752 & 0.0419 \\
      -1.3831 & 0.7378 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      0.2040 & -0.0341 \\
      -2.4316 & -0.3810 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if } 0.1839 < \text{Th} < 0.2161 \\
      \begin{smatrix} %ECT
      0.0068 \\ 0.1360
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %trend
      6.5\text{e-05} \\ 0.0014
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      0.1067 & -0.0028 \\
      -1.2594 & 0.3674 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.2407 & -0.0101 \\
      0.6007 & 0.2467 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}> 0.2161 \\
      \end{array}
      \right.
      \end{equation}
      

---

    Code
      options(show.signif.stars = FALSE)

---

    Code
      sapply(mods, function(x) toLatex(summary(x)))
    Output
      [[1]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \begin{smatrix}     %const
      0.0104(0.0125) \\ 0.3790(0.1702)
      \end{smatrix}+\begin{smatrix}     %trend
      -2.2\text{e-05}(8.7\text{e-05}) \\ -0.0007(0.0012)
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.9871(0.0089) & 0.0001(0.0003) \\
      -0.1395(0.1216) & 1.0023(0.0045) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      \end{equation}
      
      [[2]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \begin{smatrix}     %const
      0.0128(0.0082) \\ 0.4554(0.1120)
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.9871(0.0089) & 6.4\text{e-05}(5.1\text{e-05}) \\
      -0.1411(0.1215) & 0.9996(0.0007) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      \end{equation}
      
      [[3]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      +\begin{smatrix}      %Lag1
      1.0004(0.0024) & 1.2\text{e-05}(3.9\text{e-05}) \\
      0.3346(0.0333) & 0.9978(0.0005) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      \end{equation}
      
      [[4]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \begin{smatrix}     %trend
      -7.7\text{e-05}(5.7\text{e-05}) \\ -0.0027(0.0008)
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      0.9920(0.0068) & 0.0003(0.0002) \\
      0.0365(0.0929) & 1.0090(0.0033) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      \end{equation}
      
      [[5]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \begin{smatrix}     %const
      0.0136(0.0125) \\ 0.2832(0.1537)
      \end{smatrix}+\begin{smatrix}     %trend
      -2.6\text{e-05}(8.7\text{e-05}) \\ 0.0003(0.0011)
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      1.1504(0.0554) & -0.0037(0.0041) \\
      0.5332(0.6837) & 1.4575(0.0502) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1661(0.0554) & 0.0039(0.0041) \\
      -0.6049(0.6833) & -0.4591(0.0504) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      \end{equation}
      
      [[6]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \begin{smatrix}     %const
      0.0164(0.0084) \\ 0.2463(0.1031)
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      1.1507(0.0553) & -0.0037(0.0041) \\
      0.5294(0.6827) & 1.4571(0.0501) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1664(0.0553) & 0.0038(0.0041) \\
      -0.6006(0.6822) & -0.4574(0.0501) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      \end{equation}
      
      [[7]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      +\begin{smatrix}      %Lag1
      1.1610(0.0553) & -0.0019(0.0040) \\
      0.6843(0.6846) & 1.4836(0.0492) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1602(0.0555) & 0.0020(0.0040) \\
      -0.5073(0.6861) & -0.4848(0.0491) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      \end{equation}
      
      [[8]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \begin{smatrix}     %trend
      -9.6\text{e-05}(5.8\text{e-05}) \\ -0.0011(0.0007)
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      1.1535(0.0554) & -0.0032(0.0040) \\
      0.5978(0.6854) & 1.4694(0.0499) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1629(0.0553) & 0.0036(0.0041) \\
      -0.5384(0.6849) & -0.4660(0.0505) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      \end{equation}
      
      [[9]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      +\begin{smatrix}  %ECT
      -0.0151(0.0084) \\ -0.0503(0.1041)
      \end{smatrix}ECT_{-1}
      \begin{smatrix}     %const
      0.0110(0.0053) \\ 0.1931(0.0650)
      \end{smatrix}+\begin{smatrix}     %trend
      -4.5\text{e-05}(2.6\text{e-05}) \\ -0.0003(0.0003)
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.1656(0.0553) & -0.0040(0.0041) \\
      0.5872(0.6821) & 0.4569(0.0503) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      \end{equation}
      
      [[10]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      +\begin{smatrix}  %ECT
      -0.0011(0.0025) \\ 0.0466(0.0310)
      \end{smatrix}ECT_{-1}
      \begin{smatrix}     %const
      0.0021(0.0013) \\ 0.1316(0.0156)
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.1595(0.0553) & -0.0030(0.0040) \\
      0.5453(0.6807) & 0.4636(0.0498) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      \end{equation}
      
      [[11]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      +\begin{smatrix}  %ECT
      -0.0007(0.0025) \\ 0.0732(0.0340)
      \end{smatrix}ECT_{-1}
      +\begin{smatrix}      %Lag1
      0.1681(0.0552) & 0.0024(0.0024) \\
      1.0771(0.7488) & 0.8004(0.0329) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      \end{equation}
      
      [[12]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      +\begin{smatrix}  %ECT
      0.0014(0.0030) \\ 0.2386(0.0376)
      \end{smatrix}ECT_{-1}
      \begin{smatrix}     %trend
      7.6\text{e-06}(6.3\text{e-06}) \\ 0.0006(7.8\text{e-05})
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      0.1611(0.0555) & -0.0012(0.0039) \\
      0.5082(0.6900) & 0.5052(0.0481) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      \end{equation}
      
      [[13]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      +\begin{smatrix}  %ECT
      -0.0152(0.0085) \\ -0.0398(0.1052)
      \end{smatrix}ECT_{-1}
      \begin{smatrix}     %const
      0.0116(0.0054) \\ 0.1849(0.0667)
      \end{smatrix}+\begin{smatrix}     %trend
      -4.6\text{e-05}(2.6\text{e-05}) \\ -0.0003(0.0003)
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.1679(0.0562) & -0.0025(0.0046) \\
      0.6918(0.6934) & 0.4539(0.0565) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.0220(0.0564) & -0.0031(0.0046) \\
      -0.6301(0.6962) & 0.0110(0.0567) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      \end{equation}
      
      [[14]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      +\begin{smatrix}  %ECT
      -0.0010(0.0025) \\ 0.0461(0.0313)
      \end{smatrix}ECT_{-1}
      \begin{smatrix}     %const
      0.0025(0.0014) \\ 0.1298(0.0175)
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.1635(0.0563) & -0.0019(0.0046) \\
      0.6656(0.6924) & 0.4578(0.0563) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.0303(0.0564) & -0.0023(0.0046) \\
      -0.6808(0.6933) & 0.0158(0.0564) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      \end{equation}
      
      [[15]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      +\begin{smatrix}  %ECT
      -0.0008(0.0025) \\ 0.0609(0.0338)
      \end{smatrix}ECT_{-1}
      +\begin{smatrix}      %Lag1
      0.1723(0.0563) & 0.0015(0.0041) \\
      1.1289(0.7465) & 0.6385(0.0550) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.0252(0.0565) & 0.0012(0.0041) \\
      -0.4086(0.7495) & 0.2017(0.0547) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      \end{equation}
      
      [[16]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      +\begin{smatrix}  %ECT
      0.0017(0.0032) \\ 0.2304(0.0400)
      \end{smatrix}ECT_{-1}
      \begin{smatrix}     %trend
      8.4\text{e-06}(7.0\text{e-06}) \\ 0.0006(8.6\text{e-05})
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      0.1656(0.0565) & -0.0007(0.0045) \\
      0.6551(0.7006) & 0.4838(0.0561) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.0302(0.0566) & -0.0010(0.0045) \\
      -0.7619(0.7019) & 0.0441(0.0560) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      \end{equation}
      
      [[17]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix}     %const
      0.0821(0.0789) \\ -0.5559(0.9804)
      \end{smatrix}+\begin{smatrix}     %trend
      0.0007(0.0004) \\ -0.0018(0.0056)
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.9823(0.0483) & -0.0024(0.0015) \\
      0.7181(0.6005) & 1.0032(0.0191) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}< 1.1439 \\
      \begin{smatrix}     %const
      -0.0135(0.0215) \\ 2.2021(0.2674)
      \end{smatrix}+\begin{smatrix}     %trend
      -0.0002(0.0001) \\ 0.0073(0.0015)
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.9915(0.0107) & 0.0007(0.0005) \\
      -0.7468(0.1326) & 0.9702(0.0059) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}> 1.1439 \\
      \end{array}
      \right.
      \end{equation}
      
      [[18]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix}     %const
      0.0821(0.0790) \\ -0.5559(0.9641)
      \end{smatrix}+\begin{smatrix}     %trend
      0.0007(0.0004) \\ -0.0018(0.0055)
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.9823(0.0484) & -0.0024(0.0015) \\
      0.7181(0.5905) & 1.0032(0.0188) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}< 1.1439 \\
      \begin{smatrix}     %const
      -0.0213(0.0222) \\ 2.0413(0.2707)
      \end{smatrix}+\begin{smatrix}     %trend
      -0.0002(0.0001) \\ 0.0064(0.0016)
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.9952(0.0124) & 0.0008(0.0005) \\
      -0.6686(0.1518) & 0.9732(0.0063) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if } 1.1439 < \text{Th} < 1.3893 \\
      \begin{smatrix}     %const
      0.3579(0.2290) \\ -1.2736(2.7959)
      \end{smatrix}+\begin{smatrix}     %trend
      0.0014(0.0010) \\ -0.0035(0.0126)
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.8674(0.0769) & -0.0056(0.0042) \\
      0.1600(0.9391) & 1.0215(0.0519) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}> 1.3893 \\
      \end{array}
      \right.
      \end{equation}
      
      [[19]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix}     %const
      0.0077(0.0171) \\ -1.0125(0.2138)
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.9940(0.0178) & 1.4\text{e-05}(6.9\text{e-05}) \\
      1.3416(0.2218) & 0.9969(0.0009) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}< 1.2252 \\
      \begin{smatrix}     %const
      0.0305(0.0196) \\ 0.4650(0.2452)
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.9689(0.0175) & 0.0001(8.2\text{e-05}) \\
      -0.2413(0.2179) & 1.0007(0.0010) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}> 1.2252 \\
      \end{array}
      \right.
      \end{equation}
      
      [[20]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix}     %const
      -0.0298(0.0348) \\ -0.2707(0.4248)
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      1.0311(0.0373) & 8.4\text{e-06}(0.0002) \\
      0.5938(0.4552) & 0.9971(0.0019) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}< 1.1439 \\
      \begin{smatrix}     %const
      0.0672(0.0777) \\ 1.9185(0.9477)
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.9434(0.0654) & 1.1\text{e-05}(8.1\text{e-05}) \\
      -1.0550(0.7981) & 0.9957(0.0010) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if } 1.1439 < \text{Th} < 1.2252 \\
      \begin{smatrix}     %const
      0.0305(0.0197) \\ 0.4650(0.2399)
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.9689(0.0175) & 0.0001(8.2\text{e-05}) \\
      -0.2413(0.2131) & 1.0007(0.0010) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}> 1.2252 \\
      \end{array}
      \right.
      \end{equation}
      
      [[21]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      +\begin{smatrix}      %Lag1
      1.0018(0.0029) & -6\text{e-06}(5.4\text{e-05}) \\
      0.3056(0.0377) & 0.9995(0.0007) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}< 1.2252 \\
      +\begin{smatrix}      %Lag1
      0.9949(0.0049) & 9.2\text{e-05}(7.4\text{e-05}) \\
      0.1555(0.0630) & 0.9999(0.0010) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}> 1.2252 \\
      \end{array}
      \right.
      \end{equation}
      
      [[22]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      +\begin{smatrix}      %Lag1
      0.9995(0.0051) & 8.1\text{e-05}(0.0001) \\
      0.3063(0.0627) & 0.9977(0.0016) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}< 1.1439 \\
      +\begin{smatrix}      %Lag1
      0.9998(0.0049) & 1.9\text{e-05}(8.1\text{e-05}) \\
      0.5560(0.0606) & 0.9960(0.0010) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if } 1.1439 < \text{Th} < 1.2252 \\
      +\begin{smatrix}      %Lag1
      0.9949(0.0049) & 9.2\text{e-05}(7.4\text{e-05}) \\
      0.1555(0.0599) & 0.9999(0.0009) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}> 1.2252 \\
      \end{array}
      \right.
      \end{equation}
      
      [[23]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix}     %trend
       2.1\text{e-05}(0.0001) \\ 0.0039(0.0017)
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      1.0038(0.0130) & -8.5\text{e-05}(0.0005) \\
      0.6855(0.1673) & 0.9846(0.0064) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}< 1.2252 \\
      \begin{smatrix}     %trend
      -0.0002(9.4\text{e-05}) \\ -0.0014(0.0012)
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      0.9703(0.0121) & 0.0010(0.0004) \\
      -0.0072(0.1564) & 1.0059(0.0054) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}> 1.2252 \\
      \end{array}
      \right.
      \end{equation}
      
      [[24]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix}     %trend
      0.0003(0.0002) \\ 0.0010(0.0024)
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      1.0280(0.0201) & -0.0010(0.0008) \\
      0.4086(0.2491) & 0.9938(0.0094) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}< 1.1439 \\
      \begin{smatrix}     %trend
      -0.0002(0.0002) \\ 0.0017(0.0024)
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      0.9847(0.0175) & 0.0007(0.0007) \\
      0.7017(0.2173) & 0.9899(0.0088) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if } 1.1439 < \text{Th} < 1.2252 \\
      \begin{smatrix}     %trend
      -0.0002(9.4\text{e-05}) \\ -0.0014(0.0012)
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      0.9703(0.0121) & 0.0010(0.0004) \\
      -0.0072(0.1501) & 1.0059(0.0051) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}> 1.2252 \\
      \end{array}
      \right.
      \end{equation}
      
      [[25]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix}     %const
      0.0923(0.0822) \\ -0.3846(0.9608)
      \end{smatrix}+\begin{smatrix}     %trend
      0.0007(0.0005) \\ -0.0006(0.0055)
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      1.1164(0.1652) & 0.0057(0.0160) \\
      -0.0579(1.9317) & 1.2098(0.1872) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1477(0.1670) & -0.0081(0.0163) \\
      0.6665(1.9529) & -0.2107(0.1901) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}< 1.1439 \\
      \begin{smatrix}     %const
      0.0006(0.0231) \\ 1.6399(0.2701)
      \end{smatrix}+\begin{smatrix}     %trend
      -0.0001(0.0001) \\ 0.0060(0.0014)
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      1.1410(0.0596) & -0.0035(0.0047) \\
      0.5722(0.6974) & 1.3343(0.0546) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1549(0.0597) & 0.0040(0.0046) \\
      -1.1085(0.6980) & -0.3583(0.0533) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}> 1.1439 \\
      \end{array}
      \right.
      \end{equation}
      
      [[26]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix}     %const
      0.0923(0.0817) \\ -0.3846(0.9528)
      \end{smatrix}+\begin{smatrix}     %trend
      0.0007(0.0005) \\ -0.0006(0.0054)
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      1.1164(0.1642) & 0.0057(0.0159) \\
      -0.0579(1.9156) & 1.2098(0.1856) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1477(0.1660) & -0.0081(0.0162) \\
      0.6665(1.9366) & -0.2107(0.1885) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}< 1.1439 \\
      \begin{smatrix}     %const
      -0.0154(0.0238) \\ 1.5002(0.2773)
      \end{smatrix}+\begin{smatrix}     %trend
      -0.0002(0.0001) \\ 0.0054(0.0016)
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      1.1806(0.0687) & 0.0010(0.0052) \\
      0.5086(0.8014) & 1.3202(0.0602) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1865(0.0681) & -0.0003(0.0051) \\
      -0.9522(0.7944) & -0.3425(0.0589) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if } 1.1439 < \text{Th} < 1.3893 \\
      \begin{smatrix}     %const
      0.3044(0.2443) \\ 1.2815(2.8498)
      \end{smatrix}+\begin{smatrix}     %trend
      0.0012(0.0011) \\ 0.0092(0.0130)
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.9527(0.1378) & -0.0297(0.0111) \\
      2.0105(1.6081) & 1.3074(0.1297) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.0770(0.1330) & 0.0253(0.0114) \\
      -2.3912(1.5514) & -0.3393(0.1329) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}> 1.3893 \\
      \end{array}
      \right.
      \end{equation}
      
      [[27]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix}     %const
      0.0140(0.0181) \\ -0.6473(0.2154)
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      1.1438(0.0902) & 0.0026(0.0065) \\
      0.8354(1.0711) & 1.3948(0.0772) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1574(0.0911) & -0.0025(0.0065) \\
      0.0150(1.0824) & -0.3968(0.0770) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}< 1.2252 \\
      \begin{smatrix}     %const
      0.0402(0.0198) \\ 0.2830(0.2353)
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      1.1335(0.0713) & -0.0102(0.0060) \\
      0.4492(0.8467) & 1.3110(0.0713) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1713(0.0702) & 0.0104(0.0060) \\
      -0.5822(0.8334) & -0.3106(0.0714) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}> 1.2252 \\
      \end{array}
      \right.
      \end{equation}
      
      [[28]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix}     %const
      0.0140(0.0181) \\ -0.6473(0.2125)
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      1.1438(0.0900) & 0.0026(0.0065) \\
      0.8354(1.0568) & 1.3948(0.0762) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1574(0.0910) & -0.0025(0.0065) \\
      0.0150(1.0678) & -0.3968(0.0759) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}< 1.2252 \\
      \begin{smatrix}     %const
      0.0233(0.0295) \\ -0.0503(0.3461)
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      1.2136(0.0890) & -0.0041(0.0072) \\
      0.5533(1.0451) & 1.2690(0.0850) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.2374(0.0858) & 0.0042(0.0072) \\
      -0.3583(1.0068) & -0.2696(0.0851) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if } 1.2252 < \text{Th} < 1.3893 \\
      \begin{smatrix}     %const
      0.0580(0.0681) \\ -0.6535(0.7991)
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.9631(0.1377) & -0.0279(0.0110) \\
      2.0919(1.6169) & 1.3218(0.1292) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.0270(0.1244) & 0.0283(0.0110) \\
      -1.9980(1.4599) & -0.3159(0.1297) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}> 1.3893 \\
      \end{array}
      \right.
      \end{equation}
      
      [[29]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      +\begin{smatrix}      %Lag1
      1.1453(0.0902) & -0.0006(0.0059) \\
      0.6625(1.0808) & 1.4943(0.0712) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1436(0.0905) & 0.0005(0.0059) \\
      -0.5054(1.0845) & -0.4946(0.0712) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}< 1.2263 \\
      +\begin{smatrix}      %Lag1
      1.1609(0.0705) & -0.0071(0.0060) \\
      0.5957(0.8445) & 1.3145(0.0719) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1656(0.0704) & 0.0072(0.0060) \\
      -0.4870(0.8432) & -0.3146(0.0719) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}> 1.2263 \\
      \end{array}
      \right.
      \end{equation}
      
      [[30]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      +\begin{smatrix}      %Lag1
      1.1964(0.1579) & 0.0037(0.0159) \\
      -0.0369(1.8495) & 1.2243(0.1863) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1986(0.1590) & -0.0036(0.0159) \\
      0.2958(1.8622) & -0.2267(0.1863) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}< 1.1439 \\
      +\begin{smatrix}      %Lag1
      1.1265(0.1110) & 0.0035(0.0075) \\
      0.9516(1.2997) & 1.3785(0.0879) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1289(0.1112) & -0.0034(0.0075) \\
      -0.6091(1.3020) & -0.3809(0.0876) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if } 1.1439 < \text{Th} < 1.2252 \\
      +\begin{smatrix}      %Lag1
      1.1589(0.0707) & -0.0080(0.0060) \\
      0.6281(0.8282) & 1.3264(0.0697) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1628(0.0706) & 0.0081(0.0060) \\
      -0.5220(0.8266) & -0.3266(0.0697) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}> 1.2252 \\
      \end{array}
      \right.
      \end{equation}
      
      [[31]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix}     %trend
      -1.6\text{e-05}(0.0001) \\ 0.0031(0.0016)
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      1.1450(0.0897) & -0.0004(0.0060) \\
      0.7035(1.0774) & 1.4672(0.0723) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1450(0.0907) & 0.0005(0.0059) \\
      -0.2401(1.0893) & -0.4792(0.0714) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}< 1.2263 \\
      \begin{smatrix}     %trend
      -0.0002(9.6\text{e-05}) \\ -0.0006(0.0012)
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      1.1271(0.0716) & -0.0092(0.0060) \\
      0.5092(0.8599) & 1.3092(0.0725) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1574(0.0701) & 0.0102(0.0061) \\
      -0.4662(0.8415) & -0.3068(0.0734) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}> 1.2263 \\
      \end{array}
      \right.
      \end{equation}
      
      [[32]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      X_{t}^{1} \\ X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix}     %trend
      -6.3\text{e-05}(0.0001) \\ 0.0030(0.0015)
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      1.1558(0.0889) & 0.0015(0.0058) \\
      1.0534(1.0500) & 1.5076(0.0686) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1604(0.0899) & -0.0013(0.0057) \\
      -0.6132(1.0612) & -0.5193(0.0677) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}< 1.2338 \\
      \begin{smatrix}     %trend
      -0.0002(0.0002) \\ 0.0012(0.0020)
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      1.1844(0.0926) & -0.0057(0.0079) \\
      0.3427(1.0931) & 1.1578(0.0929) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.2111(0.0877) & 0.0066(0.0079) \\
      -0.0438(1.0356) & -0.1633(0.0937) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if } 1.2338 < \text{Th} < 1.3893 \\
      \begin{smatrix}     %trend
      -0.0002(0.0003) \\ 0.0036(0.0037)
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      0.9859(0.1352) & -0.0276(0.0110) \\
      2.1499(1.5968) & 1.3166(0.1297) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-1}^{1} \\ X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.0234(0.1258) & 0.0286(0.0111) \\
      -2.1654(1.4854) & -0.3255(0.1308) 
      \end{smatrix}
      \begin{smatrix}
      X_{t-2}^{1} \\ X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}> 1.3893 \\
      \end{array}
      \right.
      \end{equation}
      
      [[33]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix} %ECT
      -0.0517(0.0284) \\ 0.7379(0.3404)
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      0.0053(0.0109) \\ 0.6239(0.1304)
      \end{smatrix}+\begin{smatrix}     %trend
      -3.9\text{e-05}(3.5\text{e-05}) \\ -0.0011(0.0004)
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.1675(0.0835) & -0.0122(0.0068) \\
      -0.1020(1.0027) & 0.1395(0.0819) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}< -0.1176 \\
      \begin{smatrix} %ECT
      -0.0010(0.0139) \\ 0.0155(0.1664)
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      0.0036(0.0083) \\ 0.1013(0.0996)
      \end{smatrix}+\begin{smatrix}     %trend
      -3.3\text{e-05}(3.8\text{e-05}) \\ -3.6\text{e-05}(0.0005)
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.1313(0.0747) & 0.0023(0.0051) \\
      0.2371(0.8964) & 0.6209(0.0613) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}> -0.1176 \\
      \end{array}
      \right.
      \end{equation}
      
      [[34]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix} %ECT
      -0.0517(0.0283) \\ 0.7379(0.3362)
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      0.0053(0.0108) \\ 0.6239(0.1288)
      \end{smatrix}+\begin{smatrix}     %trend
      -3.9\text{e-05}(3.5\text{e-05}) \\ -0.0011(0.0004)
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.1675(0.0832) & -0.0122(0.0068) \\
      -0.1020(0.9904) & 0.1395(0.0809) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}< -0.1176 \\
      \begin{smatrix} %ECT
      -0.0023(0.0168) \\ -0.3376(0.2003)
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      0.0144(0.0121) \\ 0.3660(0.1444)
      \end{smatrix}+\begin{smatrix}     %trend
      -0.0001(5.9\text{e-05}) \\ -0.0012(0.0007)
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.0850(0.0929) & 0.0035(0.0067) \\
      1.1887(1.1052) & 0.4816(0.0798) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if } -0.1176 < \text{Th} < 0.3236 \\
      \begin{smatrix} %ECT
      -0.0290(0.0368) \\ 0.3350(0.4381)
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      0.0168(0.0191) \\ -0.1220(0.2279)
      \end{smatrix}+\begin{smatrix}     %trend
      3.6\text{e-06}(8.2\text{e-05}) \\ 0.0020(0.0010)
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.1850(0.1286) & -0.0078(0.0106) \\
      -1.5210(1.5305) & 0.6142(0.1265) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}> 0.3236 \\
      \end{array}
      \right.
      \end{equation}
      
      [[35]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix} %ECT
      -0.0375(0.0183) \\ 0.7085(0.2200)
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      -0.0069(0.0060) \\ 0.4007(0.0724)
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.1251(0.0787) & -0.0102(0.0064) \\
      0.1907(0.9463) & 0.1982(0.0773) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}< -0.1176 \\
      \begin{smatrix} %ECT
      0.0061(0.0062) \\ 0.0587(0.0743)
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      -0.0016(0.0023) \\ 0.0858(0.0277)
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.1644(0.0783) & 0.0029(0.0052) \\
      0.2257(0.9415) & 0.6190(0.0627) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}> -0.1176 \\
      \end{array}
      \right.
      \end{equation}
      
      [[36]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix} %ECT
      -0.0375(0.0183) \\ 0.7085(0.2171)
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      -0.0069(0.0060) \\ 0.4007(0.0715)
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.1251(0.0786) & -0.0102(0.0064) \\
      0.1907(0.9340) & 0.1982(0.0763) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}< -0.1205 \\
      \begin{smatrix} %ECT
      0.0415(0.0199) \\ -0.3566(0.2361)
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      -0.0065(0.0040) \\ 0.1845(0.0473)
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.2581(0.1399) & 0.0112(0.0113) \\
      0.4801(1.6625) & 0.2797(0.1339) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if } -0.1205 < \text{Th} < 0.1832 \\
      \begin{smatrix} %ECT
      -0.0007(0.0111) \\ -0.0967(0.1321)
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      0.0008(0.0047) \\ 0.1360(0.0559)
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.1052(0.0982) & 0.0030(0.0062) \\
      -0.7273(1.1672) & 0.6801(0.0738) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}> 0.1832 \\
      \end{array}
      \right.
      \end{equation}
      
      [[37]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix} %ECT
      -0.0220(0.0069) \\ -0.5507(0.0868)
      \end{smatrix}ECT_{-1}+
      +\begin{smatrix}      %Lag1
      0.1177(0.0790) & -0.0145(0.0061) \\
      1.0162(1.0004) & 0.3249(0.0767) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}< -0.1176 \\
      \begin{smatrix} %ECT
      0.0028(0.0042) \\ 0.1897(0.0530)
      \end{smatrix}ECT_{-1}+
      +\begin{smatrix}      %Lag1
      0.1680(0.0770) & 0.0010(0.0042) \\
      -0.2859(0.9741) & 0.7504(0.0530) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}> -0.1176 \\
      \end{array}
      \right.
      \end{equation}
      
      [[38]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix} %ECT
      -0.0252(0.0094) \\ -0.4527(0.1153)
      \end{smatrix}ECT_{-1}+
      +\begin{smatrix}      %Lag1
      0.1674(0.1299) & -0.0180(0.0101) \\
      2.2624(1.5980) & 0.2430(0.1242) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}< -0.2599 \\
      \begin{smatrix} %ECT
      -0.0155(0.0118) \\ -0.9222(0.1449)
      \end{smatrix}ECT_{-1}+
      +\begin{smatrix}      %Lag1
      0.0733(0.1023) & -0.0093(0.0079) \\
      -0.5414(1.2582) & 0.3017(0.0973) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if } -0.2599 < \text{Th} < -0.1051 \\
      \begin{smatrix} %ECT
      0.0030(0.0043) \\ 0.1947(0.0524)
      \end{smatrix}ECT_{-1}+
      +\begin{smatrix}      %Lag1
      0.1689(0.0773) & 0.0007(0.0043) \\
      -0.2687(0.9509) & 0.7443(0.0527) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}> -0.1051 \\
      \end{array}
      \right.
      \end{equation}
      
      [[39]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix} %ECT
      -0.0600(0.0389) \\ -0.5717(0.4731)
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %trend
      -4.7\text{e-05}(4.3\text{e-05}) \\  4.0\text{e-06}(0.0005)
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      0.2434(0.1212) & -0.0117(0.0093) \\
      1.4443(1.4749) & 0.1073(0.1127) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}< -0.2639 \\
      \begin{smatrix} %ECT
      0.0024(0.0034) \\ 0.1882(0.0413)
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %trend
      9.1\text{e-08}(7.8\text{e-06}) \\ 0.0006(9.5\text{e-05})
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      0.1271(0.0625) & 0.0008(0.0043) \\
      -0.1181(0.7613) & 0.5950(0.0523) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}> -0.2639 \\
      \end{array}
      \right.
      \end{equation}
      
      [[40]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix} %ECT
      -0.0659(0.0410) \\ -0.7485(0.4878)
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %trend
      -4.9\text{e-05}(4.6\text{e-05}) \\ -0.0003(0.0005)
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      0.1448(0.1360) & -0.0165(0.0103) \\
      2.2377(1.6189) & 0.2192(0.1228) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}< -0.2796 \\
      \begin{smatrix} %ECT
      0.0025(0.0057) \\ 0.0803(0.0675)
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %trend
      -3.1\text{e-07}(8.7\text{e-06}) \\ 0.0008(0.0001)
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      0.1908(0.0748) & 0.0022(0.0057) \\
      -0.4850(0.8905) & 0.3836(0.0675) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if } -0.2796 < \text{Th} < 0.2373 \\
      \begin{smatrix} %ECT
      0.0043(0.0054) \\ 0.1290(0.0642)
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %trend
      -1.6\text{e-05}(3.5\text{e-05}) \\ 0.0008(0.0004)
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      0.0840(0.1103) & 0.0010(0.0071) \\
      -0.6772(1.3130) & 0.7079(0.0840) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      & \text{if Th}> 0.2373 \\
      \end{array}
      \right.
      \end{equation}
      
      [[41]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix} %ECT
      -0.0152(0.0094) \\ -0.2511(0.1096)
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      0.0111(0.0064) \\ 0.3347(0.0753)
      \end{smatrix}+\begin{smatrix}     %trend
      -4.8\text{e-05}(3.2\text{e-05}) \\ -0.0008(0.0004)
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.1839(0.0661) & -0.0045(0.0054) \\
      0.5328(0.7715) & 0.3380(0.0632) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      0.0294(0.0644) & -0.0020(0.0054) \\
      -1.2293(0.7525) & -0.1620(0.0626) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}< 0.2123 \\
      \begin{smatrix} %ECT
      -0.0128(0.0295) \\ 0.1822(0.3444)
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      0.0099(0.0152) \\ -0.0415(0.1775)
      \end{smatrix}+\begin{smatrix}     %trend
      -3.5\text{e-05}(7.8\text{e-05}) \\ 0.0011(0.0009)
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.1255(0.1144) & 0.0058(0.0107) \\
      -1.0680(1.3358) & 0.4181(0.1248) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1980(0.1223) & -0.0106(0.0113) \\
      0.9412(1.4283) & 0.3071(0.1314) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}> 0.2123 \\
      \end{array}
      \right.
      \end{equation}
      
      [[42]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix} %ECT
      -0.0174(0.0103) \\ -0.1891(0.1190)
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      0.0112(0.0066) \\ 0.2841(0.0757)
      \end{smatrix}+\begin{smatrix}     %trend
      -4.9\text{e-05}(3.4\text{e-05}) \\ -0.0006(0.0004)
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.1571(0.0693) & -0.0103(0.0058) \\
      0.4259(0.7969) & 0.2873(0.0664) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      0.0524(0.0690) & 0.0012(0.0058) \\
      -1.1874(0.7936) & -0.1195(0.0666) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}< 0.1839 \\
      \begin{smatrix} %ECT
      -0.3473(0.2701) \\ -6.1188(3.1069)
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      0.0752(0.0659) \\ 2.7091(0.7580)
      \end{smatrix}+\begin{smatrix}     %trend
      -5.3\text{e-05}(0.0002) \\ -0.0090(0.0018)
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.1603(0.2469) & 0.0336(0.0165) \\
      -1.9196(2.8400) & 0.4390(0.1899) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      0.1327(0.2603) & -0.0304(0.0161) \\
      -4.9985(2.9937) & -0.2457(0.1852) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if } 0.1839 < \text{Th} < 0.2161 \\
      \begin{smatrix} %ECT
      -0.0119(0.0302) \\ 0.3550(0.3470)
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      0.0097(0.0154) \\ -0.1138(0.1769)
      \end{smatrix}+\begin{smatrix}     %trend
      2.2\text{e-05}(8.8\text{e-05}) \\ 0.0019(0.0010)
      \end{smatrix}t
      +\begin{smatrix}      %Lag1
      0.1143(0.1149) & -0.0028(0.0114) \\
      -1.3484(1.3217) & 0.3682(0.1309) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.2251(0.1235) & -0.0106(0.0116) \\
      0.4193(1.4206) & 0.2529(0.1340) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}> 0.2161 \\
      \end{array}
      \right.
      \end{equation}
      
      [[43]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix} %ECT
      -0.0027(0.0040) \\ -0.0309(0.0473)
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      0.0019(0.0018) \\ 0.1737(0.0206)
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.1822(0.0661) & -0.0033(0.0054) \\
      0.5033(0.7770) & 0.3591(0.0629) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      0.0229(0.0643) & -0.0008(0.0053) \\
      -1.3438(0.7562) & -0.1393(0.0622) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}< 0.2123 \\
      \begin{smatrix} %ECT
      -0.0019(0.0172) \\ -0.1715(0.2017)
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      0.0041(0.0083) \\ 0.1464(0.0973)
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.1174(0.1130) & 0.0054(0.0107) \\
      -0.8027(1.3288) & 0.4307(0.1253) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.2105(0.1192) & -0.0110(0.0112) \\
      1.3485(1.4015) & 0.3192(0.1320) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}> 0.2123 \\
      \end{array}
      \right.
      \end{equation}
      
      [[44]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix} %ECT
      -0.0462(0.0266) \\ 0.9601(0.3168)
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      -0.0101(0.0088) \\ 0.5190(0.1049)
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.1582(0.0846) & -0.0112(0.0067) \\
      -0.3723(1.0060) & 0.2068(0.0800) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      0.0184(0.0828) & 0.0035(0.0072) \\
      -1.6020(0.9849) & -0.1615(0.0859) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}< -0.2136 \\
      \begin{smatrix} %ECT
      0.0204(0.0111) \\ -0.1547(0.1321)
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      -0.0037(0.0029) \\ 0.1421(0.0345)
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.1867(0.1087) & 0.0172(0.0085) \\
      1.4258(1.2918) & 0.5895(0.1006) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      0.0198(0.1063) & -0.0100(0.0080) \\
      -1.1805(1.2642) & -0.1352(0.0949) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if } -0.2136 < \text{Th} < 0.2161 \\
      \begin{smatrix} %ECT
      -0.0176(0.0191) \\ -0.1413(0.2271)
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %const
      0.0127(0.0096) \\ 0.1422(0.1137)
      \end{smatrix}
      +\begin{smatrix}      %Lag1
      0.1203(0.1116) & -0.0022(0.0110) \\
      -0.8278(1.3273) & 0.4228(0.1310) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.2166(0.1177) & -0.0101(0.0114) \\
      1.1683(1.3996) & 0.2976(0.1355) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}> 0.2161 \\
      \end{array}
      \right.
      \end{equation}
      
      [[45]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix} %ECT
      -0.0179(0.0077) \\ -0.5653(0.0972)
      \end{smatrix}ECT_{-1}+
      +\begin{smatrix}      %Lag1
      0.1392(0.0831) & -0.0133(0.0065) \\
      0.6576(1.0509) & 0.3145(0.0822) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      0.0022(0.0814) & 0.0014(0.0071) \\
      -0.6239(1.0296) & -0.0880(0.0897) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}< -0.1898 \\
      \begin{smatrix} %ECT
      0.0065(0.0042) \\ 0.1547(0.0528)
      \end{smatrix}ECT_{-1}+
      +\begin{smatrix}      %Lag1
      0.1600(0.0768) & 0.0079(0.0062) \\
      -0.0089(0.9717) & 0.6712(0.0781) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.0661(0.0776) & -0.0107(0.0060) \\
      -0.9158(0.9817) & 0.1318(0.0758) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}> -0.1898 \\
      \end{array}
      \right.
      \end{equation}
      
      [[46]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix} %ECT
      -0.0168(0.0094) \\ -0.6306(0.1192)
      \end{smatrix}ECT_{-1}+
      +\begin{smatrix}      %Lag1
      0.2488(0.1218) & -0.0090(0.0094) \\
      0.7049(1.5386) & 0.1132(0.1192) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.1622(0.1182) & -0.0005(0.0087) \\
      0.7159(1.4931) & -0.1447(0.1102) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}< -0.2834 \\
      \begin{smatrix} %ECT
      -0.0150(0.0158) \\ -0.7096(0.1996)
      \end{smatrix}ECT_{-1}+
      +\begin{smatrix}      %Lag1
      0.0339(0.1187) & -0.0261(0.0116) \\
      -0.3325(1.4989) & 0.4394(0.1463) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      0.2262(0.1210) & 0.0154(0.0131) \\
      -1.4005(1.5279) & -0.0329(0.1650) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if } -0.2834 < \text{Th} < -0.2136 \\
      \begin{smatrix} %ECT
      0.0058(0.0039) \\ 0.1521(0.0496)
      \end{smatrix}ECT_{-1}+
      +\begin{smatrix}      %Lag1
      0.1726(0.0747) & 0.0105(0.0057) \\
      -0.1469(0.9433) & 0.6396(0.0716) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.0808(0.0749) & -0.0125(0.0057) \\
      -1.1324(0.9454) & 0.1646(0.0716) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}> -0.2136 \\
      \end{array}
      \right.
      \end{equation}
      
      [[47]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix} %ECT
      -0.0016(0.0050) \\ 0.1606(0.0603)
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %trend
      5.7\text{e-06}(8.8\text{e-06}) \\ 0.0008(0.0001)
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      0.1857(0.0662) & -0.0020(0.0052) \\
      0.5883(0.7933) & 0.4155(0.0625) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      0.0235(0.0645) & 0.0006(0.0051) \\
      -1.4090(0.7727) & -0.0821(0.0617) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}< 0.2123 \\
      \begin{smatrix} %ECT
      0.0060(0.0057) \\ 0.1033(0.0687)
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %trend
      7.3\text{e-06}(4.2\text{e-05}) \\ 0.0010(0.0005)
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      0.1190(0.1142) & 0.0061(0.0107) \\
      -1.0407(1.3685) & 0.4170(0.1283) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.2116(0.1208) & -0.0097(0.0112) \\
      0.9983(1.4473) & 0.3033(0.1342) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}> 0.2123 \\
      \end{array}
      \right.
      \end{equation}
      
      [[48]]
      %insert in the preamble and uncomment the line you want for usual /medium /small matrix
      %\usepackage{amsmath} \newenvironment{smatrix}{\begin{pmatrix}}{\end{pmatrix}} %USUAL
      %\usepackage{amsmath} \newenvironment{smatrix}{\left(\begin{smallmatrix}}{\end{smallmatrix}\right)} %SMALL
      %\usepackage{nccmath} \newenvironment{smatrix}{\left(\begin{mmatrix}}{\end{mmatrix}\right)} %MEDIUM
      \begin{equation}
      \begin{smatrix} %explained vector
      \Delta X_{t}^{1} \\ \Delta X_{t}^{2}
      \end{smatrix}=
      \left\{
      \begin{array}{ll}
      \begin{smatrix} %ECT
      -0.0030(0.0060) \\ 0.1749(0.0718)
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %trend
      6.1\text{e-06}(1\text{e-05}) \\ 0.0008(0.0001)
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      0.1616(0.0694) & -0.0076(0.0056) \\
      0.5407(0.8276) & 0.3552(0.0664) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      0.0466(0.0691) & 0.0041(0.0056) \\
      -1.3326(0.8238) & -0.0463(0.0662) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}< 0.1839 \\
      \begin{smatrix} %ECT
      -0.0599(0.0983) \\ 4.2280(1.1725)
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %trend
      6.6\text{e-05}(0.0001) \\ -0.0047(0.0015)
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      0.1752(0.2472) & 0.0419(0.0149) \\
      -1.3831(2.9475) & 0.7378(0.1772) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      0.2040(0.2533) & -0.0341(0.0158) \\
      -2.4316(3.0205) & -0.3810(0.1884) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if } 0.1839 < \text{Th} < 0.2161 \\
      \begin{smatrix} %ECT
      0.0068(0.0059) \\ 0.1360(0.0702)
      \end{smatrix}ECT_{-1}+
      \begin{smatrix}     %trend
      6.5\text{e-05}(5.5\text{e-05}) \\ 0.0014(0.0007)
      \end{smatrix}     %trend
      +\begin{smatrix}      %Lag1
      0.1067(0.1146) & -0.0028(0.0114) \\
      -1.2594(1.3661) & 0.3674(0.1360) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-1}^{1} \\ \Delta X_{t-1}^{2}
      \end{smatrix}
      +\begin{smatrix}      %Lag2
      -0.2407(0.1214) & -0.0101(0.0116) \\
      0.6007(1.4470) & 0.2467(0.1389) 
      \end{smatrix}
      \begin{smatrix}
      \Delta X_{t-2}^{1} \\ \Delta X_{t-2}^{2}
      \end{smatrix}
      & \text{if Th}> 0.2161 \\
      \end{array}
      \right.
      \end{equation}
      
