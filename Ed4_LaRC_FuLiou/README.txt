README Ed4_LaRC_FuLiou          Jan 2015  fred.g.rose@nasa.gov

Source code requires FORTRAN to compile. 

!******************************************************************************!
!  Please do not re-distribute the code. Instead please direct other potential ! 
!  users to the CAVE website so that we might be able to inform all who have   !
!  an interest about changes and or errors in the source code.                 !
!             http://www-cave.larc.nasa.gov/cgi-bin/fuliou/runfl.cgi           !
!******************************************************************************!

Once you have downloaded the Fuliou distribution..

1) Uncompress using      : gunzip Ed4_LaRC_FuLiou20150106.tar.gz
2) Untar using           : tar -xvf Ed4_LaRC_FuLiou20150106.tar
3) Change directory into : Ed4_LaRC_FuLiou
4) Set up local environment variables to compile F90 and F77 source code.

An "example" on a unix computer with a gfortran compiler would be..

 setenv F90COMP " -O2 -c "
 setenv FCOMP " -O2 -c"
 setenv F90 /usr/local/bin/gfortran
 setenv F77 /usr/local/bin/gfortran

5) Familiarize yourself with the directory structure
 
 The radiative transfer F90 source code is under ./lib/src
 Once the make file is run an object library is created as ./lib/libEd3Fu_201212.a
 F90 .mod files are created under ./lib/mod
 An example code to show how to setup inputs is under ./src/simple/simple.f90 
 Example inputs of standard atmosphere temperature and humidity ./testatms

6)From the top level directory type:  make 

This should 
  A) Compile the code library under ./lib/src
  B) Compile the example code .src/simple.f90
  C) Execute the example code executible ./src/simple

\
 > comment: compliation should produce something similar to the following:
/

>make

cd ./lib/src ;make
/usr/local/bin/gfortran  -O2 -c     extras.f90
/usr/local/bin/gfortran  -O2 -c     taucorr.f90
/usr/local/bin/gfortran  -O2 -c     fuinput.f90
/usr/local/bin/gfortran  -O2 -c     fuoutput.f90
/usr/local/bin/gfortran  -O2 -c     fuprint.f90
/usr/local/bin/gfortran  -O2 -c     entropy_lw.f90
/usr/local/bin/gfortran  -O2 -c     icedirsfc.f90
/usr/local/bin/gfortran  -O2 -c     ma_tip.f90
/usr/local/bin/gfortran  -O2 -c     qsortd.f90
/usr/local/bin/gfortran  -O2 -c     vla.f90
/usr/local/bin/gfortran  -O2 -c     uvcor_all.f90
/usr/local/bin/gfortran  -O2 -c     rad_multi_200511.f90
/usr/local/bin/gfortran  -O2 -c     seiji_k2d.f90
/usr/local/bin/gfortran  -O2 -c     seiji_solver_200511.f90
/usr/local/bin/gfortran  -O2 -c     gflq.f90
/usr/local/bin/gfortran  -O2 -c     calipso_output.f90
/usr/local/bin/gfortran  -O2 -c     zjin.f90
/usr/local/bin/gfortran  -O2 -c     ar_asy.f90
/usr/local/bin/gfortran  -O2 -c    aqua_wnflt_0404.f
/usr/local/bin/gfortran  -O2 -c    misc_200511.f
/usr/local/bin/gfortran  -O2 -c    cloud_optics.f
/usr/local/bin/gfortran  -O2 -c    seiji_twostreamsolv_sw_v20.f
seiji_twostreamsolv_sw_v20.f:168.48:

     &      af_clear,bf_clear,ef_clear,ak_clear,u1i,u1s,
                                                1
Warning: Type mismatch in argument 'u1i' at (1); passed REAL(8) to REAL(4)
/usr/local/bin/gfortran  -O2 -c    aerosols_200511.f
ar -rcv libEd3Fu_201212.a  extras.o fuinput.o fuoutput.o fuprint.o entropy_lw.o icedirsfc.o ma_tip.o qsortd.o vla.o rad_multi_200511.o seiji_k2d.o seiji_solver_200511.o taucorr.o uvcor_all.o gflq.o calipso_output.o zjin.o ar_asy.o aqua_wnflt_0404.o   misc_200511.o cloud_optics.o 	seiji_twostreamsolv_sw_v20.o aerosols_200511.o
a - extras.o
a - fuinput.o
a - fuoutput.o
a - fuprint.o
a - entropy_lw.o
a - icedirsfc.o
a - ma_tip.o
a - qsortd.o
a - vla.o
a - rad_multi_200511.o
a - seiji_k2d.o
a - seiji_solver_200511.o
a - taucorr.o
a - uvcor_all.o
a - gflq.o
a - calipso_output.o
a - zjin.o
a - ar_asy.o
a - aqua_wnflt_0404.o
a - misc_200511.o
a - cloud_optics.o
a - seiji_twostreamsolv_sw_v20.o
a - aerosols_200511.o
\cp *.mod ../mod
\cp libEd3Fu_201212.a  ../
cd ./src/simple ; make
/usr/local/bin/gfortran  -O2 -c  -I ../../lib/mod   simple.f90
/usr/local/bin/gfortran  -o simple simple.o ../../lib/libEd3Fu_201212.a
cd ./src/simple ; ./simple
\
 > comment: end of compliation, below is output from simple.f90
/
 ========================================================================================================================
 Fu-Liou Model inputs in structure fi%  Begin
 # of Model LAYERS      :          34
 Solver Config Modes    : T T T T
 Curved Earth Airmass Co: T
 nirold Ray,Ice,Wat,Gas,Kwc : F F F F F
 Solar Constant (wm-2)  :   1365.0000
 Cosine Solar Zenith    :   1.0000000
 Cosine View Zenith     :  0.80000001
 fu%txt                 :           0
      Spect Emissivity 	: 0.990 0.990 0.990 0.990 0.990 0.990 0.990 0.990 0.990 0.990 0.990 0.990
 Skin Temperture (k)    :   294.00000
 Trace Gas Concentration__________________________________________________________
 CO2 Conc (ppmv)        :   360.00000
 CH4 Conc (ppmv)        :   1.7500000
 N2O Conc (ppmv)        :  0.31000000
 CFCs Conc (ppv)        :  2.68000011E-10  5.02999975E-10  1.05000002E-10
 Option Selection_________________________________________________________________
 >4 micron solar lband6a:     T
   Continuum option sel	:     5
 # of LW bands >2200cm-1:     2
 Hybrid solver option   :     T
 Solver option          :     1
      Window instrument	:     0
 Fourstream Sol fourssl :     F
 Fourstream IR  foursir :     F
 Cloud lwc profile flag :     2
 Aerosols__________________________________________________________________________
 #Aerosol Taus		:           1
 #Aerosol Constituents	:           2
Aer.Wavelength(s)(micron)   0.641
 -Aerosol Type 		:           2
Aer. Optical Depth(s)     0.80000
 -Aerosol Type 		:           1
Aer. Optical Depth(s)     0.20000
 Profiles__________________________________________________________________________
 Level.Pres(hPa).Temp(K).H20(g/g).RH(%)..O3(g/g)...AOT%PROFILES
     1     0.10   226.21 2.60E-06   0.0 1.20E-06  0.00  0.00
     2     0.21   243.34 2.79E-06  -0.0 2.02E-06  0.00  0.00
     3     0.47   260.47 2.99E-06  -0.0 2.84E-06  0.00  0.00
     4     0.87   274.01 3.14E-06  -0.0 3.49E-06  0.00  0.00
     5     1.63   270.72 3.10E-06  -0.0 5.47E-06  0.00  0.00
     6     3.17   258.93 2.97E-06   0.0 8.86E-06  0.00  0.00
     7     6.37   245.46 2.82E-06   0.0 9.90E-06  0.00  0.00
     8    13.21   233.99 2.69E-06   0.0 1.02E-05  0.00  0.00
     9    19.03   229.07 3.33E-06   0.1 8.59E-06  0.00  0.00
    10    27.86   223.96 3.99E-06   0.4 6.94E-06  0.01  0.01
    11    41.08   220.82 3.98E-06   0.9 5.43E-06  0.04  0.04
    12    60.79   217.86 3.98E-06   1.9 3.48E-06  0.03  0.03
    13    70.00   216.95 3.98E-06   2.5 2.83E-06  0.08  0.08
    14    90.27   216.00 4.00E-06   3.7 1.75E-06  0.30  0.30
    15   134.06   216.00 4.02E-06   5.5 8.73E-07  0.81  0.81
    16   196.98   219.71 1.44E-05  17.8 4.25E-07  0.05  0.05
    17   200.00   220.30 1.57E-05  18.4 4.10E-07  2.02  2.02
    18   283.31   235.40 1.60E-04  44.2 2.14E-07  1.53  1.53
    19   326.26   242.30 2.60E-04  39.3 1.83E-07  2.16  2.16
    20   374.44   248.34 4.13E-04  38.6 1.50E-07  1.36  1.36
    21   400.00   251.75 5.28E-04  37.6 1.39E-07  1.67  1.67
    22   428.11   255.22 6.49E-04  29.7 1.28E-07  4.20  4.20
    23   487.95   261.09 9.44E-04  30.3 1.06E-07  0.96  0.96
    24   500.00   262.23 1.03E-03  30.9 1.03E-07  4.84  4.84
    25   554.50   267.04 1.39E-03  31.8 9.12E-08  7.94  7.94
    26   628.32   273.02 2.36E-03  39.1 7.98E-08 10.78 10.78
    27   709.97   279.00 3.87E-03  47.2 6.99E-08  6.86  6.86
    28   754.70   282.01 4.92E-03  51.8 6.56E-08  7.88  7.88
    29   801.14   284.95 5.95E-03  54.5 6.13E-08  9.15  9.15
    30   850.00   287.47 7.26E-03  59.7 5.83E-08 10.71 10.71
    31   901.78   289.99 8.58E-03  63.5 5.54E-08 12.41 12.41
    32   956.21   292.01 1.01E-02  70.0 5.27E-08  8.30  8.30
    33   989.89   293.20 1.11E-02  73.3 5.11E-08  2.94  2.94
    34  1001.38   293.60 1.14E-02  74.3 5.05E-08  2.97  2.97
    35  1012.76   293.99 1.17E-02  75.3 5.00E-08  0.00  0.00
 Spectral Surface Albedo WITH AEROSOLS::
    Spect Surface albedo:w/AOT CLEAR  0    0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    Spect Surface albedo:w/AOT Cloud  1    0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
 Spectral Surface Albedo WITHOUT Aerosol::
    Spect Surface albedo:NOAOT CLEAR  0    0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
    Spect Surface albedo:NOAOT Cloud  1    0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
================================================================================
  CLOUDS::    1       1
    Fractions :   1.000
    DPI   mode:       F
    #Overlaps :       1
-------------------------
    Opt Depth :  10.000
    RPhase    :    2.00
    Re        :   15.00
    De        :   60.00
  Top:Bot Lay#:  17  20
  Top:Bot Pres: 200 400
    Nu .......:     0.0
    Mn_lin_tau:    11.5
 Fu-Liou Model inputs in structure fi%  End
 ========================================================================================================================
                      SHORTWAVE Down---------------    Shortwave Up------------------
  #  Presure  Height  Clear   Prist   Total  TotNOA    Clear   Prist   Total  TotNOA
 Lev  [hPa] [meters]   Down    Down    Down    Down     Up      Up      Up      UP
   1    0.10  66295. 1365.03 1365.03 1365.03 1365.03  139.24   51.95  590.23  561.83
   2   70.00  18904. 1330.89 1330.14 1331.54 1331.35  140.04   50.45  597.33  568.19
   3  200.00  12241. 1316.23 1315.08 1320.75 1320.33  136.55   45.08  600.64  569.63
   4  500.00   5780. 1226.03 1237.43  611.71  580.64  119.38   29.94   88.69   16.12
   5  850.00   1502. 1036.95 1117.80  497.63  527.49   62.56   10.01   38.19    4.67
   6 1012.76      2.  909.80 1063.17  426.23  501.17    0.00    0.00    0.00    0.00
                     LONGWAVE Down----------------    Longwave  Up------------------
  #  Presure  Height  Clear   Prist   Total  TotNOA    Clear   Prist   Total  TotNOA
 Lev  [hPa] [meters]   Down    Down    Down    Down     Up      Up      Up      UP
   1    0.10  66295.    0.00    0.00    0.00    0.00  274.78  278.82  147.73  147.79
   2   70.00  18904.   13.33   13.31   13.33   13.31  273.43  277.56  143.86  143.92
   3  200.00  12241.   27.94   27.68   27.91   27.68  281.18  285.18  146.42  146.43
   4  500.00   5780.  141.41  139.23  238.17  238.00  332.24  334.67  332.66  334.81
   5  850.00   1502.  289.57  283.92  331.75  329.39  400.24  400.64  400.46  400.80
   6 1012.76      2.  356.46  350.81  380.75  377.83  422.72  422.66  422.96  422.93
                     WINDOW Down------------------    WINDOW  Up--------------------
  #  Presure  Height  Clear   Prist   Total  TotNOA    Clear   Prist   Total  TotNOA
 Lev  [hPa] [meters]   Down    Down    Down    Down     Up      Up      Up      UP
   1    0.10  66295.    0.00    0.00    0.00    0.00  101.87  104.72   30.18   30.20
   2   70.00  18904.    1.65    1.65    1.65    1.65  103.45  106.39   29.49   29.51
   3  200.00  12241.    2.03    1.97    2.01    1.97  107.05  109.98   29.90   29.91
   4  500.00   5780.    5.09    4.06   52.37   52.27  111.64  113.57  111.99  113.71
   5  850.00   1502.   30.80   26.15   65.34   63.36  118.92  119.25  119.12  119.41
   6 1012.76      2.   62.49   57.32   84.72   82.04  121.54  121.48  121.76  121.73
STOP  Simple.f90 normal end
\
 > comment: end of output from simple.f90
/

Spectral boundaries and number of "K"s per band for SW and LW parts of code.

SW bands 
   #  mbx#  K's   -----Micron----     -----Cm-1------
    1  1    1    0.1754    0.2247    57000.    44500.  O3
    2  1    1    0.2247    0.2439    44500.    41000.  O3
    3  1    1    0.2439    0.2857    41000.    35000.  O3
    4  1    1    0.2857    0.2985    35000.    33500.  O3
    5  1    1    0.2985    0.3225    33500.    31008.  O3
    6  1    1    0.3225    0.3575    31008.    27972.  O3
    7  1    1    0.3575    0.4375    27972.    22857.  O3
    8  1    1    0.4375    0.4975    22857.    20101.  O3 & H2O
    9  1    1    0.4975    0.5950    20101.    16807.  O3 & H2O
   10  1    1    0.5950    0.6896    16807.    14500.  O3 & H2O
   11  2    8    0.690     0.794     14500     12600   H2O & O2 &O3
   12  3    6    0.794     0.889     12600     11250   H2O 
   13  4    8    0.889     1.042     11250      9600   H2O 
   14  5    7    1.042     1.410      9600      7090   H2O
   15  6    8    1.410     1.9048     7090.     5250.  H2O & CO2
   16  7    7    1.9048    2.5000     5250.     4000.  H2O & CO2 & CH4
   17  8    8    2.5000    3.5088     4000.     2850.  H2O & CO2 & O3 & CH4
   18  9    7    3.5088    4.0000     2850.     2500.  H2O & CO2 & Ch4

LW Bands
   # mbx#  K's   -----Micron----     -----Cm-1------ 
   1   10   2	 4.54  - 5.26  2200 - 1900    H2O
   2   11   3	 5.26  - 5.88  1900 - 1700    H2O
   3   12   4	 5.88  - 7.14  1700 - 1400    H2O
   4   13   4	 7.14  - 8.00  1400 - 1250    H2O& CH4& N2O 
   5   14   3(9) 8.00  - 9.09  1250 - 1100  W H2O& CH4& N2O &Cfc
   6   15   5(8) 9.09  - 10.2  1100 -  980  W H2O & O3      &Cfc
   7   16   2(2) 10.2  - 12.5   980 -  800  W H2O           &Cfc
   8   17   10	 12.5  - 14.9   800 -  670    H2O & CO2
   9   18   12	 14.9  - 18.5   670 -  540    H2O & CO2
  10   19   7	 18.5  - 25.0   540 -  400    H2O
  11   20   7	 25.0  - 35.7   400 -  280    H2O
  12   21   8	 35.7  - ....   280 -    0    H2O

LW "Hidden" bands
  13   22   5     4.54 - 4.0        2200 - 2500 H2O & N2O &CO2
  14   23   5     4.0  - 3.5        2500 - 2850 H2O  &CO2




