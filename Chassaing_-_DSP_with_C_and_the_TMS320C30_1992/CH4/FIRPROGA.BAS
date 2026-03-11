10 REM                  *** Kaiser Window Program Module ***
20 REM
30 REM   This program module will generate the necessary FIR coefficients
40 REM   using the Kaiser window sequence.  The user can then define the
50 REM   output path to either the screen, printer, or an external file
60 REM   which can be merged with a FIR program to implement the filter.
70 REM
80 DIM C(256),H(256),FACT(256),WK(256),T(256),CHEX$(256),CARY(256)
90 PI = 3.1415927#
100 CLS
110 REM
120 REM                  *** Generate filter type menu ***
130 REM
140 LOCATE 4
150 PRINT TAB(11);"*** FIR COEFFICIENT GENERATION USING THE KAISER WINDOW ***"
160 LOCATE 8
170 PRINT TAB(22);"Selections:":PRINT
180 PRINT TAB(33);"1....LOWPASS":PRINT
190 PRINT TAB(33);"2....HIGHPASS":PRINT
200 PRINT TAB(33);"3....BANDPASS":PRINT
210 PRINT TAB(33);"4....BANDSTOP":PRINT
220 REM
230 PRINT TAB(33);"5....Exit back to Main Menu":PRINT
240 INPUT "        Enter desired filter type (number only) --> ",TYPE
250 IF TYPE = 5 THEN LOAD "firprog.bas",R
260 IF TYPE = 1 THEN GOSUB 1290:GOTO 300          'Lowpass prompts routine
270 IF TYPE = 2 THEN GOSUB 1510:GOTO 300          'Highpass prompts routine
280 IF TYPE=3 OR TYPE=4 THEN GOSUB 1730:GOTO 300  'bandpass/stop prompts routine
290 GOTO 100
300 REM
310 REM      *** Prompt for general information and output specifications ***
320 REM
330 MES=0               'flag to which error has occured
340 LOCATE 12
350 INPUT "        Enter the Sampling Frequency (Fs) in Hz --> ",FSAM
360 IF TYPE=1 THEN FPASS1=0:FSTOP1=0
370 IF TYPE=2 THEN FPASS2=0:FSTOP2=0
380 IF FSAM/2<FPASS1 OR FSAM/2<FPASS2 OR FSAM/2<FSTOP1 OR FSAM/2<FSTOP2 THEN MES=1:GOSUB 4540:GOTO 330    'display error message
390 REM
400 IF TYPE=1 OR TYPE=2 THEN LOCATE 7
410 IF TYPE=3 OR TYPE=4 THEN LOCATE 9
420 PRINT TAB(24);"Sampling Frequency (Fs) =";FSAM;"Hz"
430 LOCATE 12
440 PRINT "                                                                    "
450 LOCATE 12
460 INPUT "        Are the above specifications correct (y/n)? --> ",RES$
470 IF RES$="n" OR RES$="N" THEN 100
480 LOCATE 12
490 PRINT "                                                                    "
500 REM                   *** Calculate the Order required ***
510 FCUT2=ABS(FPASS2-FSTOP2)
520 FCUT1=ABS(FPASS1-FSTOP1)
530 IF TYPE=1 THEN DELTAF=FCUT2
540 IF TYPE=2 THEN DELTAF=FCUT1
550 IF TYPE <> 3 AND TYPE <> 4 THEN 570
560 IF FCUT2 < FCUT1 THEN DELTAF = FCUT2 ELSE DELTAF = FCUT1
570 D2=10^(-.05*ATT)
580 D1=(10^(.05*RIP)-1)/(10^(.05*RIP)+1)
590 IF D1 < D2 THEN NUM = D1 ELSE NUM = D2
600 AHP=-20*LOG(NUM)/LOG(10)
610 IF AHP <= 21 THEN D = .9222 ELSE D = (AHP-7.95)/14.36
620 COEFF=INT(2+D*FSAM/DELTAF)
630 IF COEFF/2 = INT(COEFF/2) THEN COEFF=COEFF+1
640 Q=CINT(COEFF-1)/2
650 LOCATE 12
660 PRINT "        The calculated # of coefficients required is:";COEFF
670 PRINT
680 PRINT "        Enter # of coefficients desired ONLY if greater than";COEFF
690 INPUT "        otherwise, press <Enter> to continue --> ",TEMP
700 IF TEMP = 0 THEN 830
710 IF TEMP < COEFF THEN 750
720 COEFF=TEMP
730 Q=(COEFF-1)/2
740 GOTO 830
750 PRINT
760 FOR BLINK=1 TO 10
770 LOCATE 20
780 PRINT TAB(12);"ERROR! : Order will not satisfy specifications - reenter"
790 FOR DELAY=1 TO 100:NEXT DELAY
800 LOCATE 20:PRINT "                                                                               "
810 NEXT BLINK
820 GOTO 650
830 CLS
840 LOCATE 12
850 PRINT TAB(28);"Please wait ...working"
860 REM                     *** Compute Coefficients ***
870 REM
880 IF AHP <= 21 THEN ALP = 0 ELSE ALP = (.1102*(AHP-8.7))
890 IF AHP > 21 AND AHP <= 50 THEN ALP=.5842*(AHP-21)^.4+.07886*(AHP-21)
900 FACT(1)=1
910 FOR I=2 TO 30
920 FACT[I] = FACT(I-1)*I
930 NEXT I
940 FOR I=0 TO Q
950 BESS=ALP*SQR(1-(2*I/(COEFF-1))^2)
960 IOBES=1:IOALP=1
970 FOR N = 1 TO 30
980 IOBES = IOBES+(((BESS/2)^(N))/FACT(N))^2
990 IOALP = IOALP+(((ALP/2)^(N))/FACT(N))^2
1000 NEXT N
1010 WK(I) = IOBES/IOALP
1020 NEXT I
1030 REM
1040 IF TYPE = 1 THEN GOSUB 2100:GOTO 1080        'calculate coeffs for lowpass
1050 IF TYPE = 2 THEN GOSUB 2210:GOTO 1080        'calculate coeffs for highpass
1060 IF TYPE = 3 THEN GOSUB 2320:GOTO 1080        'calculate coeffs for bandpass
1070 IF TYPE = 4 THEN GOSUB 2430                 'calculate coeffs for bandstop
1080 REM
1090 REM                       *** Rearrange Coefficients ***
1100 FOR I = 0 TO COEFF
1110 C(I)=H(I)         ' C'(i)
1120 T(I)=H(I)        ' T(i) is temporary storage for H(i)'s
1130 NEXT I
1140 FOR I=0 TO Q
1150 T(I)=H(Q-I)
1160 NEXT I
1170 FOR I=1 TO Q
1180 T(Q+I)=T(Q-I)
1190 NEXT I
1200 FOR I=0 TO COEFF
1210 H(I)=T(I)
1220 NEXT I
1230 REM                *** Convert and Output Coefficients ***
1240 REM
1250 GOSUB 4390         'convert coeffs to hex
1260 REM
1270 GOSUB 2760         'call output routine
1280 REM
1290 REM ======================= Lowpass Prompts Routine =======================
1300 TYPE$="LOWPASS"
1310 CLS
1320 REM
1330 GOSUB 2690         'output filter type
1340 REM
1350 GOSUB 2540         'prompt for AS and AP
1360 REM
1370 LOCATE 12
1380 INPUT "        Enter the passband frequency in Hz --> ",FPASS2
1390 LOCATE 12
1400 PRINT "                                                                   "
1410 LOCATE 5:PRINT TAB(24);"Passband Frequency =";FPASS2;"Hz"
1420 LOCATE 12
1430 INPUT "        Enter the stopband frequency in Hz --> ",FSTOP2
1440 IF FSTOP2 <= FPASS2 THEN GOSUB 4540:GOTO 1420       'display error message
1450 LOCATE 12
1460 PRINT "                                                                   "
1470 LOCATE 6:PRINT TAB(24);"Stopband Frequency =";FSTOP2;"Hz"
1480 RETURN
1490 REM =======================================================================
1500 REM
1510 REM ====================== Highpass Prompts Routine =======================
1520 TYPE$="HIGHPASS"
1530 CLS
1540 REM
1550 GOSUB 2710         'output filter type
1560 REM
1570 GOSUB 2540         'prompt for AS and AP
1580 REM
1590 LOCATE 12
1600 INPUT "        Enter the passband frequency in Hz --> ",FPASS1
1610 LOCATE 12
1620 PRINT "                                                                   "
1630 LOCATE 5:PRINT TAB(24);"Passband Frequency =";FPASS1;"Hz"
1640 LOCATE 12
1650 INPUT "        Enter the stopband frequency in Hz --> ",FSTOP1
1660 IF FSTOP1 >= FPASS1 THEN GOSUB 4540:GOTO 1640       'display error message
1670 LOCATE 12
1680 PRINT "                                                                   "
1690 LOCATE 6:PRINT TAB(24);"Stopband Frequency =";FSTOP1;"Hz"
1700 RETURN
1710 REM =======================================================================
1720 REM
1730 REM =================== Bandpass/stop Prompts routine =====================
1740 REM
1750 IF TYPE=3 THEN TYPE$="BANDPASS" ELSE TYPE$="BANDSTOP"
1760 CLS
1770 REM
1780 GOSUB 2710         'output filter type
1790 REM
1800 GOSUB 2540         'prompt for AS and AP
1810 REM
1820 LOCATE 12
1830 INPUT "        Enter the lower passband frequency in Hz --> ",FPASS1
1840 LOCATE 12
1850 PRINT "                                                                   "
1860 LOCATE 5:PRINT TAB(24);"Lower Passband Frequency =";FPASS1;"Hz"
1870 LOCATE 12
1880 INPUT "        Enter the upper passband frequency in Hz --> ",FPASS2
1890 IF FPASS1 >= FPASS2 THEN GOSUB 4540:GOTO 1870    'display error message
1900 LOCATE 12
1910 PRINT "                                                                   "
1920 LOCATE 6:PRINT TAB(24);"Upper Passband Frequency =";FPASS2;"Hz"
1930 LOCATE 12
1940 INPUT "        Enter the lower stopband frequency in Hz --> ",FSTOP1
1950 IF TYPE=3 THEN IF FSTOP1>=FPASS1 OR FSTOP1>=FPASS2 THEN GOSUB 4540 : GOTO 1930  'display error message
1960 IF TYPE=4 THEN IF FSTOP1<=FPASS1 OR FSTOP1>=FPASS2 THEN GOSUB 4540 : GOTO 1930  'display error message
1970 LOCATE 12
1980 PRINT "                                                                   "
1990 LOCATE 7:PRINT TAB(24);"Lower Stopband Frequency =";FSTOP1;"Hz"
2000 LOCATE 12
2010 INPUT "        Enter the upper stopband frequency in Hz --> ",FSTOP2
2020 IF TYPE=3 THEN IF FSTOP2 <= FPASS2 OR FSTOP2 <= FPASS1 OR FSTOP2 <= FSTOP1 THEN GOSUB 4540:GOTO 2000     'display error message
2030 IF TYPE=4 THEN IF FSTOP2 >= FPASS2 OR FSTOP2 <= FPASS1 OR FSTOP2 <= FSTOP1 THEN GOSUB 4540:GOTO 2000     'display error message
2040 LOCATE 12
2050 PRINT "                                                                   "
2060 LOCATE 8:PRINT TAB(24);"Upper Stopband Frequency =";FSTOP2;"Hz"
2070 RETURN
2080 REM =======================================================================
2090 REM
2100 REM ======================= Kaiser Lowpass Routine ========================
2110 REM
2120 LWCUT=(FPASS2+FSTOP2)/2
2130 H0=2*LWCUT/FSAM
2140 H(0)=H0*WK(0)          'H(0) = C(0)*W(0)
2150 FOR I=1 TO Q
2160 H(I)=H0*((SIN(2*LWCUT*I*PI/FSAM))/(2*LWCUT*I*PI/FSAM))*WK(I)
2170 NEXT I
2180 RETURN
2190 REM =======================================================================
2200 REM
2210 REM ===================== Kaiser Highpass Routine =========================
2220 REM
2230 HWCUT=(FPASS1+FSTOP1)/2
2240 H0=-2*HWCUT/FSAM
2250 H(0)=(H0+1)*WK(0)          'H(0) = C(0)*W(0)
2260 FOR I=1 TO Q
2270 H(I)=H0*((SIN(2*HWCUT*I*PI/FSAM))/(2*HWCUT*I*PI/FSAM))*WK(I)
2280 NEXT I
2290 RETURN
2300 REM =======================================================================
2310 REM
2320 REM ====================== Kaiser Bandpass Routine ========================
2330 REM
2340 FCUT1=FPASS1-DELTAF/2
2350 FCUT2=FPASS2+DELTAF/2
2360 H(0)=(2/FSAM)*(FCUT2-FCUT1)*WK(0)
2370 FOR I=1 TO Q
2380 H(I)=(1/(I*PI))*(SIN((2*PI*I*FCUT2)/FSAM)-(SIN((2*PI*I*FCUT1)/FSAM)))*WK(I)
2390 NEXT I
2400 RETURN
2410 REM =======================================================================
2420 REM
2430 REM ====================== Kaiser Bandstop Routine ========================
2440 REM
2450 FCUT1=FPASS1+DELTAF/2
2460 FCUT2=FPASS2-DELTAF/2
2470 H(0)=((2*(FCUT1-FCUT2)/FSAM)+1)*WK(0)
2480 FOR I=1 TO Q
2490 H(I)=(1/(I*PI))*(SIN((2*PI*I*FCUT1)/FSAM)-(SIN((2*PI*I*FCUT2)/FSAM)))*WK(I)
2500 NEXT I
2510 RETURN
2520 REM =======================================================================
2530 REM
2540 REM ==================== Prompt for AS and AP Routine =====================
2550 REM
2560 LOCATE 12
2570 INPUT "        Enter the passband ripple in db --> ",RIP
2580 LOCATE 12
2590 PRINT "                                                                   "
2600 LOCATE 3:PRINT TAB(24);"Passband Ripple (AP) =";RIP;"db"
2610 LOCATE 12
2620 INPUT "        Enter the minimum stopband attenuation in db --> ",ATT
2630 LOCATE 12
2640 PRINT "                                                                   "
2650 LOCATE 4:PRINT TAB(24);"Stopband Attenuation (AS) =";ATT;"db"
2660 RETURN
2670 REM =======================================================================
2680 REM
2690 REM ======================= Output Specs Header ===========================
2700 LOCATE 1
2710 PRINT "        Specifications:"
2720 PRINT TAB(24);TYPE$
2730 RETURN
2740 REM =======================================================================
2750 REM
2760 REM ========================== Output Routine =============================
2770 REM
2780 REM        This routine allows the user to define the output
2790 REM        path.  These include the terminal, line printer, or
2800 REM        data file, which includes TMS320 data format for simple
2810 REM        merging into a 'blank' filter program.
2820 REM
2830 CLS
2840 LOCATE 8
2850 PRINT "       Send coefficients to:"
2860 PRINT TAB(29);"(S)creen"
2870 PRINT TAB(29);"(P)rinter"
2880 PRINT TAB(29);"(F)ile: contains TMS320 (C25 or C30) data format"
2890 PRINT TAB(29);"(R)eturn to Filter Type Menu"
2900 PRINT TAB(29);"(E)xit to DOS"
2910 PRINT
2920 INPUT "       Enter desired path --> ",PATH$
2930 IF PATH$ = "S" OR PATH$ = "s" THEN GOSUB 3000
2940 IF PATH$ = "P" OR PATH$ = "p" THEN GOSUB 3320
2950 IF PATH$ = "F" OR PATH$ = "f" THEN GOSUB 3800
2960 IF PATH$ = "R" OR PATH$ = "r" THEN GOTO 100
2970 IF PATH$ = "E" OR PATH$ = "e" THEN CLS:SYSTEM
2980 GOTO 2830
2990 REM
3000 REM -------------------- Output Coefficients to Terminal ------------------
3010 REM
3020 CLS
3030 PRINT TAB(15);COEFF;"Coefficient ";TYPE$;" Filter Using the KAISER Window"
3040 PRINT
3050 PRINT "       Filter Coefficients:"
3060 PRINT
3070 PRINT TAB(20);"C'(n)";TAB(53);"H(n)"
3080 PRINT TAB(20);"-----";TAB(53);"----"
3090 PRINT
3100 PRINT TAB(51);"DECIMAL";TAB(70);"HEX"
3110 PRINT
3120 FOR I=0 TO Q
3130 IF INT(I/11) <> I/11 OR I = 0 THEN 3190
3140 LOCATE 23:PRINT "       Press <Enter> to view";Q-I+1;"remaining coefficients";:INPUT " ",RESPONSE
3150 LOCATE 23
3160 PRINT "                                                                   "
3170 GOSUB 4090    'erase coefficients from screen
3180 REM
3190 GOSUB 4190    'convert coeff. index to char, removing spaces in output
3200 PRINT TAB(8);
3210 PRINT "C'(";C1$;")";
3220 PRINT TAB(15);" = ";
3230 PRINT USING "+#.#####";C(I);:PRINT " = C'(";C2$;")";
3240 PRINT TAB(41);
3250 PRINT "H(";C1$;")";
3260 PRINT TAB(47);" = ";
3270 PRINT USING "+#.#####";H(I);: PRINT " = H(";C2$;")";
3280 PRINT TAB(67);" = ";CHEX$(I)
3290 NEXT I
3300 LOCATE 23:INPUT "       Press <Enter> to return to Output Path Menu ",RET
3310 RETURN
3320 REM --------------------- Output Coefficients to Printer ------------------
3330 REM
3340 REM        This nested routine allows the user to produce hardcopy of the
3350 REM        coefficients.  NOTE: If this option is specified from the
3360 REM        output path menu and the printer is NOT READY, then a
3370 REM        "device timeout" will occur via GWBASIC.  This will cause
3380 REM        exiting to GWBASIC.  To restart the FILTER DEVELOPMENT
3390 REM        PACKAGE type <F2> <Enter>.  This will restart whichever
3400 REM        window module selected previously from Main Menu.  Program
3410 REM        must be re-run to re-develop filter to obtain the hardcopy
3420 REM        initially desired.
3430 REM
3440 LPRINT TAB(26);"*** FIR DEVELOPMENT PACKAGE ***"
3450 LPRINT
3460 LPRINT TAB(15);COEFF;"Coefficient ";TYPE$;" Filter Using the KAISER Window"
3470 LPRINT
3480 LPRINT TAB(8);"Passband Ripple (AP) =";RIP;"db"
3490 LPRINT TAB(8);"Stopband Attenuation (AS) =";ATT;"db"
3500 IF TYPE<>1 THEN 3560
3510 LPRINT TAB(8);"Passband frequency =";FPASS2;"Hz"
3520 LPRINT TAB(8);"Stopband frequency =";FSTOP2;"Hz"
3530 IF TYPE<>2 THEN 3560
3540 LPRINT TAB(8);"Passband frequency =";FPASS1;"Hz"
3550 LPRINT TAB(8);"Stopband frequency =";FSTOP1;"Hz"
3560 LPRINT TAB(8);"Lower Passband Frequency =";FPASS1;"Hz"
3570 LPRINT TAB(8);"Upper Passband Frequency =";FPASS2;"Hz"
3580 LPRINT TAB(8);"Lower Stopband Frequency =";FSTOP1;"Hz"
3590 LPRINT TAB(8);"Upper Stopband Frequency =";FSTOP2;"Hz"
3600 LPRINT TAB(8);"Sampling Frequency (Fs) =";FSAM;"Hz"
3610 LPRINT
3620 LPRINT TAB(20);"C'(n)";TAB(53);"H(n)"
3630 LPRINT TAB(20);"-----";TAB(53);"----"
3640 LPRINT
3650 LPRINT TAB(51);"DECIMAL";TAB(70);"HEX"
3660 LPRINT
3670 FOR I=0 TO Q
3680 GOSUB 4190    'convert coeff. index to char, removing spaces in output
3690 LPRINT TAB(8);
3700 LPRINT "C'(";C1$;")";
3710 LPRINT TAB(15);" = ";
3720 LPRINT USING "+#.#####";C(I);:LPRINT " = C'(";C2$;")";
3730 LPRINT TAB(41);
3740 LPRINT "H(";C1$;")";
3750 LPRINT TAB(47);" = ";
3760 LPRINT USING "+#.#####";H(I);:LPRINT " = H(";C2$;")";
3770 LPRINT TAB(67);" = ";CHEX$(I)
3780 NEXT I
3790 RETURN
3800 REM -------------------- Output Coefficients to File ----------------------
3810 REM
3820 REM        This routine writes the coefficients to a file (named by the
3830 REM        user).  The file contains TMS320 (C30 or C25) data format.
3840 REM
3850 REM        This file can then be merged into a 'blank' filter program
3860 REM        via a word processor like PC-WRITE.  NOTE: A directory or
3870 REM        drive specification can be given along with the filename.
3875 LOCATE 17
3876 PRINT "                                                                   "
3877 LOCATE 17
3879 INPUT "       Enter DSP type (C25 OR C30):";DSP$
3880 IF DSP$="C30" OR DSP$="c30" THEN GOTO 4030
3885 IF DSP$="C25" OR DSP$="c25" THEN GOTO 3900
3890 GOTO 3875
3900 CLS
3910 LOCATE 12
3920 INPUT "       Enter Filename : ",NAM$
3930 LOCATE 12
3935 PRINT "                                                                   "
3940 LOCATE 12:PRINT TAB(30);"...saving ";NAM$
3950 OPEN NAM$ FOR OUTPUT AS 1
3960 FOR I=0 TO COEFF-1
3970 S$=STR$(I)
3980 IF I>9 AND I<100 THEN PRINT#1,"H";RIGHT$(S$,2);"     DATA    >";CHEX$(I) : GOTO 4010
3990 IF I>99 THEN PRINT#1,"H";RIGHT$(S$,3);"    DATA    >";CHEX$(I):GOTO 4010
4000 PRINT#1,"H";RIGHT$(S$,1);"      DATA    >";CHEX$(I)
4010 NEXT I
4020 CLOSE 1
4025 RETURN
4030 CLS
4032 LOCATE 12
4035 INPUT "       Enter Filename : ",NAM$
4038 LOCATE 12
4040 PRINT "                                                                   "
4042 LOCATE 12:PRINT TAB(30);"...saving ";NAM$
4045 OPEN NAM$ FOR OUTPUT AS 1
4048 PRINT#1, USING "COEFF   .FLOAT  +##.########       ;H(###)  ==>  H(N-1)";H(COEFF-1);COEFF-1
4050 FOR I=COEFF-2 TO 1 STEP -1
4052 PRINT#1, USING "        .FLOAT  +##.########       ;H(###)";H(I);I
4055 NEXT I
4058 PRINT#1, USING "H0      .FLOAT  +##.########       ;H(  0)";H(0)
4060 PRINT#1, "LENGTH  .SET (H0-COEFF)+1          ;LENGTH = # OF COEFFs =";COEFF
4062 CLOSE 1
4065 RETURN
4070 REM =======================================================================
4090 REM ======================== Clear Coefficients ===========================
4100 REM
4110 LOCATE 10
4120 FOR N=1 TO 12
4130 PRINT "                                                                               "
4140 NEXT N
4150 LOCATE 10
4160 RETURN
4170 REM =======================================================================
4180 REM
4190 REM ======================= Integer to Character ==========================
4200 REM
4210 REM        This routine eliminates the spaces produced in the output
4220 REM        when an Integer TYPE is printed.  The spaces before and after
4230 REM        the integer value are set aside in the event that the integer
4240 REM        value is negative.  Since the Coefficient index is always
4250 REM        non-negative, this value is converter to a character type,
4260 REM        and output from this routine as c1$ and c2$.
4270 REM
4280 S1$=STR$(I)
4290 S2$=STR$(COEFF-I-1)      'COEFF-I-1 is the symmetrical alternate of index
4300 C1$=RIGHT$(S1$,1)
4310 C2$=RIGHT$(S2$,1)
4320 IF I > 9 THEN C1$=RIGHT$(S1$,2)
4330 IF I > 99 THEN C1$=RIGHT$(S1$,3)
4340 IF COEFF-I-1 > 9 THEN C2$=RIGHT$(S2$,2)
4350 IF COEFF-I-1 > 99 THEN C2$=RIGHT$(S2$,3)
4360 RETURN
4370 REM =======================================================================
4380 REM
4390 REM ==================== Convert Coefficients to Hex ======================
4400 REM
4410 REM        This routine converts the coefficients to Hexidecimal
4420 REM        via a routine in GWBASIC called HEX$().  The HEX values
4430 REM        are accurate to + or - 1.
4440 REM
4450 FOR M=0 TO COEFF-1
4460 CARY(M)=CINT(H(M)*32768!)    'scale by 2^15
4470 IF CARY(M) >= 0 GOTO 4490
4480 CARY(M)=65536!+CARY(M)
4490 CHEX$(M)=HEX$(CARY(M))
4500 NEXT M
4510 RETURN
4520 REM =======================================================================
4530 REM
4540 REM ====================== Error Message Routine ==========================
4550 REM
4560 FOR BLINK=1 TO 10
4570 IF MES<>1 THEN 4610
4580 LOCATE 20
4590 PRINT TAB(14);"ERROR : Sampling Frequency (Fs) >= 2*Nyquist - reenter"
4600 GOTO 4630
4610 LOCATE 20
4620 PRINT TAB(15);"ERROR! : Frequency value is inconsistant - reenter"
4630 FOR DELAY=1 TO 100:NEXT DELAY
4640 LOCATE 20
4650 PRINT "                                                                   "
4660 NEXT BLINK
4670 RETURN
4680 REM =======================================================================
