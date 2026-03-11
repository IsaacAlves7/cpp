10 REM                    *** FIR DEVELOPMENT PACKAGE ***
20 REM
30 REM   This program module will generate the necessary FIR coefficients
40 REM   using a Rectangular, Hanning, Hamming, and Blackman window
50 REM   sequence.  The user can then define the output path to either
60 REM   the screen, line printer, or an external file which can then
70 REM   be merged with a FIR program to implement the filter.
80 REM
90 REM   NOTE: The Kaiser window sequence is located in a separate module.
100 REM
110 DIM C(256),H(256),CARY(256),CHEX$(256)
120 PI=3.1415927#
130 CLS
140 KEY OFF
150 REM                     *** Generate Main Menu ***
160 REM
170 LOCATE 3:PRINT TAB(27);"FIR DEVELOPMENT PACKAGE"
180 LOCATE 8
190 PRINT TAB(33);"Main Menu"
200 PRINT TAB(33);"---------"
210 PRINT
220 PRINT TAB(33);"1....RECTANGULAR":PRINT
230 PRINT TAB(33);"2....HANNING":PRINT
240 PRINT TAB(33);"3....HAMMING":PRINT
250 PRINT TAB(33);"4....BLACKMAN":PRINT
260 PRINT TAB(33);"5....KAISER":PRINT
270 PRINT TAB(33);"6....Exit to DOS"
280 PRINT
290 INPUT "        Enter window desired (number only) --> ",WIN
300 XPOS=10
310 IF WIN = 6 THEN CLS:SYSTEM
320 IF WIN = 1 THEN WIN$="RECTANGULAR":XPOS=8
330 IF WIN = 2 THEN WIN$="HANNING"
340 IF WIN = 3 THEN WIN$="HAMMING"
350 IF WIN = 4 THEN WIN$="BLACKMAN"
360 IF WIN = 5 THEN LOAD"firproga.bas",R
370 IF WIN < 1 OR WIN > 6 THEN 130
380 REM                  *** Generate filter type menu ***
390 CLS
400 LOCATE 4
410 PRINT TAB(XPOS);"***  FIR COEFFICIENT GENERATION USING THE ";WIN$;" WINDOW ***"
420 LOCATE 8
430 PRINT TAB(22);"Selections:"
440 PRINT
450 PRINT TAB(33);"1....LOWPASS"
460 PRINT
470 PRINT TAB(33);"2....HIGHPASS"
480 PRINT
490 PRINT TAB(33);"3....BANDPASS"
500 PRINT
510 PRINT TAB(33);"4 ...BANDSTOP"
520 PRINT
530 PRINT TAB(33);"5....Exit back to Main Menu"
540 PRINT
550 INPUT "        Enter desired filter type (number only) --> ",TYPE
560 IF TYPE = 5 THEN 130
570 IF TYPE = 1 THEN GOSUB 2000:GOTO 620         'Lowpass Prompts Routine
580 IF TYPE = 2 THEN GOSUB 2120:GOTO 620         'Highpass Prompts Routine
590 IF TYPE=3 OR TYPE=4 THEN GOSUB 2240:GOTO 620 'Bandpass/stop Prompts Routine
600 IF TYPE < 1 OR TYPE > 5 THEN 390
610 GOTO 130
620 REM      *** Prompt for general information and output specifications ***
630 LOCATE 12
640 PRINT "                                                                    "
650 MES=0
660 LOCATE 12
670 INPUT "        Enter the sampling frequency (Fs) in Hz --> ",FSAM
680 IF TYPE=1 THEN LCUT=0
690 IF TYPE=2 THEN HCUT=0
700 IF FSAM/2<LCUT OR FSAM/2<HCUT THEN MES=1:GOSUB 4140:GOTO 650  'display error
710 LOCATE 12
720 PRINT "                                                                    "
730 IF TYPE=1 OR TYPE=2 THEN YPOS=4 ELSE YPOS=5
740 LOCATE (YPOS):PRINT TAB(24);"Sampling Frequency (Fs) =";FSAM;"Hz"
750 LOCATE 14
760 PRINT "        Number of Coefficients = (D*Fs)+1"
770 PRINT
780 LOCATE 12
790 INPUT "        Enter the duration of the impulse response (D) in msec --> ",D
800 LOCATE (YPOS+1):PRINT TAB(24);"Impulse Duration =";D;"msec"
810 LOCATE 12
820 FOR I = 1 TO 3
830 PRINT "                                                                               "
840 NEXT I
850 LOCATE 12:INPUT "        Are the above specifications correct (y/n) ? ",RES$
860 IF RES$="n" OR RES$="N" THEN 390
870 REM
880 REM               *** calculate number of coefficients required ***
890 D = D/1000
900 NYQST=FSAM/2
910 Q=CINT((D*FSAM)/2)
920 COEFF=2*Q+1
930 LOW=LCUT/NYQST      'Nu 1
940 IF TYPE=2 THEN GOTO 960     'if highpass then high=1
950 HIGH=HCUT/NYQST     'Nu 2
960 LOCATE 12
970 PRINT "                                                                    "
980 LOCATE 12
990 PRINT "        The calculated # of coefficients for the filter is:";COEFF
1000 PRINT
1010 PRINT "        Enter # of coefficients desired ONLY if greater than";COEFF
1020 INPUT "        otherwise, press <Enter> to continue --> ",TEMP
1030 IF TEMP = 0 THEN 1160
1040 IF TEMP < COEFF THEN 1080
1050 COEFF=TEMP
1060 Q=(COEFF-1)/2
1070 GOTO 1160
1080 FOR BLINK=1 TO 10
1090 LOCATE 20
1100 PRINT TAB(12);"ERROR! : Order will not satisfy specifications - reenter"
1110 FOR DELAY=1 TO 100:NEXT DELAY
1120 LOCATE 20
1130 PRINT "                                                                   "
1140 NEXT BLINK
1150 GOTO 980
1160 CLS
1170 REM
1180 REM
1190 LOCATE 12:PRINT TAB(28)"Please wait ...working"
1200 REM
1210 GOSUB 1420       'Routine to calculate FS coefficients, C'(n)
1220 REM
1230 IF WIN = 2 THEN GOSUB 1640    'Hanning
1240 IF WIN = 3 THEN GOSUB 1760    'Hamming
1250 IF WIN = 4 THEN GOSUB 1880    'Blackman
1260 REM
1270 REM                    *** Rearrange the coefficients ***
1280 FOR N=0 TO Q
1290 H(N)=C(Q-N)
1300 NEXT N
1310 REM                   *** Generate the symmetry about Q ***
1320 FOR N=1 TO Q
1330 H(Q+N)=H(Q-N)    '(i.e., H[i] = C[q]-i)
1340 NEXT N
1350 REM                    *** Convert coefficients to Hex ***
1360 GOSUB 3990   'call hex conversion routine
1370 PRINT
1380 GOSUB 2430   'call output menu routine
1390 REM
1400 END
1410 REM ======================= FS Calculation Routine ========================
1420 REM
1430 C(0)=HIGH-LOW
1440 IF TYPE = 4 THEN C(0) = 1-C(0)   'for bandstop
1450 FOR I=1 TO Q
1460 C(I)=(SIN(HIGH*I*PI)/(I*PI))-(SIN(LOW*I*PI)/(I*PI))   'Fourier Series
1470 IF TYPE = 4 THEN C(I)=-C(I)      'for bandstop
1480 NEXT I
1490 RETURN
1500 REM =======================================================================
1510 REM
1520 REM ================== Rectangular Window Routine =========================
1530 REM
1540 REM        This trivial routine is placed here for documentation
1550 REM        purposes only, it is not called from anywhere within the
1560 REM        main program.
1570 REM
1580 REM        The Rectangular window sequence is given by:
1590 REM
1600 REM           W(n) = 1, |n| <= Q; 0, elsewhere
1610 REM
1620 REM =======================================================================
1630 REM
1640 REM ==================== Hanning Window Routine ===========================
1650 REM
1660 REM        The Hanning window sequence is given by:
1670 REM
1680 REM        W(n) = 0.5 + 0.5cos(nPI/Q), |n| <= Q; 0, elsewhere
1690 REM
1700 FOR I=0 TO Q
1710 C(I)=C(I)*(.5+.5*COS(I*PI/Q))      ' C'(n) = C(n)*W(n)
1720 NEXT I
1730 RETURN
1740 REM =======================================================================
1750 REM
1760 REM ==================== Hamming Window Routine ===========================
1770 REM
1780 REM        The Hamming window sequence is given by:
1790 REM
1800 REM        W(n) = 0.54 + 0.46cos(nPI/Q), |n| <= Q; 0, elsewhere
1810 REM
1820 FOR I=0 TO Q
1830 C(I)=C(I)*(.54+.46*COS(I*PI/Q))     ' C'(n) = C(n)*W(n)
1840 NEXT I
1850 RETURN
1860 REM =======================================================================
1870 REM
1880 REM ====================== Blackman Window Routine ========================
1890 REM
1900 REM        The Blackman window sequence is given by:
1910 REM
1920 REM        W(n) = 0.42 + 0.5cos(2nPI/2Q) + 0.08cos(4nPI/2Q)
1930 REM
1940 FOR I=0 TO Q
1950 C(I)=C(I)*(.42+.5*COS((2*I*PI)/(2*Q))+.08*COS((4*I*PI)/(2*Q)))
1960 NEXT I
1970 RETURN
1980 REM =======================================================================
1990 REM
2000 REM =================== Lowpass Prompts Routine ===========================
2010 TYPE$="LOWPASS"
2020 CLS
2030 LOCATE 1
2040 PRINT "        Specifications:"
2050 PRINT TAB(24);TYPE$
2060 LCUT=0
2070 LOCATE 12:INPUT "        Enter the 3-db cutoff frequency in Hz --> ",HCUT
2080 LOCATE 3:PRINT TAB(24);"Cutoff Frequency =";HCUT;"Hz"
2090 RETURN
2100 REM =======================================================================
2110 REM
2120 REM ==================== Highpass Prompts Routine =========================
2130 TYPE$="HIGHPASS"
2140 CLS
2150 LOCATE 1
2160 PRINT "        Specifications:"
2170 PRINT TAB(24);TYPE$
2180 HIGH=1    'for highpass normalized upper cutoff = 1
2190 LOCATE 12:INPUT "        Enter the 3-db cutoff frequency in Hz --> ",LCUT
2200 LOCATE 3:PRINT TAB(24);"Cutoff Frequency =";LCUT;"Hz"
2210 RETURN
2220 REM =======================================================================
2230 REM
2240 REM ================== Bandpass/stop Prompts Routine ======================
2250 REM
2260 TYPE$="BANDPASS":IF TYPE = 4 THEN TYPE$="BANDSTOP"
2270 CLS
2280 LOCATE 1
2290 PRINT "        Specifications:"
2300 PRINT TAB(24);TYPE$
2310 LOCATE 12
2320 INPUT "        Enter the 3-db lower cutoff frequency in Hz --> ",LCUT
2330 LOCATE 12
2340 PRINT "                                                                  "
2350 LOCATE 3:PRINT TAB(24);"Lower Cutoff Frequency =";LCUT;"Hz"
2360 LOCATE 12
2370 INPUT "        Enter the 3-db upper cutoff frequency in Hz --> ",HCUT
2380 IF HCUT <= LCUT THEN GOSUB 4140:GOTO 2360
2390 LOCATE 4:PRINT TAB(24);"Upper Cutoff Frequency =";HCUT;"Hz"
2400 RETURN
2410 REM =======================================================================
2420 REM
2430 REM ========================== Output Routine =============================
2440 REM
2450 REM        This routine allows the user to define the output
2460 REM        path.  These include the terminal, line printer, or
2470 REM        data file, which includes TMS320 data format for simple
2480 REM        merging into a 'blank' filter program.
2490 REM
2500 CLS
2510 LOCATE 8
2520 PRINT "       Send coefficients to:"
2530 PRINT TAB(29);"(S)creen"
2540 PRINT TAB(29);"(P)rinter"
2550 PRINT TAB(29);"(F)ile: contains TMS320 (C25 or C30) data format"
2560 PRINT TAB(29);"(R)eturn to Filter Type Menu"
2570 PRINT TAB(29);"(E)xit to DOS"
2580 PRINT
2590 INPUT "       Enter desired path --> ",PATH$
2600 IF PATH$ = "S" OR PATH$ = "s" THEN GOSUB 2670
2610 IF PATH$ = "P" OR PATH$ = "p" THEN GOSUB 2990
2620 IF PATH$ = "F" OR PATH$ = "f" THEN GOSUB 3400
2630 IF PATH$ = "R" OR PATH$ = "r" THEN GOTO 390
2640 IF PATH$ = "E" OR PATH$ = "e" THEN CLS:SYSTEM
2650 GOTO 2500
2660 REM
2670 REM -------------------- Output Coefficients to Terminal ------------------
2680 REM
2690 CLS
2700 PRINT TAB(13);COEFF;"Coefficient ";TYPE$;" Filter Using the ";WIN$;" Window"
2710 PRINT
2720 PRINT "       Filter Coefficients:"
2730 PRINT
2740 PRINT TAB(20);"C'(n)";TAB(53);"H(n)"
2750 PRINT TAB(20);"-----";TAB(53);"----"
2760 PRINT
2770 PRINT TAB(51);"DECIMAL";TAB(70);"HEX"
2780 PRINT
2790 FOR I=0 TO Q
2800 IF INT(I/11) <> I/11 OR I = 0 THEN 2860
2810 LOCATE 23:PRINT "       Press <Enter> to view";Q-I+1;"remaining coefficients";:INPUT " ",RESPONSE
2820 LOCATE 23
2830 PRINT "                                                                   "
2840 GOSUB 3690    'erase coefficients from screen
2850 REM
2860 GOSUB 3790    'convert coeff. index to char, removing spaces in output
2870 PRINT TAB(8);
2880 PRINT "C'(";C1$;")";
2890 PRINT TAB(15);" = ";
2900 PRINT USING "+#.#####";C(I);:PRINT " = C'(";C2$;")";
2910 PRINT TAB(41);
2920 PRINT "H(";C1$;")";
2930 PRINT TAB(47);" = ";
2940 PRINT USING "+#.#####";H(I);: PRINT " = H(";C2$;")";
2950 PRINT TAB(67);" = ";CHEX$(I)
2960 NEXT I
2970 LOCATE 23:INPUT "       Press <Enter> to return to Output Path Menu ",RET
2980 RETURN
2990 REM --------------------- Output Coefficients to Printer ------------------
3000 REM
3010 REM        This nested routine allows the user to produce hardcopy of the
3020 REM        coefficients.  NOTE: If this option is specified from the
3030 REM        output path menu and the printer is NOT READY, then a
3040 REM        "device timeout" will occur via GWBASIC.  This will cause
3050 REM        exiting to GWBASIC.  To restart the FILTER DEVELOPMENT
3060 REM        PACKAGE type <F2> <Enter>.  This will restart whichever
3070 REM        window module selected previously from Main Menu.  Program
3080 REM        must be re-run to re-develop filter to obtain the hardcopy
3090 REM        initially desired.
3100 REM
3110 LPRINT TAB(26);"*** FIR DEVELOPMENT PACKAGE ***"
3120 LPRINT
3130 LPRINT TAB(14);COEFF;"Coefficient ";TYPE$;" Filter Using the ";WIN$;" Window"
3140 LPRINT
3150 IF TYPE=1 THEN LPRINT TAB(8);"Cutoff frequency =";HCUT;"Hz":GOTO 3190
3160 IF TYPE=2 THEN LPRINT TAB(8);"Cutoff frequency =";LCUT;"Hz":GOTO 3190
3170 LPRINT TAB(8);"Lower Cutoff Frequency =";LCUT;"Hz"
3180 LPRINT TAB(8);"Upper Cutoff Frequency =";HCUT;"Hz"
3190 LPRINT TAB(8);"Sampling Frequency (Fs) =";FSAM;"Hz"
3200 LPRINT TAB(8);"Impulse Duration (D) =";D*1000;"msec"
3210 LPRINT
3220 LPRINT TAB(20);"C'(n)";TAB(53);"H(n)"
3230 LPRINT TAB(20);"-----";TAB(53);"----"
3240 LPRINT
3250 LPRINT TAB(51);"DECIMAL";TAB(70);"HEX"
3260 LPRINT
3270 FOR I=0 TO Q
3280 GOSUB 3790    'convert coeff. index to char, removing spaces in output
3290 LPRINT TAB(8);
3300 LPRINT "C'(";C1$;")";
3310 LPRINT TAB(15);" = ";
3320 LPRINT USING "+#.#####";C(I);:LPRINT " = C'(";C2$;")";
3330 LPRINT TAB(41);
3340 LPRINT "H(";C1$;")";
3350 LPRINT TAB(47);" = ";
3360 LPRINT USING "+#.#####";H(I);:LPRINT " = H(";C2$;")";
3370 LPRINT TAB(67);" = ";CHEX$(I)
3380 NEXT I
3390 RETURN
3400 REM -------------------- Output Coefficients to File ----------------------
3410 REM
3420 REM        This routine writes the coefficients to a file (named by the
3430 REM        user).  The file contains TMS320 (C30 or C25) data format.
3440 REM
3480 REM        This file can then be merged into a 'blank' filter program
3490 REM        via a word processor like PC-WRITE.  NOTE: A directory or
3500 REM        drive specification can be given along with the filename.
3501 LOCATE 17
3502 PRINT "                                                                   "
3503 LOCATE 17
3504 INPUT "      Enter DSP type (C25 or C30):";DSP$
3505 IF DSP$="C30" OR DSP$="c30" THEN GOTO 3660
3506 IF DSP$="C25" OR DSP$="c25" THEN GOTO 3510
3507 GOTO 3501
3510 CLS
3520 LOCATE 12
3530 INPUT "       Enter Filename : ",NAM$
3540 LOCATE 12
3550 PRINT "                                                                   "
3560 LOCATE 12:PRINT TAB(30);"...saving ";NAM$
3570 OPEN NAM$ FOR OUTPUT AS 1
3580 FOR I=0 TO COEFF-1
3590 S$=STR$(I)
3600 IF I>9 AND I<100 THEN PRINT#1,"H";RIGHT$(S$,2);"     DATA    >";CHEX$(I) : GOTO 3630
3610 IF I>99 THEN PRINT#1,"H";RIGHT$(S$,3);"    DATA    >";CHEX$(I):GOTO 3630
3620 PRINT#1,"H";RIGHT$(S$,1);"      DATA    >";CHEX$(I)
3630 NEXT I
3640 CLOSE 1
3650 RETURN
3660 CLS
3661 LOCATE 12
3662 INPUT "       Enter Filename : ",NAM$
3663 LOCATE 12
3664 PRINT "                                                                   "
3665 LOCATE 12:PRINT TAB(30);"...saving ";NAM$
3666 OPEN NAM$ FOR OUTPUT AS 1
3668 PRINT#1, USING "COEFF   .FLOAT  +##.########       ;H(###)  ==>  H(N-1)";H(COEFF-1);COEFF-1
3670 FOR I=COEFF-2 TO 1 STEP -1
3672 PRINT#1, USING "        .FLOAT  +##.########       ;H(###)";H(I);I
3674 NEXT I
3675 PRINT#1, USING "H0      .FLOAT  +##.########       ;H(  0)";H(0)
3677 PRINT#1,"LENGTH  .SET (H0-COEFF)+1          ;LENGTH = # OF COEFFs =";COEFF
3678 CLOSE 1
3679 RETURN
3680 REM =======================================================================
3690 REM ======================== Clear Coefficients ===========================
3700 REM
3710 LOCATE 10
3720 FOR N=1 TO 12
3730 PRINT "                                                                               "
3740 NEXT N
3750 LOCATE 10
3760 RETURN
3770 REM =======================================================================
3780 REM
3790 REM ======================= Integer to Character ==========================
3800 REM
3810 REM        This routine eliminates the spaces produced in the output
3820 REM        when an Integer TYPE is printed.  The spaces before and after
3830 REM        the integer value are set aside in the event that the integer
3840 REM        value is negative.  Since the Coefficient index is always
3850 REM        non-negative, this value is converter to a character type,
3860 REM        and output from this routine as c1$ and c2$.
3870 REM
3880 S1$=STR$(I)
3890 S2$=STR$(COEFF-I-1)      'COEFF-I-1 is the symmetrical alternate of index
3900 C1$=RIGHT$(S1$,1)
3910 C2$=RIGHT$(S2$,1)
3920 IF I > 9 THEN C1$=RIGHT$(S1$,2)
3930 IF I > 99 THEN C1$=RIGHT$(S1$,3)
3940 IF COEFF-I-1 > 9 THEN C2$=RIGHT$(S2$,2)
3950 IF COEFF-I-1 > 99 THEN C2$=RIGHT$(S2$,3)
3960 RETURN
3970 REM =======================================================================
3980 REM
3990 REM ==================== Convert Coefficients to Hex ======================
4000 REM
4010 REM        This routine converts the coefficients to Hexidecimal
4020 REM        via a routine in GWBASIC called HEX$().  The HEX values
4030 REM        are accurate to + or - 1.
4040 REM
4050 FOR M=0 TO COEFF-1
4060 CARY(M)=CINT(H(M)*32768!)    'scale by 2^15
4070 IF CARY(M) >= 0 GOTO 4090
4080 CARY(M)=65536!+CARY(M)
4090 CHEX$(M)=HEX$(CARY(M))
4100 NEXT M
4110 RETURN
4120 REM =======================================================================
4130 REM
4140 REM ====================== Error Message Routine ==========================
4150 REM
4160 FOR BLINK=1 TO 10
4170 IF MES<>1 THEN 4210
4180 LOCATE 20
4190 PRINT TAB(14);"ERROR! : Sampling Frequency (Fs) >= 2*Nyquist - reenter"
4200 GOTO 4230
4210 LOCATE 20
4220 PRINT TAB(15);"ERROR! : Frequency value is inconsistant - reenter"
4230 FOR DELAY=1 TO 100:NEXT DELAY
4240 LOCATE 20
4250 PRINT"                                                                    "
4260 NEXT BLINK
4270 RETURN
4280 REM =======================================================================
