          .title  "MR10AIC.ASM"       ;10-BAND MULTIRATE FILTER USING THE AIC
          .option X                   ;FOR XREF
          .global BEGIN, LOOP, SET_UP, PROCESS, WAIT, SWITCH  ;REF/DEF SYMBOLS
          .global RST_WBFR, CHANGE, RETURN, FILTER, INTRPLTR                 
          .global AICSEC,AICSET,AICIO_I
          .sect   "VECTORS"           ;ASSEMBLE INTO VECTOR SECTION
RESET     .word   BEGIN               ;RESET VECTOR
          .space  4                   ;SKIP 4 WORDS
XINT0VEC  .word   ISR                 ;XINT0 VECTOR LOCATION @ 5h
          .space  58                  ;REMAINDER OF VECTOR SECTION
;*****************************************************************************
;*   THIS SECTION DEFINES MEMORY FOR SCALE VALUES, COEFFS, ETC, ETC.         *
;*****************************************************************************
          .data                       ;ASSEMBLE INTO DATA SECTION
SCALE_1L  .float  0.0                 ;SCALE VALUES
SCALE_1M  .float  0.0                 ;  10 FILTERS * 3 BANDS = 30 SCALE VALUES
SCALE_1U  .float  0.0
SCALE_2L  .float  0.0
SCALE_2M  .float  0.0
SCALE_2U  .float  0.0
SCALE_3L  .float  0.0
SCALE_3M  .float  0.0
SCALE_3U  .float  0.0
SCALE_4L  .float  0.0
SCALE_4M  .float  0.0
SCALE_4U  .float  0.0
SCALE_5L  .float  0.0
SCALE_5M  .float  0.0
SCALE_5U  .float  0.0
SCALE_6L  .float  0.0
SCALE_6M  .float  0.0
SCALE_6U  .float  0.0
SCALE_7L  .float  0.0
SCALE_7M  .float  0.0
SCALE_7U  .float  0.0
SCALE_8L  .float  0.0
SCALE_8M  .float  0.0
SCALE_8U  .float  0.0
SCALE_9L  .float  0.0
SCALE_9M  .float  0.0
SCALE_9U  .float  0.0
SCALE_0L  .float  0.0
SCALE_0M  .float  1.0                 ;BAND 10 MID
SCALE_0U  .float  0.0
LFILT_40  .float  0.000000000         ;LOW 1/3 OCTAVE COEFFICIENTS #40 TO #0
          .float  0.000000000
          .float  0.001861572
          .float  0.009765625
          .float  0.004577636
          .float -0.007446289
          .float -0.005493164
          .float -0.000274658
          .float -0.015899658
          .float -0.019592285
          .float  0.037536621
          .float  0.089355468
          .float  0.005004882
          .float -0.159637451
          .float -0.143707275
          .float  0.124420166
          .float  0.299743652
          .float  0.065185546
          .float -0.322631835
          .float -0.308929443
          .float  0.140655517
          .float  0.421447753
          .float  0.140655517
          .float -0.308929443
          .float -0.322631835
          .float  0.065185546
          .float  0.299743652
          .float  0.124420166
          .float -0.143707275
          .float -0.159637451
          .float  0.005004882
          .float  0.089355468
          .float  0.037536621
          .float -0.019592285
          .float -0.015899658
          .float -0.000274658
          .float -0.005493164
          .float -0.007446289
          .float  0.004577636
          .float  0.009765625
LFILT_0   .float  0.001861572
MFILT_40  .float -0.001342773         ;MID 1/3 OCTAVE COEFFICIENTS #40 T0 #0
          .float -0.000488281
          .float -0.003692626
          .float  0.003509521
          .float  0.009277343
          .float -0.004486083
          .float -0.007781982
          .float  0.000396728
          .float -0.011413574
          .float  0.009765625
          .float  0.057342529
          .float -0.023529052
          .float -0.131866455
          .float  0.035247802
          .float  0.225738525
          .float -0.038635253
          .float -0.319488525
          .float  0.030212402
          .float  0.389099121
          .float -0.011444091
          .float -0.414855957
          .float -0.011444091
          .float  0.389099121
          .float  0.030212402
          .float -0.319488525
          .float -0.038635253
          .float  0.225738525
          .float  0.035247802
          .float -0.131866455
          .float -0.023529052
          .float  0.057342529
          .float  0.009765625
          .float -0.011413574
          .float  0.000396728
          .float -0.007781982
          .float -0.004486083
          .float  0.009277343
          .float  0.003509521
          .float -0.003692626
          .float -0.000488281
MFILT_0   .float -0.001342773
UFILT_40  .float  0.001953125        ;UPPER 1/3 OCTAVE COEFFICIENTS #40 TO #0
          .float -0.003295898
          .float -0.000396728
          .float  0.001861572
          .float  0.001007080
          .float  0.003356933
          .float -0.010589599
          .float  0.002441406
          .float  0.006927490
          .float  0.003692626
          .float  0.006896972
          .float -0.064880371
          .float  0.055999755
          .float  0.110687255
          .float -0.229766845
          .float  0.020507812
          .float  0.334716796
          .float -0.306884764
          .float -0.166229248
          .float  0.463836669
          .float -0.166229248
          .float -0.306884764
          .float  0.334716796
          .float  0.020507812
          .float -0.229766845
          .float  0.110687255
          .float  0.055999755
          .float -0.064880371
          .float  0.006896972
          .float  0.003692626
          .float  0.006927490
          .float  0.002441406
          .float -0.010589599
          .float  0.003356933
          .float  0.001007080
          .float  0.001861572
          .float -0.000396728
          .float -0.003295898
          .float  0.001953125
          .float  0.000000000
UFILT_0   .float  0.000000000

INTR_01   .float  0.012542724         ;INTRPLTR COEFFICIENTS #1,3,..,21
INTR_03   .float -0.034423828
INTR_05   .float  0.079803466
INTR_07   .float -0.180267334
INTR_09   .float  0.625427246
INTR_11   .float  0.625427246
INTR_13   .float -0.180267334
INTR_15   .float  0.079803466
INTR_17   .float -0.034423828
INTR_19   .float  0.012542724
INTR_21   .float  0.000000000         ;set to zero
INTR_00   .float  0.007263183         ;INTRPLTR COEFFICIENTS #0,2,...,20
INTR_02   .float -0.012023925
INTR_04   .float  0.019012451
INTR_06   .float -0.025939941
INTR_08   .float  0.030700683
INTR_10   .float  0.967346191
INTR_12   .float  0.030700683
INTR_14   .float -0.025939941
INTR_16   .float  0.019012451
INTR_18   .float -0.012023925
INTR_20   .float  0.007263183

FLT1_IN   .usect  "F_DATA",64         ;START OF FILTER'S INPUT SAMPLES
FLT2_IN   .usect  "F_DATA",64
FLT3_IN   .usect  "F_DATA",64
FLT4_IN   .usect  "F_DATA",64
FLT5_IN   .usect  "F_DATA",64
FLT6_IN   .usect  "F_DATA",64
FLT7_IN   .usect  "F_DATA",64
FLT8_IN   .usect  "F_DATA",64
FLT9_IN   .usect  "F_DATA",64
FLT0_IN   .usect  "F_DATA",64

I1_DATA   .usect  "I_DATA",16         ;START OF INTRPLTR'S INPUT SAMPLES
I2_DATA   .usect  "I_DATA",16
I3_DATA   .usect  "I_DATA",16
I4_DATA   .usect  "I_DATA",16
I5_DATA   .usect  "I_DATA",16
I6_DATA   .usect  "I_DATA",16
I7_DATA   .usect  "I_DATA",16
I8_DATA   .usect  "I_DATA",16
I9_DATA   .usect  "I_DATA",16
          .bss    F1_COEFF,41         ;START OF FILTER'S ADJUSTED COEFFS
          .bss    F2_COEFF,41
          .bss    F3_COEFF,41
          .bss    F4_COEFF,41
          .bss    F5_COEFF,41
          .bss    F6_COEFF,41
          .bss    F7_COEFF,41
          .bss    F8_COEFF,41
          .bss    F9_COEFF,41
          .bss    F0_COEFF,41
          .bss    IC_IRAM,I_LENGTH+I_LENGTH

          .data
PNTR_F1   .word   0                   ;INIT EACH FLTR'S INPUT BUFFER POINTER
PNTR_F2   .word   0
PNTR_F3   .word   0
PNTR_F4   .word   0
PNTR_F5   .word   0
PNTR_F6   .word   0
PNTR_F7   .word   0
PNTR_F8   .word   0
PNTR_F9   .word   0
PNTR_F0   .word   0
PNTR_I1   .word   0                   ;INIT EACH INTRPR'S INPUT BUFFER POINTER
PNTR_I2   .word   0
PNTR_I3   .word   0
PNTR_I4   .word   0
PNTR_I5   .word   0
PNTR_I6   .word   0
PNTR_I7   .word   0
PNTR_I8   .word   0
PNTR_I9   .word   0

BUFFER_1  .float  0.0                 ;START OF BUFFER #1
          .space  511                 ;RESERVE SPACE IN BUFFER
BUFFER_2  .float  0.0                 ;START OF BUFFER #2
          .space  511                 ;RESERVE SPACE IN BUFFER
BUFFER_3  .float  0.0                 ;START OF BUFFER #3
          .space  255                 ;RESERVE SPACE IN BUFFER

I_LENGTH  .set    INTR_21-INTR_01+1   ;# OF INTRPLTR COEFFS
F_LENGTH  .set    UFILT_0-UFILT_40+1  ;# OF FILTER COEFFS
STACKS    .word   809F00H             ;INIT STACK POINTER DATA
SRAM_CNTL .word   808064H             ;PRIMARY BUS CONTROL REGISTER
IE_REG    .word   10H                 ;ENABLE SP TRANSMIT 0 INTERRUPT
ST_REG    .word   2000H               ;SET STATUS REG
PLUS      .word   000001000H          ;POSITIVE "NOISE" LEVEL 
MINUS     .word   0FFFFF000H          ;NEGATIVE "NOISE" LEVEL 
SEED      .word   7E521603H           ;INITIAL SEED VALUE
AICSEC  .word  3H,0E1CH,3H,4286H,3H,67H,0H  ;Fs=16234 kHz USING AIC MODULE
SC_ADDR   .word   SCALE_1L            ;1ST ADDR OF SCALE VALUES
LF_ADDR   .word   LFILT_40            ;1ST ADDR OF LOWER 1/3 FILTER COEFF 40
MF_ADDR   .word   MFILT_40            ;1ST ADDR OF MID 1/3 FILTER COEFF 40
UF_ADDR   .word   UFILT_40            ;1ST ADDR OF UPPER 1/3 FILTER COEFF 40
F1C_ADDR  .word   F1_COEFF            ;FILTER 1 COEFF ADDRESS
F2C_ADDR  .word   F2_COEFF            ;FILTER 2 COEFF ADDRESS
F3C_ADDR  .word   F3_COEFF            ;FILTER 3 COEFF ADDRESS
F4C_ADDR  .word   F4_COEFF            ;FILTER 4 COEFF ADDRESS
F5C_ADDR  .word   F5_COEFF            ;FILTER 5 COEFF ADDRESS
F6C_ADDR  .word   F6_COEFF            ;FILTER 6 COEFF ADDRESS
F7C_ADDR  .word   F7_COEFF            ;FILTER 7 COEFF ADDRESS
F8C_ADDR  .word   F8_COEFF            ;FILTER 8 COEFF ADDRESS
F9C_ADDR  .word   F9_COEFF            ;FILTER 9 COEFF ADDRESS
F0C_ADDR  .word   F0_COEFF            ;FILTER 10 COEFF ADDRESS
IC_ADDR   .word   IC_IRAM             ;FIRST INTRPLTR COEFF ADDR,INTERNAL RAM
IC_SRAM   .word   INTR_01             ;FIRST INTRPLTR COEFF ADDR,SRAM
F1D_ADDR  .word   FLT1_IN             ;1ST FILTER INPUT DATA ADDRESS
F2D_ADDR  .word   FLT2_IN             ;2ND FILTER INPUT DATA ADDRESS
F3D_ADDR  .word   FLT3_IN             ;3RD FILTER INPUT DATA ADDRESS
F4D_ADDR  .word   FLT4_IN             ;4TH FILTER INPUT DATA ADDRESS
F5D_ADDR  .word   FLT5_IN             ;5TH FILTER INPUT DATA ADDRESS
F6D_ADDR  .word   FLT6_IN             ;6TH FILTER INPUT DATA ADDRESS
F7D_ADDR  .word   FLT7_IN             ;7TH FILTER INPUT DATA ADDRESS
F8D_ADDR  .word   FLT8_IN             ;8TH FILTER INPUT DATA ADDRESS
F9D_ADDR  .word   FLT9_IN             ;9TH FILTER INPUT DATA ADDRESS
F0D_ADDR  .word   FLT0_IN             ;10TH FILTER INPUT DATA ADDRESS
I1D_ADDR  .word   I1_DATA             ;1ST INTRPLTR DATA ADDRESS
I2D_ADDR  .word   I2_DATA             ;2ND INTRPLTR DATA ADDRESS
I3D_ADDR  .word   I3_DATA             ;3RD INTRPLTR DATA ADDRESS
I4D_ADDR  .word   I4_DATA             ;4TH INTRPLTR DATA ADDRESS
I5D_ADDR  .word   I5_DATA             ;5TH INTRPLTR DATA ADDRESS
I6D_ADDR  .word   I6_DATA             ;6TH INTRPLTR DATA ADDRESS
I7D_ADDR  .word   I7_DATA             ;7TH INTRPLTR DATA ADDRESS
I8D_ADDR  .word   I8_DATA             ;8TH INTRPLTR DATA ADDRESS
I9D_ADDR  .word   I9_DATA             ;9TH INTRPLTR DATA ADDRESS
F1P_ADDR  .word   PNTR_F1             ;FILTER 1 INPUT BUFFER POINTER ADDRESS
F2P_ADDR  .word   PNTR_F2             ;FILTER 2 INPUT BUFFER POINTER ADDRESS
F3P_ADDR  .word   PNTR_F3             ;FILTER 3 INPUT BUFFER POINTER ADDRESS
F4P_ADDR  .word   PNTR_F4             ;FILTER 4 INPUT BUFFER POINTER ADDRESS
F5P_ADDR  .word   PNTR_F5             ;FILTER 5 INPUT BUFFER POINTER ADDRESS
F6P_ADDR  .word   PNTR_F6             ;FILTER 6 INPUT BUFFER POINTER ADDRESS
F7P_ADDR  .word   PNTR_F7             ;FILTER 7 INPUT BUFFER POINTER ADDRESS
F8P_ADDR  .word   PNTR_F8             ;FILTER 8 INPUT BUFFER POINTER ADDRESS
F9P_ADDR  .word   PNTR_F9             ;FILTER 9 INPUT BUFFER POINTER ADDRESS
F0P_ADDR  .word   PNTR_F0             ;FILTER 10 INPUT BUFFER POINTER ADDRESS
I1P_ADDR  .word   PNTR_I1             ;INTRPLTR 1 INPUT BUFFER POINTER ADDRESS
I2P_ADDR  .word   PNTR_I2             ;INTRPLTR 2 INPUT BUFFER POINTER ADDRESS
I3P_ADDR  .word   PNTR_I3             ;INTRPLTR 3 INPUT BUFFER POINTER ADDRESS
I4P_ADDR  .word   PNTR_I4             ;INTRPLTR 4 INPUT BUFFER POINTER ADDRESS
I5P_ADDR  .word   PNTR_I5             ;INTRPLTR 5 INPUT BUFFER POINTER ADDRESS
I6P_ADDR  .word   PNTR_I6             ;INTRPLTR 6 INPUT BUFFER POINTER ADDRESS
I7P_ADDR  .word   PNTR_I7             ;INTRPLTR 7 INPUT BUFFER POINTER ADDRESS
I8P_ADDR  .word   PNTR_I8             ;INTRPLTR 8 INPUT BUFFER POINTER ADDRESS
I9P_ADDR  .word   PNTR_I9             ;INTRPLTR 9 INPUT BUFFER POINTER ADDRESS
B1_ADDR   .word   BUFFER_1            ;BUFFER #1 ADDRESS
B2_ADDR   .word   BUFFER_2            ;BUFFER #2 ADDRESS
B3_ADDR   .word   BUFFER_3            ;BUFFER #3 ADDRESS
;*****************************************************************************
;*                             INITIALIZATION                                *
;*****************************************************************************
          .text                       ;ASSEMBLE INTO TEXT SECTION
BEGIN     LDP     STACKS              ;INIT DATA PAGE
          LDI     @STACKS, SP         ;SP-> 809F00H
;FOR 0 WAIT STATE **********************
          LDI     -249,R0             ;0 WAIT STATE DATA
          LDI     @SRAM_CNTL,AR0      ;PRIMARY BUS CNTL REG => AR0
          AND     *AR0,R0             ;"AND" CONTROL REGISTER WITH 0 WAIT DATA
          STI     R0,*AR0             ;INITIATE SRAM TO 0 WAIT STATE
;*****************************************************************************
;*                           CREATE ADJUSTED COEFFS                          *
;*    AR0 = SCALE VALUES ADDRESS        AR3 = UPR 1/3 FLTR COEFF ADDR        *
;*    AR1 = LOW 1/3 FLTR COEFF ADDR     AR4 = ADJUSTED COEFF ADDR (10 SETS)  *
;*    AR2 = MID 1/3 FLTR COEFF ADDR     R7 = LOOP COUNTER 1 TO 10(10 FILTERS)*
;*    REGISTERS USED FOR CALCULATIONS: R0, R1, R2, R3, R4, R5, R6            *
;*****************************************************************************
          LDI     @SC_ADDR, AR0       ;SCALE VALUE ADDR -> AR0
          LDI     @F1C_ADDR, AR4      ;ADJUSTED FLTR COEFFS ADDR -> AR4
          LDI     10, R7              ;SET LOOP COUNTER TO 10 (# OF FLTRS)
LOOP      LDI     @LF_ADDR, AR1       ;LOW 1/3 FLTR COEFFS ADDR -> AR1
          LDI     @MF_ADDR, AR2       ;MID 1/3 FLTR COEFFS ADDR -> AR2
          LDI     @UF_ADDR, AR3       ;UPR 1/3 FLTR COEFFS ADDR -> AR3
          LDF     *AR0++(1), R0       ;LOW 1/3 SCALE VALUE -> R0
          LDF     *AR0++(1), R1       ;MID 1/3 SCALE VALUE -> R1
          LDF     *AR0++(1), R2       ;UPR 1/3 SCALE VALUE -> R2
          LDI     F_LENGTH-1, RC      ;LENGTH-1 -> RC
          RPTB    RPTB1_E             ;REPEAT BLOCK ("LENGTH" TIMES)
            LDF     0.0, R6             ;CLEAR R6
            LDF     *AR1++(1), R3       ;LOW 1/3 FLTR COEFF -> R3
            MPYF    R0, R3              ;COEFF * SCALE VALUE -> R3
            LDF     *AR2++(1), R4       ;MID 1/3 FLTR COEFF -> R4
            MPYF    R1, R4              ;COEFF * SCALE VALUE -> R4
            LDF     *AR3++(1), R5       ;UPR 1/3 FLTR COEFF -> R5
            MPYF    R2, R5              ;COEFF * SCALE VALUE -> R5
            ADDF3   R3, R4, R6          ;ADD R3 + R4 -> R6
            ADDF    R5, R6              ;ADD R5 + R6 -> R6
RPTB1_E     STF     R6, *AR4++(1)       ;STORE RESULTS
          SUBI    1, R7               ;DECREMENT LOOP COUNTER (R7-1)
          BNZ     LOOP                ;BRANCH ON NOT ZERO TO LOOP
;*********** MOVE INTRPLTR COEFF FROM SRAM TO INTERNAL RAM *****************
            LDI  @IC_SRAM,AR0         ;AR0=ADDR INTR COEFF IN SRAM
            LDI  @IC_ADDR,AR1         ;AR1=ADDR INTR COEFF IN RAM
            LDF  *AR0,R0              
            RPTS  I_LENGTH+I_LENGTH-2 
            LDF  *++AR0,R0
||          STF  R0,*AR1++
            STF  R0,*AR1
;*****************************************************************************
;*               CLEAR INPUT BUFFER FOR ALL FILTERS AND INTRPLTRS            *
;*    AR0 = INPUT BUFFER ADDRESS            R0 = 0                           *
;*    NOTE: CLEAR ALL INPUT DATA LOCATIONS FOR EACH FILTER/INTRPLTR          *
;*           10 FLTRS x 64 LOCATIONS + 9 INTPLTRS x 16 LOCATIONS             *
;*****************************************************************************
          LDF     0.0, R0             ;CLEAR R0
          LDI     @F1D_ADDR, AR0      ;FLTR 1 DATA ADDR -> AR0
          RPTS    10*64+9*16          ;REPEAT CNTR = SIZE OF ALL INPUT BUFFERS
            STF     R0, *AR0++(1)        ;STORE R0 INTO CURRENT ADDRESS
;clear all registers
          LDF     0.0, R0             ;CLEAR REG
          LDF     0.0, R1             ;CLEAR REG
          LDF     0.0, R2             ;CLEAR REG
          LDF     0.0, R3             ;CLEAR REG
          LDF     0.0, R4             ;CLEAR REG
          LDF     0.0, R5             ;CLEAR REG
          LDF     0.0, R6             ;CLEAR REG
          LDF     0.0, R7             ;CLEAR REG
;*****************************************************************************
;*  STORE EACH FLTR/INTPLTR STARTING INPUT DATA ADDRESS INTO POINTER ADDRESS *
;*   AR0 = FLTR/INTPLR PNTR ADDR         AR1 = FLTR/INTPLR INPUT DATA ADDR   *
;*   NOTE: STORE STARTING DATA ADDR INTO POINTER ADDR FOR ALL FLTRS/INTPLTRS *
;*****************************************************************************
          LDI     @F1P_ADDR, AR0      ;FILTER 1 POINTER ADDRESS -> AR0
          LDI     @F1D_ADDR, AR1      ;FILTER 1 DATA ADDRESS -> AR1
          STI     AR1, *AR0++         ;STORE AR1 -> AR0++
          LDI     @F2D_ADDR, AR1      ;FILTER 2 DATA ADDRESS -> AR1
          STI     AR1, *AR0++         ;STORE AR1 -> AR0++
          LDI     @F3D_ADDR, AR1      ;FILTER 3 DATA ADDRESS -> AR1
          STI     AR1, *AR0++         ;STORE AR1 -> AR0++
          LDI     @F4D_ADDR, AR1      ;FILTER 4 DATA ADDRESS -> AR1
          STI     AR1, *AR0++         ;STORE AR1 -> AR0++
          LDI     @F5D_ADDR, AR1      ;FILTER 5 DATA ADDRESS -> AR1
          STI     AR1, *AR0++         ;STORE AR1 -> AR0++
          LDI     @F6D_ADDR, AR1      ;FILTER 6 DATA ADDRESS -> AR1
          STI     AR1, *AR0++         ;STORE AR1 -> AR0++
          LDI     @F7D_ADDR, AR1      ;FILTER 7 DATA ADDRESS -> AR1
          STI     AR1, *AR0++         ;STORE AR1 -> AR0++
          LDI     @F8D_ADDR, AR1      ;FILTER 8 DATA ADDRESS -> AR1
          STI     AR1, *AR0++         ;STORE AR1 -> AR0++
          LDI     @F9D_ADDR, AR1      ;FILTER 9 DATA ADDRESS -> AR1
          STI     AR1, *AR0++         ;STORE AR1 -> AR0++
          LDI     @F0D_ADDR, AR1      ;FILTER 10 DATA ADDRESS -> AR1
          STI     AR1, *AR0++         ;STORE AR1 -> AR0++
          LDI     @I1D_ADDR, AR1      ;INTRPLTR 1 DATA ADDRESS -> AR1
          STI     AR1, *AR0++         ;STORE AR1 -> AR0++
          LDI     @I2D_ADDR, AR1      ;INTRPLTR 2 DATA ADDRESS -> AR1
          STI     AR1, *AR0++         ;STORE AR1 -> AR0++
          LDI     @I3D_ADDR, AR1      ;INTRPLTR 3 DATA ADDRESS -> AR1
          STI     AR1, *AR0++         ;STORE AR1 -> AR0++
          LDI     @I4D_ADDR, AR1      ;INTRPLTR 4 DATA ADDRESS -> AR1
          STI     AR1, *AR0++         ;STORE AR1 -> AR0++
          LDI     @I5D_ADDR, AR1      ;INTRPLTR 5 DATA ADDRESS -> AR1
          STI     AR1, *AR0++         ;STORE AR1 -> AR0++
          LDI     @I6D_ADDR, AR1      ;INTRPLTR 6 DATA ADDRESS -> AR1
          STI     AR1, *AR0++         ;STORE AR1 -> AR0++
          LDI     @I7D_ADDR, AR1      ;INTRPLTR 7 DATA ADDRESS -> AR1
          STI     AR1, *AR0++         ;STORE AR1 -> AR0++
          LDI     @I8D_ADDR, AR1      ;INTRPLTR 8 DATA ADDRESS -> AR1
          STI     AR1, *AR0++         ;STORE AR1 -> AR0++
          LDI     @I9D_ADDR, AR1      ;INTRPLTR 9 DATA ADDRESS -> AR1
          STI     AR1, *AR0++         ;STORE AR1 -> AR0++
;*****************************************************************************
;*                      MAIN PROCESSING INITIALIZATION                        *
;*****************************************************************************
SET_UP    LDI     @SEED, R0           ;NOISE SEED VALUE -> R0
          LDI     @B1_ADDR, AR5       ;BUFFER #1 STARTING ADDR -> AR5
          LDI     0, R6               ;CLEAR BUFFER STATUS REG
          LDI     0, IR0              ;CLEAR OUTPUT COUNTER         
          CALL    AICSET              ;SET UP AIC
          LDI     @IE_REG, IE         ;ENABLE XINT0 INTERRUPT
          OR      @ST_REG, ST         ;SET STATUS REG (ENABLE GIE BIT)
                                      ;NOTE: 1st time thru output buffer => 0
;*****************************************************************************
;*                              MAIN PROCESSING                              *
;*    AR0 = OUTPUT ADDRESS               R0 = SEED FOR RND NOISE GENERATOR   *
;*    AR1 = FLTR & INTRPLTR COEFF PNTR   R1 = WORKING REG (USED IN FLTR SUB) *
;*    AR2 = INPUT SMPL BUFFER PNTR       R2 = WORKING REG (USED IN FLTR SUB) *
;*    AR3 = STORAGE ADDR FOR SMPL PNTR   R3 = OUTPUT SMPL FROM FILT/INTRPLTR *
;*    AR4 = NOT USED                     R4 = OUTPUT SAMPLE (FP) USED IN ISR *
;*    AR5 = BUFFER [1-2] POINTER (512)   R5 = WORKING REG                    *
;*    AR6 = BUFFER [1-2] POINTER (512)   R6 = BUFFER STATUS REG              *
;*    AR7 = BUFFER 3 POINTER     (256)   R7 = OUTPUT SAMPLE (INT) USED IN ISR*
;*                                      IR0 = OUTPUT CNTR (0-512)            *
;*                                      IR1 = INDICATES CURRENT BAND #       *
;*****************************************************************************
;band 1
PROCESS   LDI     1, IR1              ;BAND # -> IR1           ==DEBUG==
          CALL    RST_WBFR            ;RESET BUFFERS TO TOP
          LDI     F_LENGTH, BK        ;FILTER LENGTH -> BK
          LDI     @F1C_ADDR, AR1      ;FILTER 1'S COEFF ADDRESS -> AR1
          LDI     @PNTR_F1, AR2       ;FILTER 1'S LAST DATA ADDR -> AR2
          CALL    FILTER              ;RUN FIR FILTER SUBROUTINE -> R3
          LDI     @F1P_ADDR, AR3      ;FILTER 1'S POINTER ADDRESS -> AR3
          STI     AR2, *AR3           ;STORE FILTER 1'S CURRENT DATA ADDRESS
          LDI     I_LENGTH, BK        ;INTRPLTR LENGTH -> BK
          LDI     @IC_ADDR, AR1       ;INTRPLTR COEFF ADDRESS  -> AR1
          LDI     @PNTR_I1, AR2       ;INTRPLTR 1'S LAST DATA ADDR -> AR2
          STF     R3, *AR2++%         ;STORE SAMPLE IN INTRPLTR 1'S INPUT BFFR
          CALL    INTRPLTR            ;RUN 2:1 INTERPOLATION -> R3
          STF     R3, *AR6++          ;STORE OUTPUT TO BUFFER
          CALL    INTRPLTR            ;RUN 2:1 INTERPOLATION -> R3
          STF     R3, *AR6++          ;STORE OUTPUT TO BUFFER
          LDI     @I1P_ADDR, AR3      ;INTRPLTR 1'S PNTR ADDRESS -> AR3
          STI     AR2, *AR3           ;STORE CURRENT INTRPLTR DATA ADDRESS
;band 2
          LDI     2, IR1              ;BAND # -> IR1            ==DEBUG==
          CALL    RST_WBFR            ;RESET BUFFERS TO TOP
          LDI     F_LENGTH, BK        ;FILTER LENGTH -> BK
          LDI     @PNTR_F2, AR2       ;FILTER 2'S LAST DATA ADDR -> AR2
          LDI     1, RC               ;1 -> RC
          RPTB    RPTB2_E             ;REPEAT BLOCK
            LDI     @F2C_ADDR, AR1      ;FILTER 2'S COEFF ADDRESS -> AR1
            CALL    FILTER              ;RUN FIR FILTER SUBROUTINE -> R3
            ADDF    *AR6, R3            ;ADD LAST BAND'S SMPL TO LATEST SMPL
RPTB2_E     STF     R3, *AR6++          ;STORE SUMMED SAMPLE BACK IN BUFFER
          LDI     @F2P_ADDR, AR3      ;FILTER 2'S POINTER ADDRESS -> AR3
          STI     AR2, *AR3           ;STORE FILTER 2'S CURRENT DATA ADDRESS
          CALL    RST_WBFR            ;RESET BUFFER TO TOP
          LDI     I_LENGTH, BK        ;INTRPLTR LENGTH -> BK
          LDI     @PNTR_I2, AR2       ;INTRPLTR 2'S LAST DATA ADDR -> AR2
          LDI     1, RC               ;1 -> RC
          RPTB    RPTB3_E             ;REPEAT BLOCK
            LDI     @IC_ADDR, AR1       ;INTRPLTR COEFF ADDRESS  -> AR1
            LDF     *AR6++, R3          ;SUMMED SAMPLE -> R3
            STF     R3, *AR2++%         ;STORE SMPL IN INTRPLTR 2'S INPUT BFFR
            CALL    INTRPLTR            ;RUN 2:1 INTERPOLATION -> R3
            STF     R3, *AR7++          ;STORE OUTPUT TO BUFFER
            CALL    INTRPLTR            ;RUN 2:1 INTERPOLATION -> R3
RPTB3_E     STF     R3, *AR7++          ;STORE OUTPUT TO BUFFER
          LDI     @I2P_ADDR, AR3      ;INTRPLTR 2'S PNTR ADDRESS -> AR3
          STI     AR2, *AR3           ;STORE CURRENT INTRPLTR DATA ADDRESS
;band 3
          LDI     3, IR1              ;BAND # -> IR1            ==DEBUG==
          CALL    RST_WBFR            ;RESET BUFFERS TO TOP
          LDI     F_LENGTH, BK        ;FILTER LENGTH -> BK
          LDI     @PNTR_F3, AR2       ;FILTER 3'S LAST DATA ADDR -> AR2
          LDI     3, RC               ;3 -> RC
          RPTB    RPTB4_E             ;REPEAT BLOCK
            LDI     @F3C_ADDR, AR1      ;FILTER 3'S COEFF ADDRESS -> AR1
            CALL    FILTER              ;RUN FIR FILTER SUBROUTINE -> R3
            ADDF    *AR7, R3            ;ADD LAST BAND'S SMPL TO LATEST SMPL
RPTB4_E     STF     R3, *AR7++          ;STORE SUMMED SAMPLE BACK IN BUFFER
          LDI     @F3P_ADDR, AR3      ;FILTER 3'S POINTER ADDRESS -> AR3
          STI     AR2, *AR3           ;STORE FILTER 3'S CURRENT DATA ADDRESS
          CALL    RST_WBFR            ;RESET BUFFER TO TOP
          LDI     I_LENGTH, BK        ;INTRPLTR LENGTH -> BK
          LDI     @PNTR_I3, AR2       ;INTRPLTR 3'S LAST DATA ADDR -> AR2
          LDI     3, RC               ;3 -> RC
          RPTB    RPTB5_E             ;REPEAT BLOCK
            LDI     @IC_ADDR, AR1       ;INTRPLTR COEFF ADDRESS  -> AR1
            LDF     *AR7++, R3          ;SUMMED SAMPLE -> R3
            STF     R3, *AR2++%         ;STORE SMPL IN INTRPLTR 3'S INPUT BFFR
            CALL    INTRPLTR            ;RUN 2:1 INTERPOLATION -> R3
            STF     R3, *AR6++          ;STORE OUTPUT TO BUFFER
            CALL    INTRPLTR            ;RUN 2:1 INTERPOLATION -> R3
RPTB5_E     STF     R3, *AR6++          ;STORE OUTPUT TO BUFFER
          LDI     @I3P_ADDR, AR3      ;INTRPLTR 3'S PNTR ADDRESS -> AR3
          STI     AR2, *AR3           ;STORE CURRENT INTRPLTR DATA ADDRESS
;band 4
          LDI     4, IR1              ;BAND # -> IR1            ==DEBUG==
          CALL    RST_WBFR            ;RESET BUFFERS TO TOP
          LDI     F_LENGTH, BK        ;FILTER LENGTH -> BK
          LDI     @PNTR_F4, AR2       ;FILTER 4'S LAST DATA ADDR -> AR2
          LDI     7, RC               ;7 -> RC
          RPTB    RPTB6_E             ;REPEAT BLOCK
            LDI     @F4C_ADDR, AR1      ;FILTER 4'S COEFF ADDRESS -> AR1
            CALL    FILTER              ;RUN FIR FILTER SUBROUTINE -> R3
            ADDF    *AR6, R3            ;ADD LAST BAND'S SMPL TO LATEST SMPL
RPTB6_E     STF     R3, *AR6++          ;STORE SUMMED SAMPLE BACK IN BUFFER
          LDI     @F4P_ADDR, AR3      ;FILTER 4'S POINTER ADDRESS -> AR3
          STI     AR2, *AR3           ;STORE FILTER 4'S CURRENT DATA ADDRESS
          CALL    RST_WBFR            ;RESET BUFFER TO TOP
          LDI     I_LENGTH, BK        ;INTRPLTR LENGTH -> BK
          LDI     @PNTR_I4, AR2       ;INTRPLTR 4'S LAST DATA ADDR -> AR2
          LDI     7, RC               ;7 -> RC
          RPTB    RPTB7_E             ;REPEAT BLOCK
            LDI     @IC_ADDR, AR1       ;INTRPLTR COEFF ADDRESS  -> AR1
            LDF     *AR6++, R3          ;SUMMED SAMPLE -> R3
            STF     R3, *AR2++%         ;STORE SMPL IN INTRPLTR 4'S INPUT BFFR
            CALL    INTRPLTR            ;RUN 2:1 INTERPOLATION -> R3
            STF     R3, *AR7++          ;STORE OUTPUT TO BUFFER
            CALL    INTRPLTR            ;RUN 2:1 INTERPOLATION -> R3
RPTB7_E     STF     R3, *AR7++          ;STORE OUTPUT TO BUFFER
          LDI     @I4P_ADDR, AR3      ;INTRPLTR 4'S PNTR ADDRESS -> AR3
          STI     AR2, *AR3           ;STORE CURRENT INTRPLTR DATA ADDRESS
;band 5
          LDI     5, IR1              ;BAND # -> IR1            ==DEBUG==
          CALL    RST_WBFR            ;RESET BUFFERS TO TOP
          LDI     F_LENGTH, BK        ;FILTER LENGTH -> BK
          LDI     @PNTR_F5, AR2       ;FILTER 5'S LAST DATA ADDR -> AR2
          LDI     15, RC              ;15 -> RC
          RPTB    RPTB8_E             ;REPEAT BLOCK
            LDI     @F5C_ADDR, AR1      ;FILTER 5'S COEFF ADDRESS -> AR1
            CALL    FILTER              ;RUN FIR FILTER SUBROUTINE -> R3
            ADDF    *AR7, R3            ;ADD LAST BAND'S SMPL TO LATEST SMPL
RPTB8_E     STF     R3, *AR7++          ;STORE SUMMED SAMPLE BACK IN BUFFER
          LDI     @F5P_ADDR, AR3      ;FILTER 5'S POINTER ADDRESS -> AR3
          STI     AR2, *AR3           ;STORE FILTER 5'S CURRENT DATA ADDRESS
          CALL    RST_WBFR            ;RESET BUFFER TO TOP
          LDI     I_LENGTH, BK        ;INTRPLTR LENGTH -> BK
          LDI     @PNTR_I5, AR2       ;INTRPLTR 5'S LAST DATA ADDR -> AR2
          LDI     15, RC              ;15 -> RC
          RPTB    RPTB9_E             ;REPEAT BLOCK
            LDI     @IC_ADDR, AR1       ;INTRPLTR COEFF ADDRESS  -> AR1
            LDF     *AR7++, R3          ;SUMMED SAMPLE -> R3
            STF     R3, *AR2++%         ;STORE SMPL IN INTRPLTR 5'S INPUT BFFR
            CALL    INTRPLTR            ;RUN 2:1 INTERPOLATION -> R3
            STF     R3, *AR6++          ;STORE OUTPUT TO BUFFER
            CALL    INTRPLTR            ;RUN 2:1 INTERPOLATION -> R3
RPTB9_E     STF     R3, *AR6++          ;STORE OUTPUT TO BUFFER
          LDI     @I5P_ADDR, AR3      ;INTRPLTR 5'S PNTR ADDRESS -> AR3
          STI     AR2, *AR3           ;STORE CURRENT INTRPLTR DATA ADDRESS
;band 6
          LDI     6, IR1              ;BAND # -> IR1            ==DEBUG==
          CALL    RST_WBFR            ;RESET BUFFERS TO TOP
          LDI     F_LENGTH, BK        ;FILTER LENGTH -> BK
          LDI     @PNTR_F6, AR2       ;FILTER 6'S LAST DATA ADDR -> AR2
          LDI     31, RC              ;31 -> RC
          RPTB    RPTB10_E            ;REPEAT BLOCK
            LDI     @F6C_ADDR, AR1      ;FILTER 6'S COEFF ADDRESS -> AR1
            CALL    FILTER              ;RUN FIR FILTER SUBROUTINE -> R3
            ADDF    *AR6, R3            ;ADD LAST BAND'S SMPL TO LATEST SMPL
RPTB10_E    STF     R3, *AR6++          ;STORE SUMMED SAMPLE BACK IN BUFFER
          LDI     @F6P_ADDR, AR3      ;FILTER 6'S POINTER ADDRESS -> AR3
          STI     AR2, *AR3           ;STORE FILTER 6'S CURRENT DATA ADDRESS
          CALL    RST_WBFR            ;RESs4 BUFFER TO TOP
          LDI     I_LENGTH, BK        ;INTRPLTR LENGTH -> BK
          LDI     @PNTR_I6, AR2       ;INTRPLTR 6'S LAST DATA ADDR -> AR2
          LDI     31, RC              ;31 -> RC
          RPTB    RPTB11_E            ;REPEAT BLOCK
            LDI     @IC_ADDR, AR1       ;INTRPLTR COEFF ADDRESS  -> AR1
            LDF     *AR6++, R3          ;SUMMED SAMPLE -> R3
            STF     R3, *AR2++%         ;STORE SMPL IN INTRPLTR 6'S INPUT BFFR
            CALL    INTRPLTR            ;RUN 2:1 INTERPOLATION -> R3
            STF     R3, *AR7++          ;STORE OUTPUT TO BUFFER
            CALL    INTRPLTR            ;RUN 2:1 INTERPOLATION -> R3
RPTB11_E    STF     R3, *AR7++          ;STORE OUTPUT TO BUFFER
          LDI     @I6P_ADDR, AR3      ;INTRPLTR 6'S PNTR ADDRESS -> AR3
          STI     AR2, *AR3           ;STORE CURRENT INTRPLTR DATA ADDRESS
;band 7
          LDI     7, IR1              ;BAND # -> IR1            ==DEBUG==
          CALL    RST_WBFR            ;RESET BUFFERS TO TOP
          LDI     F_LENGTH, BK        ;FILTER LENGTH -> BK
          LDI     @PNTR_F7, AR2       ;FILTER 7'S LAST DATA ADDR -> AR2
          LDI     63, RC              ;63 -> RC
          RPTB    RPTB12_E            ;REPEAT BLOCK
            LDI     @F7C_ADDR, AR1      ;FILTER 7'S COEFF ADDRESS -> AR1
            CALL    FILTER              ;RUN FIR FILTER SUBROUTINE -> R3
            ADDF    *AR7, R3            ;ADD LAST BAND'S SMPL TO LATEST SMPL
RPTB12_E    STF     R3, *AR7++          ;STORE SUMMED SAMPLE BACK IN BUFFER
          LDI     @F7P_ADDR, AR3      ;FILTER 7'S POINTER ADDRESS -> AR3
          STI     AR2, *AR3           ;STORE FILTER 7'S CURRENT DATA ADDRESS
          CALL    RST_WBFR            ;RESET BUFFER TO TOP
          LDI     I_LENGTH, BK        ;INTRPLTR LENGTH -> BK
          LDI     @PNTR_I7, AR2       ;INTRPLTR 7'S LAST DATA ADDR -> AR2
          LDI     63, RC              ;63 -> RC
          RPTB    RPTB13_E            ;REPEAT BLOCK
            LDI     @IC_ADDR, AR1       ;INTRPLTR COEFF ADDRESS  -> AR1
            LDF     *AR7++, R3          ;SUMMED SAMPLE -> R3
            STF     R3, *AR2++%         ;STORE SMPL IN INTRPLTR 7'S INPUT BFFR
            CALL    INTRPLTR            ;RUN 2:1 INTERPOLATION -> R3
            STF     R3, *AR6++          ;STORE OUTPUT TO BUFFER
            CALL    INTRPLTR            ;RUN 2:1 INTERPOLATION -> R3
RPTB13_E    STF     R3, *AR6++          ;STORE OUTPUT TO BUFFER
          LDI     @I7P_ADDR, AR3      ;INTRPLTR 7'S PNTR ADDRESS -> AR3
          STI     AR2, *AR3           ;STORE CURRENT INTRPLTR DATA ADDRESS
;band 8
          LDI     8, IR1              ;BAND # -> IR1            ==DEBUG==
          CALL    RST_WBFR            ;RESET BUFFERS TO TOP
          LDI     F_LENGTH, BK        ;FILTER LENGTH -> BK
          LDI     @PNTR_F8, AR2       ;FILTER 8'S LAST DATA ADDR -> AR2
          LDI     127, RC             ;127 -> RC
          RPTB    RPTB14_E            ;REPEAT BLOCK
            LDI     @F8C_ADDR, AR1      ;FILTER 8'S COEFF ADDRESS -> AR1
            CALL    FILTER              ;RUN FIR FILTER SUBROUTINE -> R3
            ADDF    *AR6, R3            ;ADD LAST BAND'S SMPL TO LATEST SMPL
RPTB14_E    STF     R3, *AR6++          ;STORE SUMMED SAMPLE BACK IN BUFFER
          LDI     @F8P_ADDR, AR3      ;FILTER 8'S POINTER ADDRESS -> AR3
          STI     AR2, *AR3           ;STORE FILTER 8'S CURRENT DATA ADDRESS
          CALL    RST_WBFR            ;RESET BUFFER TO TOP
          LDI     I_LENGTH, BK        ;INTRPLTR LENGTH -> BK
          LDI     @PNTR_I8, AR2       ;INTRPLTR 8'S LAST DATA ADDR -> AR2
          LDI     127, RC             ;127 -> RC
          RPTB    RPTB15_E            ;REPEAT BLOCK
            LDI     @IC_ADDR, AR1       ;INTRPLTR COEFF ADDRESS  -> AR1
            LDF     *AR6++, R3          ;SUMMED SAMPLE -> R3
            STF     R3, *AR2++%         ;STORE SMPL IN INTRPLTR 8'S INPUT BFFR
            CALL    INTRPLTR            ;RUN 2:1 INTERPOLATION -> R3
            STF     R3, *AR7++          ;STORE OUTPUT TO BUFFER
            CALL    INTRPLTR            ;RUN 2:1 INTERPOLATION -> R3
RPTB15_E    STF     R3, *AR7++          ;STORE OUTPUT TO BUFFER
          LDI     @I8P_ADDR, AR3      ;INTRPLTR 8'S PNTR ADDRESS -> AR3
          STI     AR2, *AR3           ;STORE CURRENT INTRPLTR DATA ADDRESS
;band 9
          LDI     9, IR1              ;BAND # -> IR1            ==DEBUG==
          CALL    RST_WBFR            ;RESET BUFFERS TO TOP
          LDI     F_LENGTH, BK        ;FILTER LENGTH -> BK
          LDI     @PNTR_F9, AR2       ;FILTER 9'S LAST DATA ADDR -> AR2
          LDI     255, RC             ;255 -> RC
          RPTB    RPTB16_E            ;REPEAT BLOCK
            LDI     @F9C_ADDR, AR1      ;FILTER 9'S COEFF ADDRESS -> AR1
            CALL    FILTER              ;RUN FIR FILTER SUBROUTINE -> R3
            ADDF    *AR7, R3            ;ADD LAST BAND'S SMPL TO LATEST SMPL
RPTB16_E    STF     R3, *AR7++          ;STORE SUMMED SAMPLE BACK IN BUFFER
          LDI     @F9P_ADDR, AR3      ;FILTER 9'S POINTER ADDRESS -> AR3
          STI     AR2, *AR3           ;STORE FILTER 9'S CURRENT DATA ADDRESS
          CALL    RST_WBFR            ;RESET BUFFER TO TOP
          LDI     I_LENGTH, BK        ;INTRPLTR LENGTH -> BK
          LDI     @PNTR_I9, AR2       ;INTRPLTR 9'S LAST DATA ADDR -> AR2
          LDI     255, RC             ;255 -> RC
          RPTB    RPTB17_E            ;REPEAT BLOCK
            LDI     @IC_ADDR, AR1       ;INTRPLTR COEFF ADDRESS  -> AR1
            LDF     *AR7++, R3          ;SUMMED SAMPLE -> R3
            STF     R3, *AR2++%         ;STORE SMPL IN INTRPLTR 8'S INPUT BFFR
            CALL    INTRPLTR            ;RUN 2:1 INTERPOLATION -> R3
            STF     R3, *AR6++          ;STORE OUTPUT TO BUFFER
            CALL    INTRPLTR            ;RUN 2:1 INTERPOLATION -> R3
RPTB17_E    STF     R3, *AR6++          ;STORE OUTPUT TO BUFFER
          LDI     @I9P_ADDR, AR3      ;INTRPLTR 9'S PNTR ADDRESS -> AR3
          STI     AR2, *AR3           ;STORE CURRENT INTRPLTR DATA ADDRESS
;band 10
          LDI     10, IR1             ;BAND # -> IR1            ==DEBUG==
          CALL    RST_WBFR            ;RESET BUFFERS TO TOP
          LDI     F_LENGTH, BK        ;FILTER LENGTH -> BK
          LDI     @PNTR_F0, AR2       ;FILTER 10'S LAST DATA ADDR -> AR2
          LDI     511, RC             ;511 -> RC
          RPTB    RPTB18_E            ;REPEAT BLOCK
            LDI     @F0C_ADDR, AR1      ;FILTER 10'S COEFF ADDRESS -> AR1
            CALL    FILTER              ;RUN FIR FILTER SUBROUTINE -> R3
            ADDF    *AR6, R3            ;ADD LAST BAND'S SMPL TO LATEST SMPL
RPTB18_E    STF     R3, *AR6++          ;STORE SUMMED SAMPLE BACK IN BUFFER
          LDI     @F0P_ADDR, AR3      ;FILTER 10'S POINTER ADDRESS -> AR3
          STI     AR2, *AR3           ;STORE FILTER 10'S CURRENT DATA ADDRESS
;finished filling up next output buffer
;wait for current output buffer to empty
          LDF     0.0, R1             ;CLEAR REG                ==DEBUG==
          LDF     0.0, R2             ;CLEAR REG                ==DEBUG==
          LDF     0.0, R3             ;CLEAR REG                ==DEBUG==
          LDF     0.0, R4             ;CLEAR REG                ==DEBUG==
          LDI     512, R5             ;LOAD REG                 ==DEBUG==
          LDI     15, IR1             ;INDICATES FINISHED       ==DEBUG==
WAIT      IDLE                        ;WAIT FOR INTERRUPT
          SUBI3   R5, IR0, R7         ;CHECK IF OUTPUT CNTR = 256
          BNZ     WAIT                ;IF NOT GOTO WAIT
;current output buffer empty,
;switch output sample pointer to new output buffer
          LDI     0, IR0              ;RESET OUTPUT COUNTER
          LDI     1, R5               ;1 --> R5
          AND     R6, R5              ;MASK BIT 1 OF R6 (NONDESTRUCTIVE)
          BZ      SWITCH              ;IF ZERO BRANCH TO SWITCH
          ADDI    1, R6               ;INCREMENT BUFFER STATUS REG
          BRD     PROCESS             ;DELAYED BRANCH BACK TO PROCESS
          LDI     @B1_ADDR, AR5       ;BUFFER #1 STARTING ADDR -> AR5
          LDI     @B2_ADDR, AR6       ;BUFFER #2 STARTING ADDR -> AR6
          LDI     @B3_ADDR, AR7       ;BUFFER #3 STARTING ADDR -> AR7
SWITCH    ADDI    1, R6               ;INCREMENT BUFFER STATUS REG
          BRD     PROCESS             ;DELAYED BRANCH BACK TO PROCESS
          LDI     @B2_ADDR, AR5       ;BUFFER #2 STARTING ADDR -> AR5
          LDI     @B1_ADDR, AR6       ;BUFFER #1 STARTING ADDR -> AR6
          LDI     @B3_ADDR, AR7       ;BUFFER #3 STARTING ADDR -> AR7
;*****************************************************************************
;*                                SUBROUTINES                                *
;*****************************************************************************

;*****************************************************************************
;*      RST_WBFR               RESETS BUFFERS TO TOP ADDRESSES               *
;*****************************************************************************
;reset working buffers
RST_WBFR  PUSH    ST                  ;PUSH ST REG TO TOP OF STACK
          LDI     1, R5               ;1 --> R5
          AND     R6, R5              ;MASK BIT 1 OF R6 (NONDESTRUCTIVE)
          BNZ     CHANGE              ;IF NOT ZERO BRANCH TO CHANGE
          BRD     RETURN              ;DELAYED BRANCH TO RETURN
          LDI     @B2_ADDR, AR6       ;BUFFER #2 STARTING ADDR -> AR6
          LDI     @B3_ADDR, AR7       ;BUFFER #3 STARTING ADDR -> AR7
          POP     ST                  ;POP LAST ST REG FROM STACK
CHANGE    LDI     @B1_ADDR, AR6       ;BUFFER #1 STARTING ADDR -> AR6
          LDI     @B3_ADDR, AR7       ;BUFFER #3 STARTING ADDR -> AR7
          POP     ST                  ;POP LAST ST REG FROM STACK
RETURN    OR      2000h, ST           ;MAKE SURE GIE BIT REMAINS SET
          RETS
;*****************************************************************************
;*      FILTER                 GENERATES RANDOM NOISE --> FIR FILTER --> R3  *
;*****************************************************************************
;generate noise sample, store into filter's input buffer, run fir filter
FILTER    PUSH    ST                  ;PUSH CURRENT ST REG TO TOP OF STACK
          PUSH    RC                  ;PUSH CURRENT RC REG TO TOP OF STACK
          PUSH    RS                  ;PUSH CURRENT RS REG TO TOP OF STACK
          PUSH    RE                  ;PUSH CURRENT RE REG TO TOP OF STACK
;generate noise sample
          LDI     0, R1               ;INIT R1=0
          LDI     R0, R2              ;PUT SEED IN R2
          LSH     -17, R2             ;MOVE BIT 17 TO LSB     ->R2
          ADDI    R2, R1              ;ADD BIT (17)           ->R1
          LSH     -11, R2             ;MOVE BIT 28 TO LSB     ->R2
          ADDI    R2, R1              ;ADD BITS (28+17)       ->R1
          LSH     -2, R2              ;MOVE BIT 30 TO LSB     ->R2
          ADDI    R2, R1              ;ADD BITS (30+28+17)    ->R1
          LSH     -1, R2              ;MOVE BIT 31 TO LSB     ->R2
          ADDI    R2, R1              ;ADD BITS (31+30+28+17) ->R1
          AND     1, R1               ;MASK LSB OF R1
          LDIZ    @MINUS, R2          ;IF R1 = 0, R2 = @MINUS
          LDINZ   @PLUS, R2           ;IF R1 = 1, R2 = @PLUS
          LSH     1, R0               ;SHIFT SEED LEFT BY 1
          OR      R1, R0              ;R1 (FEEDBACK BIT) --> LSB OF R0 (SEED)
;store noise sample into filter's input buffer
          FLOAT   R2, R1              ;FLOAT R2 -> R1
          STF     R1, *AR2++%         ;STORE R1 INTO FILTER'S INPUT BUFFER
;run fir filter
          LDF     0.0, R1             ;INIT R1=0
          LDF     0.0, R3             ;INIT R3=0
          LDI     F_LENGTH-1, RC      ;LOAD FILTER LENGTH INTO REPEAT COUNTER
          RPTB    RPT_END1            ;NEXT 2 INSTRUCTIONS F_LENGTH TIMES
RPT_END1    MPYF3   *AR1++, *AR2++%, R1 ;H(N-1-i)*x(n-(N-1-i))),i=0..N-1
||          ADDF3   R1, R3, R3          ;IN // WITH ACC -> R3
          ADDF3   R1, R3, R3          ;ADD LAST PRODUCT = Y(n) -> R3
          POP     RE                  ;POP LAST RE REG FROM STACK
          POP     RS                  ;POP LAST RS REG FROM STACK
          POP     RC                  ;POP LAST RC REG FROM STACK
          POP     ST                  ;POP LAST ST REG FROM STACK
          OR      2000h, ST           ;MAKE SURE GIE BIT REMAINS SET
          RETS
;*****************************************************************************
;*      INTRPLTR               INTERPOLATOR FIR FILTER --> R3                *
;*****************************************************************************
INTRPLTR  PUSH    ST                  ;PUSH CURRENT ST REG TO TOP OF STACK
          PUSH    RC                  ;PUSH CURRENT RC REG TO TOP OF STACK
          PUSH    RS                  ;PUSH CURRENT RS REG TO TOP OF STACK
          PUSH    RE                  ;PUSH CURRENT RE REG TO TOP OF STACK
          LDF     0.0, R1             ;INIT R1=0
          LDF     0.0, R3             ;INIT R3=0
          LDI     I_LENGTH-1, RC      ;LOAD INTRPLTR LENGTH INTO REPEAT COUNTER
          RPTB    RPT_END2            ;NEXT 2 INSTRUCTIONS I_LENGTH TIMES
RPT_END2    MPYF3   *AR1++, *AR2++%, R1 ;H(N-1-i)*x(n-(N-1-i))),i=0..N-1
||          ADDF3   R1, R3, R3          ;IN // WITH ACC -> R3
          ADDF3   R1, R3, R3          ;ADD LAST PRODUCT = Y(n) -> R3
          POP     RE                  ;POP LAST RE REG FROM STACK
          POP     RS                  ;POP LAST RS REG FROM STACK
          POP     RC                  ;POP LAST RC REG FROM STACK
          POP     ST                  ;POP LAST ST REG FROM STACK
          OR      2000h, ST           ;MAKE SURE GIE BIT REMAINS SET
          RETS
;*****************************************************************************
;*      ISR                    INTERRUPT SERVICE ROUTINE - OUTPUTS SAMPLE    *
;*****************************************************************************
ISR       PUSH    ST                  ;PUSH ST REG TO STACK
          PUSH    IE                  ;PUSH IE REG TO STACK
          LDF     *AR5++, R4          ;GET SAMPLE FROM OUTPUT BUFFER
          FIX     R4, R7              ;R7=INTEGER(R4)
          PUSH    R6                  ;SAVE R6 (AICIO_I CHANGES R6)
          CALL    AICIO_I             ;OUTPUT R7 AND INPUT R6 (AIC)
          POP     R6                  ;RESTORE R6
          ADDI    1, IR0              ;INCREMENT OUTPUT COUNTER
          POP     IE                  ;POP LAST IE REG FROM STACK
          POP     ST                  ;POP LAST ST REG FROM STACK
          OR      2000h, ST           ;MAKE SURE GIE BIT REMAINS SET
          RETI                        ;RETURN FROM INTERRUPT
          .end                        ;END


