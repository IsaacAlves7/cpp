;MR10MCR.ASM-10-BAND MULTIRATE FILTER USING THE AIC
          .title  "MR10MCR.ASM"       ;10-BAND WITH MACROS
          .option X                   ;FOR XREF
          .global NOISE, FILTER, INTRPLTR, RST_WBFR
          .global BEGIN, LOOP, SET_UP, PROCESS, SETAR7_1, SETAR7_2
          .global PROCLOOP, F_REPEAT, I_REPEAT, EXIT, WAIT
          .global AICSEC,AICSET_I,AICIO_I
          .sect   "VECTORS"           ;ASSEMBLE INTO VECTOR SECTION
RESET     .word   BEGIN               ;RESET VECTOR
          .space  4                   ;SKIP 4 WORDS
TIMER0    .word   ISR                 ;TINT0 VECTOR LOCATION @ 5h
          .space  58                  ;REMAINDER OF VECTOR SECTION
;*****************************************************************************
;*   MACRO DEFINITIONS BELOW                                                 *
;*****************************************************************************
RST_WBFR  $MACRO                      ;MACRO TO PING-PONG WORKING BUFFERS
          LDI     1, R5
          AND     R6, R5
          BNZ     SWITCH?
          LDI     1, R5
          AND     IR1, R5
          BNZ     CHNG_1?
          BD      CONT?
          LDI     @B3_ADR, AR5
          LDI     @B2_ADR, AR6
          NOP
CHNG_1?   BD      CONT?
          LDI     @B2_ADR, AR5
          LDI     @B3_ADR, AR6
          NOP
SWITCH?   LDI     1, R5
          AND     IR1, R5
          BNZ     CHNG_2?
          BD      CONT?
          LDI     @B3_ADR, AR5
          LDI     @B1_ADR, AR6
          NOP
CHNG_2?   BD      CONT?
          LDI     @B1_ADR, AR5
          LDI     @B3_ADR, AR6
          NOP
CONT?     NOP
          $ENDM

NOISE     $MACRO                      ;****RANDOM NOISE GENERATOR****
          LDI     0,R4                ;CLEAR R4
          LDI     R0,R2               ;PUT SEED IN R2
          LSH     -17,R2              ;MOVE BIT 17 TO LSB
          ADDI    R2,R4               ;ADD TO R4 (17)
          LSH     -11,R2              ;MOVE BIT 28 TO LSB
          ADDI    R2,R4               ;ADD TO R4 (28+17)
          LSH     -2,R2               ;MOVE BIT 30 TO LSB
          ADDI    R2,R4               ;ADD TO R4 (30+(28+17))
          LSH     -1,R2               ;MOVE BIT 31 TO LSB
          ADDI    R2,R4               ;ADD TO R4 (31+(30+(28+17)))
          AND     1,R4                ;MASK LSB OF R4
          LDIZ    @MINUS,R2           ;IF R4 = 0     R2 = MINUS VALUE
          LDINZ   @PLUS,R2            ;IF R4 = 1     R2 = PLUS VALUE
          LSH     1,R0                ;SHIFT SEED LEFT 1
          OR      R4,R0               ;PUT R4 INTO LSB OF R0 (SEED)
          $ENDM
;*****************************************************************************
;*   THIS SECTION DEFINES MEMORY FOR SCALE VALUES, COEFFS, ETC.              *
;*****************************************************************************
          .data                       ;ASSEMBLE INTO DATA SECTION
SCALE_1   .float  0.0, 0.0, 0.0       ;FILTER SCALE VALUES  (LOW, MID, UPPER)
SCALE_2   .float  0.0, 0.0, 0.0
SCALE_3   .float  0.0, 0.0, 0.0
SCALE_4   .float  0.0, 0.0, 0.0
SCALE_5   .float  0.0, 0.0, 0.0
SCALE_6   .float  0.0, 0.0, 0.0
SCALE_7   .float  0.0, 0.0, 0.0
SCALE_8   .float  0.0, 0.0, 0.0
SCALE_9   .float  0.0, 0.0, 0.0
SCALE_10  .float  0.0, 1.0, 0.0
                                      ;LOW 1/3 OCTAVE COEFFICIENTS #40 TO #0
LFILT_40  .float  0.000000000,  0.000000000,  0.001861572,  0.009765625
          .float  0.004577636, -0.007446289, -0.005493164, -0.000274658
          .float -0.015899658, -0.019592285,  0.037536621,  0.089355468
          .float  0.005004882, -0.159637451, -0.143707275,  0.124420166
          .float  0.299743652,  0.065185546, -0.322631835, -0.308929443
          .float  0.140655517,  0.421447753,  0.140655517, -0.308929443
          .float -0.322631835,  0.065185546,  0.299743652,  0.124420166
          .float -0.143707275, -0.159637451,  0.005004882,  0.089355468
          .float  0.037536621, -0.019592285, -0.015899658, -0.000274658
          .float -0.005493164, -0.007446289,  0.004577636,  0.009765625
LFILT_0   .float  0.001861572
                                      ;MID 1/3 OCTAVE COEFF 40
MFILT_40  .float -0.001342773, -0.000488281, -0.003692626,  0.003509521
          .float  0.009277343, -0.004486083, -0.007781982,  0.000396728
          .float -0.011413574,  0.009765625,  0.057342529, -0.023529052
          .float -0.131866455,  0.035247802,  0.225738525, -0.038635253
          .float -0.319488525,  0.030212402,  0.389099121, -0.011444091
          .float -0.414855957, -0.011444091,  0.389099121,  0.030212402
          .float -0.319488525, -0.038635253,  0.225738525,  0.035247802
          .float -0.131866455, -0.023529052,  0.057342529,  0.009765625
          .float -0.011413574,  0.000396728, -0.007781982, -0.004486083
          .float  0.009277343,  0.003509521, -0.003692626, -0.000488281
MFILT_0   .float -0.001342773
                                     ;UPPER 1/3 OCTAVE COEFFICIENTS #40 TO #0
UFILT_40  .float  0.001953125, -0.003295898, -0.000396728,  0.001861572
          .float  0.001007080,  0.003356933, -0.010589599,  0.002441406
          .float  0.006927490,  0.003692626,  0.006896972, -0.064880371
          .float  0.055999755,  0.110687255, -0.229766845,  0.020507812
          .float  0.334716796, -0.306884764, -0.166229248,  0.463836669
          .float -0.166229248, -0.306884764,  0.334716796,  0.020507812
          .float -0.229766845,  0.110687255,  0.055999755, -0.064880371
          .float  0.006896972,  0.003692626,  0.006927490,  0.002441406
          .float -0.010589599,  0.003356933,  0.001007080,  0.001861572
          .float -0.000396728, -0.003295898,  0.001953125,  0.000000000
UFILT_0   .float  0.000000000
                                      ;INTRPLTR COEFFICIENTS 1,3,...,21
INTR_01   .float  0.012542724, -0.034423828,  0.079803466, -0.180267334
          .float  0.625427246,  0.625427246, -0.180267334,  0.079803466
          .float -0.034423828,  0.012542724,  0.000000000
                                      ;INTRPLTR COEFFICIENTS 0,2,...,20
INTR_00   .float  0.007263183, -0.012023925,  0.019012451, -0.025939941
          .float  0.030700683,  0.967346191,  0.030700683, -0.025939941
          .float  0.019012451, -0.012023925,  0.007263183
FLT1_IN   .usect  "F_DATA", 64        ;START OF FILTER'S INPUT SAMPLES
FLT2_IN   .usect  "F_DATA", 64
FLT3_IN   .usect  "F_DATA", 64
FLT4_IN   .usect  "F_DATA", 64
FLT5_IN   .usect  "F_DATA", 64
FLT6_IN   .usect  "F_DATA", 64
FLT7_IN   .usect  "F_DATA", 64
FLT8_IN   .usect  "F_DATA", 64
FLT9_IN   .usect  "F_DATA", 64
FLT10_IN  .usect  "F_DATA", 64
INT1_IN   .usect  "I_DATA", 16        ;START OF INTRPLTR'S INPUT SAMPLES
INT2_IN   .usect  "I_DATA", 16
INT3_IN   .usect  "I_DATA", 16
INT4_IN   .usect  "I_DATA", 16
INT5_IN   .usect  "I_DATA", 16
INT6_IN   .usect  "I_DATA", 16
INT7_IN   .usect  "I_DATA", 16
INT8_IN   .usect  "I_DATA", 16
INT9_IN   .usect  "I_DATA", 16
          .bss    F1_COEFF, 41
          .bss    F2_COEFF, 41
          .bss    F3_COEFF, 41
          .bss    F4_COEFF, 41
          .bss    F5_COEFF, 41
          .bss    F6_COEFF, 41
          .bss    F7_COEFF, 41
          .bss    F8_COEFF, 41
          .bss    F9_COEFF, 41
          .bss    F10_COEFF, 41
          .bss    I_COEFF, 22
          .bss    BUFFER_3, 256
          .bss    FCP_STRG, 10
          .bss    FIP_STRG, 10
          .bss    ICP_STRG, 1
          .bss    IIP_STRG, 9
BUFFER_1  .float  0.0
          .space  511
BUFFER_2  .float  0.0
          .space  511
I_LENGTH  .set    11                  ;# OF INTRPLTR COEFFS
F_LENGTH  .set    41                  ;# OF FILTER COEFFS
STACKS    .word   809F00H             ;INIT STACK POINTER DATA
PLUS      .word   000001000H          ;POSITIVE LEVEL (+2^12)
MINUS     .word   0FFFFF000H          ;NEGATIVE LEVEL (-2^12)
SEED      .word   7E521603H           ;INITIAL SEED VALUE
AICSEC    .word   0E1CH,1h,4286H,67H  ;Fs=16234 Hz
SRAM_CNTL .word   808064H             ;PRIMARY BUS CONTROL REGISTER
FLTC_ADR  .word   F1_COEFF            ;POINTS TO START OF FILTER COEFFS
FLTI_ADR  .word   FLT1_IN             ;POINTS TO START OF FILTER INPUTS
INTC_ADR  .word   I_COEFF             ;POINTS TO START OF INTRPLTR COEFFS
INTI_ADR  .word   INT1_IN             ;POINTS TO START OF INTRPLTR INPUTS
B1_ADR    .word   BUFFER_1            ;POINTS TO START OF BUFFER #1
B2_ADR    .word   BUFFER_2            ;POINTS TO START OF BUFFER #2
B3_ADR    .word   BUFFER_3            ;POINTS TO START OF BUFFER #3
FC_PNTR   .word   FCP_STRG            ;POINTS TO START OF FLTR COEFF PNTRS
FI_PNTR   .word   FIP_STRG            ;POINTS TO START OF FLTR INPUT PNTRS
IC_PNTR   .word   ICP_STRG            ;POINTS TO START OF INTRPLTR COEFF PNTRS
II_PNTR   .word   IIP_STRG            ;POINTS TO START OF INTRPLTR INPUT PNTRS
SCALER    .float  0.25                ;OUTPUT SCALE FACTOR
;*****************************************************************************
;*                             INITIALIZATION                                *
;*****************************************************************************
          .text                       ;ASSEMBLE INTO TEXT SECTION
BEGIN     LDP     STACKS              ;INIT DATA PAGE
          LDI     @STACKS, SP         ;SP-> 809F00H
          LDI     -249, R0            ;0 WAIT STATE DATA
          LDI     @SRAM_CNTL, AR0     ;PRIMARY BUS CNTL REG => AR0
          AND     *AR0, R0            ;SET DESIRED BITS
          STI     R0, *AR0            ;REPLACE CONTENTS PRIMARY BUS CNTL REG
;*****************************************************************************
;*                           CREATE ADJUSTED COEFFS                           *
;*    AR0 = SCALE VALUES ADDRESS        AR3 = UPR 1/3 FLTR COEFF ADDR         *
;*    AR1 = LOW 1/3 FLTR COEFF ADDR     AR4 = ADJUSTED COEFF ADDR (10 SETS)   *
;*    AR2 = MID 1/3 FLTR COEFF ADDR     R7 = LOOP COUNTER 1 TO 10(10 FILTERS) *
;*    REGISTERS USED FOR CALCULATIONS: R0, R1, R2, R3, R4, R5, R6             *
;*****************************************************************************
          LDI     SCALE_1, AR0        ;1ST SCALE VALUE ADDR -> AR0
          LDI     @FLTC_ADR, AR4      ;ADJUSTED FLTR COEFFS STARTING ADDR -> AR4
          LDI     10, R7              ;SET LOOP COUNTER TO 10 (10 -> R7)
LOOP      LDI     LFILT_40, AR1       ;LOW 1/3 FLTR COEFFS ADDR -> AR1
          LDI     MFILT_40, AR2       ;MID 1/3 FLTR COEFFS ADDR -> AR2
          LDI     UFILT_40, AR3       ;UPR 1/3 FLTR COEFFS ADDR -> AR3
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

                                      ;MOVE INTRPLTR COEFFS FROM
                                      ;SRAM TO INTERNAL RAM
          LDI     INTR_01, AR0        ;AR0=ADDR OF INTR COEFF IN SRAM
          LDI     @INTC_ADR, AR1      ;AR1=ADDR OF INTR COEFF IN RAM
          LDF     *AR0, R0
          RPTS    I_LENGTH+I_LENGTH-2
            LDF     *++AR0, R0
||          STF     R0, *AR1++
          STF     R0, *AR1
;*****************************************************************************
;*                  CLEAR INPUTS FOR ALL FILTERS AND INTRPLTRS               *
;*    AR0 = FILTER & INTRPLTR INPUT ADDRESS           R0 = 0                 *
;*    NOTE: CLEAR ALL INPUT DATA LOCATIONS FOR EACH FILTER/INTRPLTR          *
;*          (10 FLTRS x 64 LOCATIONS) + (9 INTPLTRS x 16 LOCATIONS)          *
;*****************************************************************************
          LDF     0.0, R0             ;CLEAR R0
          LDI     @FLTI_ADR, AR0      ;ADDR OF FLTRS INPUT DATA -> AR0
          RPTS    10*64+9*16          ;REPEAT CNTR = SIZE OF ALL INPUT BUFFERS
            STF     R0, *AR0++(1)        ;STORE R0 INTO CURRENT ADDRESS
;*****************************************************************************
;*     STORE STARTING DATA ADDRS INTO POINTER ADDR FOR ALL FLTRS/INTPLTRS    *
;*   AR0 = FILTER AND INTRPLTR INPUT DATA AND COEFFS ADDRESSES               *
;*   AR1 = FILTER AND INTRPLTR POINTER STORAGE ADDRESSES                     *
;*****************************************************************************
          LDI     @FC_PNTR, AR1       ;PNTR STORAGE ADDRESS
          LDI     @FLTC_ADR, AR0      ;FILTER COEFF ADDRESS
          LDI     9, RC
          RPTB    RPTB2_E
            STI     AR0, *AR1++
RPTB2_E     ADDI    41, AR0
          LDI     @FLTI_ADR, AR0      ;FILTER INPUT ADDRESS
          LDI     9, RC
          RPTB    RPTB3_E
            STI     AR0, *AR1++
RPTB3_E     ADDI    64, AR0
          LDI     @INTC_ADR, AR0      ;INTRPLTR COEFF ADDRESS
          STI     AR0, *AR1++
          LDI     @INTI_ADR, AR0      ;INTRPLTR INPUT ADDRESS
          LDI     8, RC
          RPTB    RPTB3
            STI     AR0, *AR1++
RPTB3       ADDI    16, AR0
;*****************************************************************************
;*                      MAIN PROCESSING INITIALIZATION                       *
;*****************************************************************************
SET_UP    LDI     @SEED, R0           ;NOISE SEED VALUE -> R0
          LDI     0,R6                ;CLEAR BUFFER
          CALL    AICSET_I            ;INIT AIC
                                      ;NOTE: 1st time thru output buffer => 0
;*****************************************************************************
;*                              MAIN PROCESSING                              *
;*    AR0 = D/A OUTPUT ADDRESS           R0 = SEED FOR RND NOISE GENERATOR   *
;*    AR1 = FLTR & INTRP INPUT DATA PNTR R1 = WORKING REG (USED IN FLTR SUB) *
;*    AR2 = FLTR & INTRP COEFF PNTR      R2 = WORKING REG (USED IN FLTR SUB) *
;*    AR3 = STORAGE PNTR                 R3 = OUTPUT SMPL FROM FILT/INTRPLTR *
;*    AR4 = USED FOR REPEAT LOOPS        R4 = NOT USED                       *
;*    AR5 = WORKING BUFFER POINTER       R5 = WORKING REG                    *
;*    AR6 = WORKING BUFFER POINTER       R6 = BUFFER STATUS REG              *
;*    AR7 = OUTPUT BUFFER POINTER        R7 = OUTPUT SAMPLE                  *
;*                                      IR0 = OUTPUT CNTR (0-512)            *
;*                                      IR1 = INDICATES LAST BAND# COMPLETED *
;*****************************************************************************
PROCESS   LDI     0, IR0              ;****BEGIN PROCESSING****
          LDI     0, IR1              ;SET REGS TO DEFAULT VALUES
          LDI     1, AR4
          LDI     0, RC
	  RST_WBFR
	  LDF     0.0, R5
	  STF     R5, *AR5
          LDI     1, R5
          AND     R6, R5
          BNZ     SETAR7_2
SETAR7_1  LDI     @B1_ADR, AR7
          BR      PROCLOOP
SETAR7_2  LDI     @B2_ADR, AR7

PROCLOOP  RST_WBFR                    ;BEGIN PROCESSING LOOP
          LDI     @FI_PNTR, AR3
          LDI     *+AR3(IR1), AR1
          LDI     F_LENGTH, BK
          LDI     1, R4
          SUBI3   R4, AR4, RC
F_REPEAT  RPTB    RPTB_E1
            LDI     @FC_PNTR, AR3
            LDI     *+AR3(IR1), AR2
            NOISE                       ;RETURNS NOISE SAMPLE IN R2
            FLOAT   R2, R1
            STF     R1, *AR1++%         ;STORE NOISE SAMPLE IN FILTER INPUT
            CALL    FILTER              ;RETURNS FILTERED DATA IN R3
            ADDF    *AR5, R3            ;SUM PREVIOUS BAND'S SAMPLE
RPTB_E1     STF     R3, *AR5++          ;STORE RESULTS
          LDI     @FI_PNTR, AR3
          STI     AR1, *+AR3(IR1) 
          RST_WBFR
          CMPI    9, IR1                ;IF THIS IS BAND 10 THEN BRANCH TO EXIT
          BZ      EXIT                  ;     (INTERPOLATION NOT NEEDED)
          LDI     @II_PNTR, AR3
          LDI     *+AR3(IR1), AR1
          LDI     I_LENGTH, BK
          LDI     1, R4
          SUBI3   R4, AR4, RC
I_REPEAT  RPTB    RPTB_E2
            LDI     @IC_PNTR, AR3
            LDI     *AR3, AR2
            LDF     *AR5++, R3
            STF     R3, *AR1++%
            CALL    INTRPLTR
            STF     R3, *AR6++
            CALL    INTRPLTR
RPTB_E2     STF     R3, *AR6++
          LDI     @II_PNTR, AR3
          BD      PROCLOOP
          STI     AR1, *+AR3(IR1) 
          ADDI    1, IR1
          LSH     1, AR4

EXIT      ADDI    1, IR1
          LDI     512, R5
                                     ;finished filling up next output buffer
                                     ;wait for current output buffer to empty
WAIT      IDLE
          CMPI    R5, IR0
          BNZ     WAIT
          ADDI    1, R6
          B       PROCESS
;*****************************************************************************
;*                        SUBROUTINES DEFINED BELOW                          *
;*****************************************************************************
FILTER    PUSH    ST
	  PUSH    RC
	  PUSH    RS
	  PUSH    RE
          LDF     0.0,R1
          LDF     0.0,R3
          LDI     F_LENGTH-1, RC
          RPTB    RPT_FLT
RPT_FLT     MPYF3   *AR1++%,*AR2++,R1
||          ADDF3   R1,R3,R3
          ADDF3   R1,R3,R3
	  POP     RE
	  POP     RS
	  POP     RC
	  POP     ST
	  OR      2000h, ST
          RETS

INTRPLTR  PUSH    ST
	  PUSH    RC
	  PUSH    RS
	  PUSH    RE
          LDF     0.0, R1
          LDF     0.0, R3
          LDI     I_LENGTH-1, RC
          RPTB    RPTB_I?
RPTB_I?     MPYF3   *AR1++%, *AR2++, R1
       ||   ADDF3   R1, R3, R3
          ADDF3   R1, R3, R3
	  POP     RE
	  POP     RS
	  POP     RC
	  POP     ST
          OR      2000h, ST
          RETS
;*****************************************************************************
;*      ISR                    INTERRUPT SERVICE ROUTINE - OUTPUTS SAMPLE    *
;*****************************************************************************
ISR       PUSH    ST                  ;PUSH ST REG TO STACK
          PUSH    IE                  ;PUSH IE REG TO STACK
          PUSH    R4                  ;PUSH R4 REG TO STACK
          LDF     *AR7++, R4          ;GET SAMPLE FROM OUTPUT BUFFER
          MPYF    @SCALER,R4          ;SCALE OUTPUT SAMPLE
          FIX     R4, R7              ;R7=INTEGER(R4)
          PUSH    R6                  ;SAVE R6
          CALL    AICIO_I             ;CALL AIC I/O ROUTINE
          POP     R6                  ;RESTORE R6
          ADDI    1, IR0              ;INCREMENT OUTPUT COUNTER
          POP     R4                  ;POP LAST R4 REG FROM STACK
          POP     IE                  ;POP LAST IE REG FROM STACK
          POP     ST                  ;POP LAST ST REG FROM STACK
          OR      2000h, ST           ;MAKE SURE GIE BIT REMAINS SET
          RETI                        ;RETURN FROM INTERRUPT
          .end





