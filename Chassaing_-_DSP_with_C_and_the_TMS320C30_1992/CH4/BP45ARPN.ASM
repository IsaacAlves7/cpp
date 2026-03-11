;BP45ARPN.ASM-FIR BANDPASS FILTER WITH NOISE GENERATOR
        .TITLE  "BP45ARPN.ASM"   ;BP FIR, Fc=Fs/10, Fs=8 kHz
        .GLOBAL MAIN,BEGIN,AICSET_I,AICSEC,SPSET,AICIO_I,FILT
        .DATA                    ;ASSEMBLE INTO DATA SECTION
XN_ADDR .WORD   XN+LENGTH-1      ;(LAST) SAMPLE ADDRESS
HN_ADDR .WORD   COEFF            ;ADDR OF COEFF H(N-1)  
AICSEC  .WORD   1428h,1h,5EBEh,67h   ;AIC CONFIG DATA 
PLUS    .WORD   1000H            ;POSITIVE NOISE LEVEL
MINUS   .WORD   0FFFFF000H       ;NEGATIVE NOISE LEVEL
SEED    .WORD   7E521603H        ;INITIAL SEED VALUE
SCALER  .FLOAT  0.25             ;OUTPUT SCALE FACTOR
XN      .USECT  "XN_BUFF",LENGTH ;BUFFER SIZE OF SAMPLES  
        .SECT   "VECTORS"        ;VECTOR SECTION
MAIN    .WORD   BEGIN            ;BEGIN @ RESET (0H)    
        .SPACE  4                ;SKIP 4 SPACES
        .WORD   ISR              ;INTERRUPT @ 5H
        .SPACE  58               ;REMAINDER OF VECTOR/TRAP 
        .TEXT                    ;ASSEMBLE INTO TEXT
BEGIN   LDP     SPSET            ;INIT DATA PAGE
        CALL    AICSET_I         ;INIT AIC
        LDI     LENGTH,BK        ;SIZE OF CIRCULAR BUFFER
        LDI     @XN_ADDR,AR1     ;LAST SAMPLE ADDR ->AR1
        LDI     @SEED,R0         ;R0=INITIAL SEED VALUE
WAIT    IDLE                     ;WAIT FOR INTERRUPT
        BR      WAIT             ;BRANCH BACK TO WAIT
;INTERRUPT VECTOR, START OF NOISE GENERATOR 
ISR     LDI     0,R4           ;INIT R4=0
        LDI     R0,R2          ;PUT SEED IN R2
        LSH     -17,R2         ;MOVE BIT 17 TO LSB    =>R2
        ADDI    R2,R4          ;ADD BIT (17)          =>R4
        LSH     -11,R2         ;MOVE BIT 28 TO LSB    =>R2
        ADDI    R2,R4          ;ADD BITS (28+17)      =>R4
        LSH     -2,R2          ;MOVE BIT 30 TO LSB    =>R2
        ADDI    R2,R4          ;ADD BITS (30+28+17)   =>R4
        LSH     -1,R2          ;MOVE BIT 31 TO LSB    =>R2
        ADDI    R2,R4          ;ADD BITS (31+30+28+17)=>R4
        AND     1,R4           ;MASK LSB OF R4
        LDIZ    @MINUS,R7      ;IF R4 = 0, R7 = @MINUS
        LDINZ   @PLUS,R7       ;IF R4 = 1, R7 = @PLUS
        LSH     1,R0           ;SHIFT SEED LEFT BY 1
        OR      R4,R0          ;PUT R4 INTO LSB OF R0 
;MAIN SECTION
        FLOAT  R7,R3           ;INPUT NEW NOISE SAMPLE
        STF     R3,*AR1++%     ;STORE NEWEST SAMPLE
        LDI     @HN_ADDR,AR0   ;AR0 POINTS TO COEFF H(N-1)
        CALL    FILT           ;CALL FILTER ROUTINE  
        MPYF    @SCALER,R3     ;SCALE OUTPUT SAMPLE
        FIX     R3,R7          ;R7=INTEGER(R3)
        CALL    AICIO_I        ;AIC I/O ROUTINE FOR OUTPUT
        RETI                   ;RETURN FROM INTERRUPT
;FILTER ROUTINE
FILT    LDF     0,R1           ;INIT R1=0
        LDF     0,R3           ;INIT R3=0
        RPTS    LENGTH-1       ;NEXT 2 INSTR. LENGTH TIMES
        MPYF    *AR0++,*AR1++%,R1   ;R1=HN*XN
||      ADDF    R1,R3,R3       ;R3 IS ACCUMULATOR
        ADDF    R1,R3,R3       ;LAST VALUE ACCUMULATED
        RETS                   ;RETURN FROM SUBROUTINE
        .DATA                  ;ASSEMBLE INTO DATA SECTION
COEFF   .FLOAT  -1.839E-3      ;H44
 .FLOAT  -2.657E-3,-1.437E-7, 3.154E-3, 2.595E-3,-4.159E-3,-1.540E-2
 .FLOAT  -2.507E-2,-2.547E-2,-1.179E-2, 1.392E-2, 4.206E-2, 5.888E-2
 .FLOAT   5.307E-2, 2.225E-2,-2.410E-2,-6.754E-2,-8.831E-2,-7.475E-2
 .FLOAT  -2.956E-2, 3.030E-2, 8.050E-2, 1.000E-1, 8.050E-2, 3.030E-2
 .FLOAT  -2.956E-2,-7.475E-2,-8.831E-2,-6.754E-2,-2.410E-2, 2.225E-2
 .FLOAT   5.307E-2, 5.888E-2, 4.206E-2, 1.392E-2,-1.179E-2,-2.547E-2
 .FLOAT  -2.507E-2,-1.540E-2,-4.159E-3, 2.595E-3, 3.154E-3,-1.437E-7
        .FLOAT  -2.657E-3        ;H1
H0      .FLOAT  -1.839E-3        ;H0
LENGTH  .SET    H0-COEFF+1       ;LENGTH = 45
        .END                     ;END




