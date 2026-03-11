;LP11.ASM-LOWPASS FIR WITH 11 COEFFICIENTS
          .TITLE "LP11.ASM"             ;LP @ 1 kHz
          .GLOBAL MAIN,BEGIN,FILT,COEFF ;REF/DEF SYMBOLS
          .DATA                    ;ASSEMBLE -> DATA SECT
IMP_ADDR  .WORD   IMPULSE          ;INIT VALUE IMPULSE
IO_OUT    .WORD   804001H          ;OUTPUT ADDRESS
XN_ADDR   .WORD   XN+LENGTH-1      ;(LAST) SAMPLE ADDR
HN_ADDR   .WORD   COEFF            ;COEFF TABLE ADDR     
XN        .USECT  "XN_BUFF",LENGTH ;BUFFER SIZE OF SAMPLES
          .SECT   "VECTORS"        ;ASSEMBLE -> VECT SECT
MAIN      .WORD   BEGIN            ;BEGIN @ RESET (0H)
          .TEXT                    ;ASSEMBLE -> TEXT SECT
BEGIN     LDI     @IMP_ADDR,AR5    ;IMPULSE ADDR -> AR5
          LDI     @IO_OUT,AR6      ;OUTPUT ADDR  -> AR6
          LDI     @XN_ADDR,AR1     ;"LAST"SAMPLE ADDR->AR1
          LDI     LENGTH,BK        ;BK=SIZE OF CIRC BUFFER
          LDF     0,R0             ;INIT R0=0
          LDF     0,R7             ;INIT R7=0
          LDI     LENGTH,R4        ;COUNTER FOR FILT SUB
          LDF     *AR5,R3          ;INPUT IMPULSE VALUE
          STF     R3,*AR1++(1)%    ;STORE IMPULSE
LOOP      LDI     @HN_ADDR,AR0     ;COEFF H(N-1) ADDR->AR0
          CALL    FILT             ;GO TO SUBROUTINE FILT
          FIX     R2,R1            ;R1=INTEGER(R2)
          STI     R1,*AR6          ;OUTPUT INTEGER VALUE
          SUBI    1,R4             ;DECREMENT R4
          STF     R7,*AR1++(1)%    ;ALL OTHER SAMPLES=0
          BNZ     LOOP             ;BRANCH UNTIL R4 < 0
WAIT      BR      WAIT             ;WAIT
      ;SUBROUTINE FILT          
FILT      LDF     0,R0             ;INIT R0=0
          LDF     0,R2             ;INIT R2=0
          RPTS    LENGTH-1         ;N MULTIPLY
          MPYF    *AR0++,*AR1++%,R0 ; HN*XN -> R0
||        ADDF    R0,R2            ;// WITH ACC -> R2
          ADDF    R0,R2            ;LAST ACC -> R2
          RETS                     ;RETURN FROM SUBROUTINE
          .DATA                    ;ASSEMBLE -> DATA SECT
COEFF     .FLOAT   0.0             ;H10
          .FLOAT   0.0468          ;H9
          .FLOAT   0.1009          ;H8
          .FLOAT   0.1514          ;H7
          .FLOAT   0.1872          ;H6
          .FLOAT   0.2             ;H5
          .FLOAT   0.1872          ;H4
          .FLOAT   0.1514          ;H3
          .FLOAT   0.1009          ;H2
          .FLOAT   0.0468          ;H1
H0        .FLOAT   0.0             ;H0
LENGTH    .SET    H0-COEFF+1       ;LENGTH=11
IMPULSE   .FLOAT  10000            ;IMPULSE VALUE
          .END                     ;END




