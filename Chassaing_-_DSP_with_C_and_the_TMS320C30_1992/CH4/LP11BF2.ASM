;LP11BF2.ASM-LOWPASS FIR WITH 11 COEFF AND 2 CIRCULAR BUFFERS
          .TITLE "LP11BF2.ASM"     ;LP @ 1 kHz, 2 BUFFERS
          .GLOBAL MAIN,BEGIN,FILT,COEFF ;REF/DEF SYMBOLS
          .DATA                    ;ASSEMBLE -> DATA SECT
IO_IN     .WORD   804000H          ;INPUT ADDRESS
IO_OUT    .WORD   804001H          ;OUTPUT ADDRESS
XN_ADDR   .WORD   XN+LENGTH-1      ;(LAST) SAMPLE ADDR
HN1_ADDR  .WORD   COEFF            ;COEFF TABLE ADDR IN SRAM
HN2_ADDR  .WORD   HN               ;COEFF TABLE ADDR IN RAM    
XN        .USECT  "XN_BUFF",LENGTH ;BUFFER SIZE OF SAMPLES 
HN        .USECT  "HN_BUFF",LENGTH ;COEFF BUFFER SIZE IN RAM
          .SECT   "VECTORS"        ;ASSEMBLE -> VECT SECT
MAIN      .WORD   BEGIN            ;BEGIN @ RESET (0H)
          .TEXT                    ;ASSEMBLE -> TEXT SECT
BEGIN     LDI     @IO_IN,AR5       ;INPUT ADDR -> AR5
          LDI     @IO_OUT,AR6      ;OUTPUT ADDR  -> AR6
          LDI     @XN_ADDR,AR1     ;"LAST"SAMPLE ADDR->AR1
          LDI     LENGTH,BK        ;BK=SIZE OF CIRC BUFFER
          LDF     0,R0             ;INIT R0=0
          RPTS    LENGTH-1         ;INIT ALL SAMPLES TO ZERO
          STF     R0,*AR1--%       ;START/END AT LAST ADDR
          LDI     LENGTH,R4        ;COUNTER FOR FILT SUB     
          LDI     LENGTH-1,RC      ;REPEAT COUNTER=LENGTH
          LDI     @HN1_ADDR,AR0    ;1ST COEFF ADDR IN SRAM
          LDI     @HN2_ADDR,AR2    ;1ST COEFF ADDR IN RAM
          RPTB    EXCH             ;TRANSFER BLOCK SRAM->RAM
          LDF     *AR0++,R5        ;COEFF IN SRAM->R5
EXCH      STF     R5,*AR2++        ;STORE IN RAM
          LDI     @HN2_ADDR,AR0    ;1ST COEFF ADDR->AR0    
LOOP      FLOAT   *AR5,R3          ;INPUT NEW SAMPLE
          STF     R3,*AR1++%       ;STORE TO MODEL DELAY
          CALL    FILT             ;GO TO SUBROUTINE FILT
          FIX     R2,R1            ;R1=INTEGER(R2)
          STI     R1,*AR6          ;OUTPUT INTEGER VALUE
          SUBI    1,R4             ;DECREMENT R4
          BNZ     LOOP             ;BRANCH UNTIL R4 < 0
WAIT      BR      WAIT             ;WAIT
      ;SUBROUTINE FILT          
FILT      LDF     0,R0             ;INIT R0=0
          LDF     0,R2             ;INIT R2=0
          RPTS    LENGTH-1         ;N MULTIPLY
          MPYF    *AR0++%,*AR1++%,R0 ; HN*XN -> R0
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
          .END                     ;END




