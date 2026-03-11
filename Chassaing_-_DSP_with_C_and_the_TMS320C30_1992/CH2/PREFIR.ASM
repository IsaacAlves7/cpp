;PREFIR.ASM - BACKGROUND FOR FILTER,NOISE(4 SAMPLES)INPUT
          .TITLE  "PREFIR.ASM"       ;BACKGROUND FOR CONVO 
          .GLOBAL BEGIN,MAIN         ;REF/DEF SYMBOLS   
          .DATA                      ;ASSEMBLE INTO DATA SECTION
HN        .FLOAT  1                  ;1ST VALUE IN HN ARRAY (H3)
          .FLOAT  2                  ;2ND HN VALUE (H2)
          .FLOAT  0                  ;3RD HN VALUE (H1)
H0        .FLOAT  3                  ;4TH HN VALUE (H0)
H_ADDR    .WORD   HN                 ;STARTING ADDR OF HN ARRAY
X0_ADDR   .WORD   XN+LENGTH-1        ;ADDR OF BOTTOM SAMPLE X(n) 
IN_ADDR   .WORD   IN                 ;INPUT PORT ADDR
OUT_ADDR  .WORD   OUT                ;OUTPUT PORT ADDR
LENGTH    .SET    H0-HN+1            ;# OF HN VALUES
XN        .USECT  "XN_BUFFER",LENGTH ;BUFFER SIZE OF SAMPLES XN
          .SECT   "VECTORS"          ;ASSEMBLE INTO VECT SECTION
MAIN      .WORD   BEGIN              ;START OF CODE
IN        .USECT  "IN_PORT",1        ;1 SPACE FOR INPUT PORT ADDR
OUT       .USECT  "OUT_PORT",1       ;1 SPACE FOR OUTPUT PORT ADDR
;           +--------+       +------------+
; LOW ADDR  | H3 = 1 |       |   X(n-3)   | 
;           +--------+       +------------+
;           | H2 = 2 |       |   X(n-2)   | 
;           |--------+       +------------+ 
;           | H1 = 0 :       |   X(n-1)   |              
;           +--------+       +------------+
; HIGH ADDR | H0 = 3 |       |   X(n)     |
;           +--------+       +------------+
          .TEXT                      ;ASSEMBLE INTO TEXT SECTION
BEGIN     LDP     IN_ADDR            ;INIT DATA PAGE (0)
          LDI     @X0_ADDR,AR1       ;BOTTOM SAMPLE X(n) ADDR->AR1     
          LDI     LENGTH,BK          ;BK = 4 (LENGTH OF CIRC BUFFER)
          LDF     0,R0               ;INIT R0=0
          RPTS    LENGTH-1           ;NEXT INST. LENGTH TIMES
          STF     R0,*AR1--%         ;INIT SAMPLES=0,AR1->X(n) AGAIN
          LDI     @IN_ADDR,AR5       ;INPUT PORT ADDR -> AR5=804000h
          LDI     @OUT_ADDR,AR6      ;OUTPUT PORT ADDR-> AR6=804002h
          LDI     LENGTH,R4          ;R4=4,LOOP COUNTER 
LOOP      FLOAT   *AR5,R3            ;INPUT FROM PORT ADDR->R3=4096
          STF     R3,*AR1++%         ;NEWEST SAMPLE @ 809803
          LDI     @H_ADDR,AR0        ;START ADDR OF HN->AR0
          CALL    FILTER             ;GO TO SUBROUTINE FILTER
          FIX     R2,R1              ;R1=INTEGER(R2)
          STI     R1,*AR6            ;OUTPUT FROM @ 804002H
          SUBI    1,R4               ;DECREMENT R4
          BNZ     LOOP               ;BRANCH BACK UNTIL R4=0
WAIT      BR      WAIT               ;WAIT INDEFINITELY
;SUBROUTINE FILTER
FILTER    LDF     0,R0               ;INIT R0=0
          LDF     0,R2               ;INIT R2=0
          RPTS    LENGTH-1           ;EXECUTE NEXT 2 INSTR 4 TIMES
          MPYF    *AR0++,*AR1++%,R0  ;(AR0)*(AR1) -> R0
||        ADDF    R0,R2              ;IN PARALLEL WITH ACC IN R2
          ADDF    R0,R2              ;LAST ACC->R2=H3*X(n-3)+..+H0*X(n)
          RETS                       ;RETURN FROM SUBROUTINE
          .END                       ;END





