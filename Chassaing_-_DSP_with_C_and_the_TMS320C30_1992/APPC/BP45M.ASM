;BP45M.ASM-MAIN PROGRAM TO CALL FIR SUBROUTINE FROM HYPERSIGNAL
                .TITLE "BP45M.ASM"      ;MAIN FILE FOR FIR
                .GLOBAL MAIN,BEGIN,FILT,COEFF ;REF/DEF SYMBOLS
                .DATA                   ;ASSEMBLE INTO DATA SECT
IN_ADDR         .WORD   IN              ;INIT VALUE IN
OUT_ADDR        .WORD   OUT             ;INIT VALUE OUT
XN_ADDR         .WORD   XN+LENGTH-1     ;NEWEST(LAST) IN SAMPLE 
HN_ADDR         .WORD   COEFF           ;COEFF TABLE ADDR
LENGTH          .SET    45              ;FILTER LENGTH 
XN              .USECT  "XN_BUFF",LENGTH ;SAMPLES BUFFER SIZE
                .SECT   "VECS"          ;ASSEMBLE INTO VECT SECT
MAIN            .WORD   BEGIN
IN              .USECT  "IN_PORT",1     ;1 SPACE FOR IN
OUT             .USECT  "OUT_PORT",1    ;1 SPACE FOR OUT  
                .TEXT                   ;ASSEMBLE INTO TEXT SECT
BEGIN           LDP     IN_ADDR         ;INIT DATA PAGE
                LDI     @IN_ADDR,AR5    ;INPUT ADDR -> AR5
                LDI     @OUT_ADDR,AR6   ;OUTPUT ADDR -> AR6
                LDI     @XN_ADDR,AR1    ;"LAST"SAMPLE ADDR
                LDI     LENGTH,BK       ;BK = FILTER LENGTH
                LDF     0.0,R0          ;INIT R0=0
                RPTS    LENGTH-1        ;INIT INPUT SAMPLES TO 0
                STF     R0,*AR1--(1)%   ;START/END AT LAST SAMPLE
                LDI     LENGTH,R4       ;COUNTER TO CALL FILT SUB
LOOP            FLOAT   *AR5,R3         ;INPUT NEW SAMPLE -> R3
                LDI     @HN_ADDR,AR0    ;"FIRST" COEFF H(N-1)->AR0
                CALL    FILT            ;GO TO SUBROUTINE FILT
                FIX     R3,R1           ;R1=INTEGER(R3)
                STI     R1,*AR6         ;OUTPUT INTEGER VALUE
                SUBI    1,R4            ;DECREMENT R4
                BNZ     LOOP            ;BRANCH UNTIL R4 <> 0
WAIT            BR      WAIT            ;WAIT
                .END

