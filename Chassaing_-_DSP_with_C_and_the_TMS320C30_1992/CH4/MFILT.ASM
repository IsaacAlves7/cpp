*MFILT.ASM-MACRO FOR FIR FILTER
MFILT   $MACRO                   ;FILTER MACRO
FILT    LDI     LENGTH,AR4       ;LENGTH IN AR4
        SUBI    1,AR4            ;AR4 = LENGTH - 1
LOOP    CALL    AICIO_P          ;AICIO ROUTINE,IN->R6 OUT->R7
        FLOAT   R6,R3            ;INPUT NEW SAMPLE
        STF     R3,*AR1++        ;STORE NEWEST SAMPLE
        LDI     @HN_ADDR,AR0     ;AR0 POINTS TO COEFF H(N-1)
        LDF     0,R0             ;INITIALIZE R0
        LDF     0,R2             ;INITIALIZE R2
        RPTS    LENGTH-1         ;REPEAT LENGTH-1 TIMES
        MPYF    *AR0++,*AR1++%,R0  ;R0 = HN*XN 
||      ADDF    R0,R2            ;R2 IS THE ACCUMULATOR
        DBNZD   AR4,LOOP         ;DELAYED BRANCH UNTIL AR4<0
        ADDF    R0,R2            ;LAST VALUE ACCUMULATED
        FIX     R2,R7            ;FLOAT(R2) TO INTEGER(R7)
        NOP                      ;ADDED DUE TO DELAY BRANCH
        BR      FILT             ;BRANCH BACK TO FILTER ROUTINE
        $ENDM                    ;END OF FILTER MACRO





