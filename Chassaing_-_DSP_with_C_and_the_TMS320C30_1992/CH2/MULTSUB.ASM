;MULTSUB.ASM - SUBROUTINE FOR MULTIPLICATION PROGRAM
        .TITLE "MULTSUB.ASM"  ;MULTIPLY SUBROUTINE
        .GLOBAL MULT          ;REF/DEF SYMBOL
        .TEXT                 ;ASSEMBLE INTO TEXT SECTION
MULT    LDF     0,R0          ;INIT R0=0
        LDF     0,R2          ;INIT R2=0
        RPTS    RC            ;EXEC NEXT 2 INSTR RC+1 TIMES
        MPYF    *AR0++,*AR1++,R0    ;(AR0)*(AR1)->R0
||      ADDF    R0,R2         ;IN // WITH ACCUM INTO R2
        ADDF    R0,R2         ;LAST MULT RESULT ADDED TO R2
        RETS                  ;RETURN FROM SUBROUTINE
        .END                  ;END
