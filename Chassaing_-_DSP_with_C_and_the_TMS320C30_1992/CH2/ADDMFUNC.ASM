*ADDMFUNC.ASM-ASSEMBLY FUNCTION CALLED FROM C PROGRAM
FP      .set    AR3         ;FRAME POINTER IN AR3
        .global _addmfunc   ;GLOBAL REF/DEF
        .global _temp       ;GLOBAL REF/DEF
_addmfunc                   ;FUNCTION IN ASSEMBLY
        PUSH    FP          ;SAVE FP ONTO STACK
        LDI     SP,FP       ;POINT TO START OF STACK
        LDI     *-FP(2),R0  ;1ST ARGUMENT INTO R0
        ADDI    @_temp,R0   ;ADD GLOBAL VARIABLE WITH R0
        POP     FP          ;RESTORE FP
        RETS                ;RETURN FROM SUBROUTINE
