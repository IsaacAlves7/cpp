*MATRIXMF.ASM-ASSEMBLY FUNCTION CALLED FROM C PROGRAM
FP       .SET    AR3          ;FRAME POINTER IN AR3
         .GLOBAL _matrixmf    ;GLOBAL REF/DEF
_matrixmf                     ;FUNCTION IN ASSEMBLY
         PUSH    FP           ;SAVE OLD FRAME POINTER
         LDI     SP,FP        ;POINT TO START OF STACK
         PUSH    R4           ;R4 IS A DEDICATED C REGISTER
         PUSH    AR4          ;AR4 IS A DEDICATED C REGISTER
         LDI     *-FP(2),AR0  ;A ARRAY
         LDI     *-FP(3),AR1  ;B ARRAY
         LDI     *-FP(4),AR2  ;POINTER TO IO_OUTPUT
         LDI     3,R4         ;R4 IS LOOPI COUNTER
LOOPI    LDF     0,R0         ;INITIALIZE R0
         LDI     2,AR4        ;AR4 IS LOOPJ COUNTER
LOOPJ    MPYF    *AR0++,*AR1++,R1  ;A[I,J] * B[J] = R1
         ADDF    R1,R0        ;R0 ACCUMULATES RESULT
         DB      AR4,LOOPJ    ;DEC AR4 AND BRANCH TIL AR4<0
         FIX     R0           ;FLOAT TO INTEGER
         STI     R0,*AR2      ;OUTPUT RESULT TO IO_OUT
         LDI     *-FP(3),AR1  ;RELOAD START ADDR OF B ARRAY
         SUBI    1,R4         ;DECREMENT R4
         BNZ     LOOPI        ;BRANCH WHILE  R4 <> 0
         POP     AR4          ;RESTORE AR4
         POP     R4           ;RESTORE R4
         POP     FP           ;RESTORE FP
         RETS                 ;RETURN FROM SUBROUTINE

