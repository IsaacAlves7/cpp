*BP45MCFR.ASM-REAL- TIME FIR FUNCTION IN ASSEMBLY CALLED FROM BP45MCR.C
FP     .SET     AR3               ;FRAME POINTER IN AR3
       .GLOBAL  _filt             ;GLOBAL REF/DEF FILTER ROUTINE
_filt  PUSH     FP                ;SAVE FRAME POINTER
       LDI      SP,FP             ;LOAD STACK POINTER INTO FP
       PUSH     AR4               ;SAVE AR4
       PUSH     AR5               ;SAVE AR5
       PUSH     AR6               ;SAVE AR6
       LDI      *-FP(2),AR0       ;ADDR OF HN POINTER->AR0
       LDI      *-FP(3),AR1       ;ADDR OF XN POINTER->AR1
       LDI      *-FP(4),AR5       ;ADDR OF IO_INPUT POINTER->AR5
       LDI      *-FP(5),AR6       ;ADDR OF IO_OUTPUT POINTER->AR6
       LDI      *-FP(6),AR2       ;FILTER LENGTH->AR2
       LDI      AR2,BK            ;SIZE OF CIRCULAR BUFFER->BK
       SUBI     1,AR2             ;DECREMENT AR2
       ADDI     AR2,AR1           ;AR1=XN ADDR+LENGTH-1(BOTTOM)
       LDI      AR2,AR4           ;AR4 IS LOOP COUNTER
LOOP   IDLE                       ;WAIT FOR INTERRUPT
       FLOAT    *AR5,R3           ;INPUT NEW SAMPLE
       STF      R3,*AR1++%        ;STORE NEWEST SAMPLE
       LDI      *-FP(2),AR0       ;AR0 POINTS TO COEFF H(N-1)
       LDF      0,R0              ;INITIALIZE R0
       LDF      0,R2              ;INITIALIZE R2
       RPTS     AR2               ;REPEAT LENGTH-1 TIMES
       MPYF     *AR0++,*AR1++%,R0 ;R0 = HN*XN
||     ADDF     R0,R2             ;R2 = ACCUMULATOR
       DBNZD    AR4,LOOP          ;DELAYED BRANCH UNTIL AR4<0
       ADDF     R0,R2             ;LAST VALUE ACCUMULATED
       FIX      R2,R0             ;FLOAT(R2) TO INTEGER(R0)
       STI      R0,*AR6           ;OUTPUT R0 TO IO_OUTPUT
       POP      AR6               ;RESTORE THE CONTENTS OF AR6
       POP      AR5               ;RESTORE THE CONTENTS OF AR5
       POP      AR4               ;RESTORE THE CONTENTS OF AR4
       POP      FP                ;RESTORE THE FRAME POINTER
       RETS                       ;RETURN TO C PROGRAM







