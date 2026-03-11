;MATRIX.ASM-MATRIX MULTIPLICATION IN ASSEMBLY CODE
           .GLOBAL BEGIN
           .DATA               ;ASSEMBLE INTO DATA SECTION
A          .FLOAT  1,2,3,4,5,6,7,8,9 ;VALUES FOR MATRIX A
B          .FLOAT  1,2,3             ;VALUES FOR MATRIX B
A_ADDR     .WORD   A           ;ADDR OF MATRIX A
B_ADDR     .WORD   B           ;ADDR OF MATRIX B
IO_OUT     .WORD   804001H     ;ADDR OF OUTPUT PORT
           .TEXT               ;ASSEMBLE INTO TEXT SECTION
BEGIN      LDI     @A_ADDR,AR0 ;AR0=START ADDR OF A
           LDI     @B_ADDR,AR1 ;AR1=START ADDR OF B
           LDI     @IO_OUT,AR2 ;AR2=OUTPUT PORT ADDR
           LDI     3,R4        ;R4 IS LOOPI COUNTER
LOOPI      LDF     0,R0        ;INITIALIZE R0
           LDI     2,AR4       ;AR4 IS LOOPJ COUNTER
LOOPJ      MPYF    *AR0++,*AR1++,R1  ;A[I,J]*B[J]=R1
           ADDF    R1,R0       ;R0 ACCUMULATES RESULT
           DB      AR4,LOOPJ   ;DEC AR4 & BRANCH TIL AR4<0
           FIX     R0          ;FLOAT TO INTEGER
           STI     R0,*AR2     ;OUTPUT RESULT TO IO_OUT
           LDI     @B_ADDR,AR1 ;RELOAD START ADDR OF B
           SUBI    1,R4        ;DECREMENT R4
           BNZ     LOOPI       ;BRANCH WHILE  R4 <> 0
WAIT       BR      WAIT        ;WAIT INDEFINITELY


