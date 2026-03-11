;SINEASM.ASM-SINE GENERATION  Y(N)=A*Y(N-1)+B*Y(N-2)+C*X(N-1)
;                Y(N)=A*Y(N-1)+B*Y(N-2)  ,  FOR N >= 2
         .TITLE    "SINEASM"    ;SINE GENERATOR, Fd = 1 kHz
         .GLOBAL   RESET, BEGIN ;REF/DEF SYMBOLS
         .SECT     "VECTORS"    ;ASSEMBLE INTO VECTOR SECTION
RESET    .WORD     BEGIN        ;RESET VECTOR
         .DATA                  ;ASSEMBLE INTO DATA SECTION
STACKS   .WORD     809F00H      ;INIT STACK POINTER DATA
IO_ADDR  .WORD     804002H      ;OUTPUT ADDRESS
A        .FLOAT    1.618034     ;A = 2(COS wT), Fs=10 kHz
B        .FLOAT    -1.0         ;B = -1
Y1       .FLOAT    0.587785     ;INITIALLY Y(1)=C=SIN(WT)=.587785
Y0       .FLOAT    0.0          ;INITIALLY Y(0)=0
SCALER   .FLOAT    1000         ;SCALING FACTOR
         .TEXT                  ;ASSEMBLE INTO TEXT SECTION
BEGIN    LDP       STACKS       ;INIT DATA PAGE
         LDI       @STACKS,SP   ;SP-> 0809F00H
         LDI       @IO_ADDR,AR0 ;OUT ADDRESS -> AR0
         LDF       @Y0,R1       ;INITIALLY R1=Y(0)=0
         CALL      OUT          ;OUTPUT Y(0)
         LDF       @Y1,R1       ;INITIALLY R1=Y(1)
         CALL      OUT          ;OUTPUT Y(1)
         LDF       @A,R3        ;R3=A
         MPYF      R3,R1,R1     ;R1=A*Y1
         CALL      OUT          ;OUTPUT Y(2)=A*C
         LDF       @Y1,R0       ;R0=Y2(PREVIOUSLY Y1)DUE TO DELAY
         LDF       @B,R4        ;R4=B     
;Y(N) FOR N=>3                
LOOP     LDF       R1,R2        ;R2=A*Y1
         MPYF      R3,R1,R1     ;R1=A(A*Y1)
         MPYF      R4,R0,R0     ;R0=B*Y2
         ADDF      R0,R1        ;R1=OUTPUT
         CALL      OUT          ;GO TO SUB FOR OUTPUT
         LDF       R2,R0        ;R0=A*(Y1)   (FOR NEXT N)
         BR        LOOP         ;CONTINUE FOR EACH N
;OUTPUT SUBROUTINE
OUT      LDF       R1,R5        ;SAVE R1
         MPYF      @SCALER,R5   ;SCALE OUTPUT
         FIX       R5,R6        ;R6=INTEGER(R5)
         STI       R6,*AR0      ;OUTPUT @804002H
         RETS                   ;RETURN FROM SUBROUTINE
         .END                   ;END




