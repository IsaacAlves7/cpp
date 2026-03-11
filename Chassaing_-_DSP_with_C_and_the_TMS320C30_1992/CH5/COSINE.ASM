;COSINE GENERATOR  Y(N)=A*Y(N-1)+B*Y(N-2)+X(N)-(A/2)*X(N-1)
;              Y(N)=A*Y(N-1)+B*Y(N-2)  ,  FOR N >= 2
         .TITLE    "COSINE"           ;COSINE GENERATOR PROGRAM"
         .OPTION   X                  ;FOR SYMBOLS XREF
         .GLOBAL   RESET, BEGIN, LOOP ;REF/DEF SYMBOLS
         .SECT     "VECTORS"          ;ASSEMBLE INTO VECTOR SECTION
RESET    .WORD     BEGIN              ;RESET VECTOR
         .DATA                        ;ASSEMBLE INTO DATA SECTION
STACKS   .WORD     809F00H            ;INIT STACK POINTER DATA
IO_ADDR  .WORD     804000H            ;I/O BASE ADDRESS
A        .FLOAT    1.618034           ;A = 2(COS wT), Fs=10 kHz
B        .FLOAT    -1.0               ;B = -1
YN_1     .FLOAT    0.809017           ;INITIALLY Y(1)=A-(A/2)
YN_2     .FLOAT    1.0                ;INITIALLY Y(0)=1.0
SCALER   .FLOAT    1000               ;SCALING FACTOR
         .TEXT                        ;ASSEMBLE INTO TEXT SECTION
BEGIN    LDP       STACKS             ;INIT DATA PAGE
         LDI       @STACKS,SP         ;SP-> 0809F00H
         LDI       @IO_ADDR,AR0       ;1ST ADDRESS OF EXP BUS->AR0
         LDF       @YN_2,R1           ;INITIALLY R1=Y(0)=1.0
         CALL      OUT                ;OUTPUT Y(0)
         LDF       @YN_1,R1           ;INITIALLY R1=Y(1)
         CALL      OUT                ;OUTPUT Y(1)
         LDF       @A,R3              ;R3=A
         MPYF      R3,R1,R1           ;R1=A*YN_1           
         ADDF      @B,R1              ;R1=R1+B
         CALL      OUT                ;OUTPUT Y(2)=A*(A-A/2) + B
         LDF       @YN_1,R0           ;R0=YN_2(PREVIOUSLY YN_1)DUE TO DELAY
         LDF       @B,R4              ;R4=B     
;Y(N) FOR N=>3
LOOP     LDF       R1,R2              ;R2=A*YN_1
         MPYF      R3,R1,R1           ;R1=A(A*YN_1)
         MPYF      R4,R0,R0           ;R0=B*YN_2
         ADDF      R0,R1              ;R1=A*YN_1+B*YN_2
         CALL      OUT                ;GO TO SUB FOR OUTPUT
         LDF       R2,R0              ;R0=A*Y(N_1)   (FOR NEXT N)
         BR        LOOP               ;CONTINUE FOR EACH N
;OUTPUT SUBROUTINE
OUT      LDF       R1,R5              ;SAVE R1
         MPYF      @SCALER,R5         ;SCALE OUTPUT
         FIX       R5,R6              ;R6=INTEGER(R5)
         STI       R6,*+AR0(2)        ;OUTPUT @804002H
         RETS                         ;RETURN FROM SUBROUTINE
         .END                         ;END




