;ALARMR.ASM - PASS/FAIL ALARM GENERATOR USING EVM/AIC
        .TITLE  "ALARMR"      ;PASS/FAIL ALARM GEN.
        .GLOBAL RESET,BEGIN,AICSET_I,AICSEC,AICIO_I,SPSET
        .SECT   "VECTORS"     ;ASSEMBLE INTO VECTOR SECT
RESET   .WORD   BEGIN         ;RESET VECTOR
        .SPACE  4             ;SKIP 4 WORDS
TIMER0  .WORD   LOOP_S        ;TINTO VECTOR LOCATION @ 5H
        .SPACE  58            ;REMAINDER OF VECTOR SECTION
        .DATA                 ;ASSEMBLE INTO DATA SECTION
AICSEC  .WORD   1428H,1H,4A96H,67H
LENGTH  .SET    8             ;# OF REPETITIONS OF PROGRAM
TIME_P  .SET    3000          ;LENGTH OF PASS SIGNAL             
TIME_R  .SET    4             ;# OF REPS OF FAIL SIGNAL
TIME_F  .SET    1500          ;LENGTH OF FAIL SIGNAL
SEED    .WORD   7E521603H     ;INITIAL SEED VALUE
A1      .FLOAT  +1.618034     ;A COEFFICIENT FOR 1-kHz  
A2      .FLOAT  +0.618034     ;A COEFFICIENT FOR 2-kHz  
A4      .FLOAT  -1.618034     ;A COEFFICIENT FOR 4-kHz  
Y1      .FLOAT  +0.5877853    ;C COEFFICIENT FOR 1-kHz   
Y2      .FLOAT  +0.9510565    ;C COEFFICIENT FOR 2-kHz   
Y4      .FLOAT  +0.5877853    ;C COEFFICIENT FOR 4-kHz
B       .FLOAT  -1.0          ;B COEFFICIENT 
Y0      .FLOAT  0.0           ;INITIAL CONDITION
SCALER  .FLOAT  1000          ;SCALING FACTOR
        .TEXT                 ;ASSEMBLE INTO TEXT SECTION
BEGIN   LDP     SPSET         ;INIT DATA PAGE
        CALL    AICSET_I      ;INIT AIC
        LDI     @SEED,R0      ;R0 = INITIAL SEED VALUE
        PUSH    R0            ;PUSH R0 ONTO STACK
        LDI     LENGTH+1,R5   ;R5 = # OF REPS + 1 
LOOP_N  SUBI    1,R5          ;DECREMENT R5
        BZ      END           ;END AFTER LENGTH SAMPLES
        POP     R0            ;POP R0 FROM STACK
        LDI     R0,R4         ;PUT SEED IN R4
        LSH     -31,R4        ;MOVE BIT 31 TO LSB   =>R4
        LDI     R0,R2         ;R2=R0=SEED         
        LSH     -30,R2        ;MOVE BIT 30 TO LSB   =>R2
        ADDI    R2,R4         ;ADD BITS (31+30)     =>R4
        LDI     R0,R2         ;R2=R0=SEED
        LSH     -28,R2        ;MOVE BIT 28 TO LSB   =>R2
        ADDI    R2,R4         ;ADD BITS (31+30+28)  =>R4
        LDI     R0,R2         ;R2=R0=SEED
        LSH     -17,R2        ;MOVE BIT 17 TO LSB   =>R2
        ADDI    R2,R4         ;ADD BITS(31+30+28+17)=>R4
        AND     1,R4          ;MASK LSB OF R4
        LSH     1,R0          ;SHIFT SEED LEFT BY 1
        OR      R4,R0         ;PUT R4 INTO LSB OF R0 
        PUSH    R0            ;PUSH R0 ONTO STACK
        LDI     R4,R4         ;STORE INTEGER R4 
        BNZ     LOOP_P        ;TO PASS LOOP IF # 0
        BZ      LOOP_F        ;TO FAIL LOOP IF 0
;SEQUENCE FOR PASS SIGNAL => 4 kHz
LOOP_P  PUSH    R4            ;PUSH R4 ONTO STACK
        LDI     TIME_P,R6     ;LENGTH OF PASS SIGNAL        
        PUSH    R6            ;PUSH R6 ONTO STACK
        LDF     @Y0,R1        ;INITIALLY R1=Y(0)=0
        LDF     @Y4,R1        ;INITIALLY R1=Y(1)
        LDF     @A4,R3        ;R3=A
        MPYF    R3,R1,R1      ;R1=A x Y(1)
        LDF     @Y4,R0        ;R0=Y2 (PREV Y1) DUE TO DELAY
        LDF     @B,R4         ;R4=B
        BR      WAIT          ;GO TO OUTPUT ROUTINE
;SEQUENCE FOR 2-kHz FAIL SIGNAL 
LOOP_F  PUSH    R4            ;PUSH R4 ONTO STACK
        LDI     TIME_R,R6     ;# OF REPS OF FAIL SIGNAL 
        PUSH    R6            ;PUSH R6 ONTO STACK
LOOP_F2 LDI     TIME_F,R6     ;LENGTH OF FAIL SIGNAL
        LDF     @Y0,R1        ;INITIALLY R1=Y(0)=0
        LDF     @Y2,R1        ;INITIALLY R1=Y(1)
        LDF     @A2,R3        ;R3=A
        MPYF    R3,R1,R1      ;R1=A x Y(1)
        LDF     @Y2,R0        ;R0=Y2 (PREVY1) DUE TO DELAY
        LDF     @B,R4         ;R4=B
        BR      WAIT          ;GO TO OUTPUT ROUTINE
;SEQUENCE FOR 1-kHz FAIL SIGNAL
LOOP_F1 LDI     TIME_F,R6     ;LENGTH OF FAIL SIGNAL
        LDF     @Y0,R1        ;INITIALLY R1=Y(0)=0
        LDF     @Y1,R1        ;INITIALLY R1=Y(1)
        LDF     @A1,R3        ;R3=A
        MPYF    R3,R1,R1      ;R1=A x Y(1)
        LDF     @Y1,R0        ;R0=Y2 (PREV Y1) DUE TO DELAY
        LDF     @B,R4         ;R4=B
;Y(N) FOR N >= 3
WAIT    IDLE                  ;WAIT FOR INTERRUPT
        BR      WAIT          ;BRANCH TO WAIT          
LOOP_S  LDF     R1,R2         ;R2=A x Y1
        MPYF    R3,R1,R1      ;R1=A(A x Y1)
        MPYF    R4,R0,R0      ;R0=B x Y2
        ADDF    R0,R1         ;R1=OUTPUT
;OUTPUT ROUTINE                 
        PUSH    R6            ;SAVE R6
        LDF     R1,R7         ;STORE R1 INTO R7
        MPYF    @SCALER,R7    ;SCALE OUTPUT
        FIX     R7            ;CONVERT R7 INTO INTEGER
        CALL    AICIO_I       ;AIC I/O ROUTINE.OUTPUT R7
        LDF     R2,R0         ;R0=A x Y1   (FOR NEXT N) 
        POP     R6            ;RESTORE R6
        SUBI    1,R6          ;DECREMENT TIME COUNTER
        BZ      CONT          ;CONTINUE IF TIME_( ) = 0  
        RETI                  ;RETURN TO INTERRUPT
CONT    POP     R6            ;POP INTER VALUE FROM STACK
        POP     R6            ;POP NEXT STACK VALUE TO R6
        POP     R4            ;POP NEXT STACK VALUE TO R4
        BNZ     LOOP_N        ;BRANCH FOR NEXT SAMPLE
        PUSH    R4            ;PUSH R4 ONTO STACK
        SUBI    1,R6          ;DECREMENT REPS COUNTER
        LDI     R6,R1         ;LOAD R6 INTO R1
        PUSH    R6            ;PUSH R6 ONTO STACK 
        AND     1,R1          ;LOGICAL AND OF 1 & R1
        BNZ     LOOP_F1       ;GO TO 1kHz FAIL LOOP
        POP     R6            ;POP R6 FROM STACK       
        BNZD    LOOP_F2       ;BRANCH TO 2 kHz FAIL SIG.
        PUSH    R6            ;PUSH R6 ONTO STACK
        NOP                   ;NO OPERATION
        NOP                   ;NO OPERATION
        POP     R6            ;POP R6 FROM STACK
        POP     R4            ;POP R4 FROM STACK
        BZ      LOOP_N        ;GET NEXT SAMPLE
END     .END                  ;END







