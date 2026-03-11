;NOISEINT.ASM-PSEUDORANDOM NOISE GENERATOR USING INTERRUPT
         .TITLE    "NOISEINT.ASM" ;INTERRUPT WITH NOISE PROGRAM
         .SECT     "VECTORS"      ;ASSEMBLE INTO VECTOR SECTION
RESET    .WORD     BEGIN          ;RESET VECTOR
         .SPACE    8              ;SKIP 8 WORDS
TIMER0   .WORD     TIM_INT        ;TINT0 VECTOR LOCATION @ 9h
         .SPACE    54             ;REMAINDER OF VECTOR SECTION
         .DATA                    ;ASSEMBLE INTO DATA SECTION
STACKS   .WORD     809F00H        ;INIT STACK POINTER DATA
PERIOD   .WORD     30H            ;INTERRUPT RATE=7.5MHz/(2*PERIOD)
IE_REG   .WORD     100H           ;ENABLE TIMER 0 (TINT0)INTERRUPT
PER_ADDR .WORD     808028H        ;(TINT0) PERIOD REG LOCATION
TCNTL    .WORD     2C1H           ;CONTROL REGISTER VALUE
ST_REG   .WORD     2000H          ;SET STATUS REG
IO_ADDR  .WORD     804000H        ;I/O ADDRESS
SEED     .WORD     7E521603H      ;INITIAL SEED VALUE
MINUS    .WORD     0FFFFF000H     ;NEGATIVE LEVEL
PLUS     .WORD     1000H          ;POSITIVE LEVEL
         .TEXT                    ;ASSEMBLE INTO TEXT SECTION
BEGIN    LDP       STACKS         ;INIT DATA PAGE
         LDI       @STACKS,SP     ;SP-> 809F00H
         LDI       @PER_ADDR,AR0  ;TINT0 PERIOD REG =>AR0
         LDI       @PERIOD,R0     ;PERIOD VALUE => R0
         STI       R0,*AR0--(8)   ;SET TINT0 PERIOD
         LDI       @TCNTL,R0      ;CONTROL REGISTER VALUE =>R0
         STI       R0,*AR0        ;SET TINT0 GLOBAL CNTRL @ 808020H
         LDI       @SEED,R0       ;R0 = INITIAL SEED VALUE
         LDI       0,R7           ;INIT R7 (OUTPUT) TO 0
         LDI       @IE_REG,IE     ;ENABLE TINT0 INTERRUPT(BIT 8)
         OR        @ST_REG,ST     ;SET STATUS REG (ENABLE GIE BIT)
WAIT     IDLE                     ;WAIT FOR INTERRUPT
         BR        WAIT           ;BRANCH TO WAIT TIL INTERRUPT
;  INTERRUPT VECTOR
TIM_INT  LDI       0,R4           ;INIT R4=0
         LDI       R0,R2          ;PUT SEED IN R2
         LSH       -17,R2         ;MOVE BIT 17 TO LSB      =>R2
         ADDI      R2,R4          ;ADD BIT (17)            =>R4
         LSH       -11,R2         ;MOVE BIT 28 TO LSB      =>R2
         ADDI      R2,R4          ;ADD BITS (28+17)        =>R4
         LSH       -2,R2          ;MOVE BIT 30 TO LSB      =>R2
         ADDI      R2,R4          ;ADD BITS (30+28+17)     =>R4
         LSH       -1,R2          ;MOVE BIT 31 TO LSB      =>R2
         ADDI      R2,R4          ;ADD BITS (31+30+28+17)  =>R4
         AND       1,R4           ;MASK LSB OF R4
         LDIZ      @MINUS,R7      ;IF R4 = 0, R7 = @MINUS
         LDINZ     @PLUS,R7       ;IF R4 = 1, R7 = @PLUS
         LSH       1,R0           ;SHIFT SEED LEFT BY 1
         OR        R4,R0          ;PUT R4 INTO LSB OF R0 
         LDI       @IO_ADDR,AR0   ;1ST ADDRESS OF EXP BUS=>AR0
         STI       R7,*AR0        ;OUTPUT NOISE SAMPLE 
         RETI                     ;RETURN FROM INTERRUPT
         .END                     ;END


