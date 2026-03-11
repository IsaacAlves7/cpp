         .TITLE  "SQUARE"             ;TO GENERATE A SQUARE WAVE WITH INTERRUPT
         .SECT   "VECTORS"            ;ASSEMBLE INTO VECTOR SECTION
RESET    .WORD   BEGIN                ;RESET VECTOR
         .SPACE  8                    ;SKIP 8 WORDS
TIMER0   .WORD   TIM_INT              ;TINT0 VECTOR LOCATION @ 9h
         .SPACE  54                   ;REMAINDER OF VECTOR SECTION
         .DATA                        ;ASSEMBLE INTO DATA SECTION
STACKS   .WORD   809F00H              ;INIT STACK POINTER DATA 
PERIOD   .WORD   190H                 ;INTER RATE=8MHz/(2*PERIOD)=10 kHz
IE_REG   .WORD   100H                 ;ENABLE (TINT0) TIMER 0 INTERRUPT
PER_ADDR .WORD   808028H              ;(TINT0) TIMER 0 PERIOD REG LOCATION
TCNTL    .WORD   3C3H                 ;CONTROL REGISTER VALUE
ST_REG   .WORD   2000H                ;SET STATUS REG
IO_ADDR  .WORD   804000H              ;AIB I/O ADDRESS
         .TEXT                        ;ASSEMBLE INTO TEXT SECTION
BEGIN    LDP     STACKS               ;INIT DATA PAGE
         LDI     @STACKS,SP           ;SP-> 809F00H
         LDI     @PER_ADDR,AR0        ;TINT0 PERIOD REG =>AR0
         LDI     @PERIOD,R0           ;INTERRUPT RATE =>R0
         STI     R0,*AR0--(8)         ;SET TINT0 PERIOD
         LDI     @TCNTL,R0            ;CONTROL REGISTER VALUE =>R0
         STI     R0,*AR0              ;SET TINT0 GLOBAL CNTRL REG
         LDI     @IE_REG,IE           ;ENABLE TINT0 INTERRUPT
         OR      @ST_REG,ST           ;SET STATUS REG (ENABLE GIE BIT)
         LDI     @IO_ADDR,AR0         ;1st ADDRESS OF EXP BUS -> AR0
         LDI     7FFFH,R5             ;POSITIVE VALUE -> R5
WAIT     IDLE                         ;WAIT FOR INTERRUPT
         BR      WAIT                 ;BRANCH BACK TO WAIT
;  INTERRUPT VECTOR
TIM_INT  NEGI    R5,R7                ;R7=-R5, NEGATIVE PEAK OF WAVE
         LDI     0,R5                 ;REINITIALIZE R5
         ADDI    R7,R5                ;NEW VALUE OF R5
         LSH     16,R7                ;SHIFT OUTPUT INTO UPPER 16 BITS
         STI     R7,*+AR0(2)          ;OUTPUT TO AIB
         RETI                         ;RETURN FROM SUBROUTINE
         .END                         ;END             



