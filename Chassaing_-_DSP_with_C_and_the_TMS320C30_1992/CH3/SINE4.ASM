;SINE4.ASM-GENERATES A SINE WITH 4 POINTS.USE OF INTERRUPTS
           .TITLE  "SINE4.ASM"        ;SINE GENERATION,F=Fs/4 
           .SECT   "VECTORS"          ;ASSEMBLE INTO VECTOR SECTION
RESET      .WORD   BEGIN              ;RESET VECTOR
           .SPACE  8                  ;SKIP 8 WORDS
TIMER0     .WORD   TIM_INT            ;TINT0 VECTOR LOCATION @ 9h
           .SPACE  54                 ;REMAINDER OF VECTOR SECTION   
SINE_TABLE .USECT  "SINE_BUFF",LENGTH ;SINE TABLE SIZE
           .DATA                      ;ASSEMBLE INTO DATA SECTION
STACKS     .WORD   809F00H            ;INIT STACK POINTER DATA
PERIOD     .WORD   177H               ;INTER RATE=7.5MHz/(2*PERIOD)=10 kHz
IE_REG     .WORD   100H               ;ENABLE TIMER 0 (TINT0)INTERRUPT
PER_ADDR   .WORD   808028H            ;(TINT0) PERIOD REG LOCATION
TCNTL      .WORD   2C1H               ;CONTROL REGISTER VALUE
ST_REG     .WORD   2000H              ;SET STATUS REG
IO_ADDR    .WORD   804000H            ;I/O ADDRESS
SINE_ADDR  .WORD   SINE_TABLE         ;ADDRESS OF SINE TABLE
LENGTH     .SET    4                  ;LENGTH OF CIRCULAR BUFFER
           .TEXT                      ;ASSEMBLE INTO TEXT SECTION
BEGIN      LDP     STACKS             ;INIT DATA PAGE
           LDI     @STACKS,SP         ;SP-> 809F00H
           LDI     @PER_ADDR,AR0      ;TINT0 PERIOD REG =>AR0
           LDI     @PERIOD,R0         ;PERIOD VALUE => R0
           STI     R0,*AR0--(8)       ;SET TINT0 PERIOD @ 808028H
           LDI     @TCNTL,R0          ;CONTROL REGISTER VALUE =>R0
           STI     R0,*AR0            ;SET TINT0 GLOBAL CNTRL @ 808020H   
           LDI     LENGTH,BK          ;SIZE OF BUFFER -> BK
           LDI     @SINE_ADDR,AR1     ;SINE_ADDR -> AR1
           LDI     0,R0               ;R0=0, FIRST SINE VALUE
           STI     R0,*AR1++%         ;STORE IN FIRST TABLE ADDRESS
           LDI     16384,R0           ;R0=SECOND VALUE
           STI     R0,*AR1++%         ;STORE SECOND VALUE
           LDI     0,R0               ;R0=0, THIRD SINE VALUE
           STI     R0,*AR1++%         ;STORE THIRD VALUE
           LDI     -16384,R0          ;R0=FOURTH VALUE
           STI     R0,*AR1++%         ;STORE FOURTH VALUE
           LDI     @IE_REG,IE         ;ENABLE TINT0 INTERRUPT(BIT 8)
           OR      @ST_REG,ST         ;SET STATUS REG (ENABLE GIE BIT)
WAIT       IDLE                       ;WAIT FOR INTERRUPT
           BR      WAIT               ;BRANCH TO WAIT TIL INTERRUPT
;  INTERRUPT VECTOR
TIM_INT    LDI     @IO_ADDR,AR0       ;I/O BASE ADDRESS -> AR0
           LDI     *AR1++%,R7         ;R7=TABLE VALUE
           STI     R7,*+AR0(2)        ;OUTPUT TO PORT ADDRESS 804002H
           RETI                       ;RETURN FROM INTERRUPT
           .END                       ;END


