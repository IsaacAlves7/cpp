;BBLOOPI.ASM-TMS320C30/BURR-BROWN LOOP FOR CHANNEL A WITH INTERRUPT
          .title "BBLOOPI.ASM" ;INTERRUPT-DRIVEN LOOP FOR CHANNEL A
          .global RESET, BEGIN, S_RATE, BB_SET, BB_IO_I
          .sect   "vectors"  ;VECTOR SECTION 
RESET     .word   BEGIN      ;RESET VECTOR             
          .space  8          ;SKIP 8 WORDS
          .word   ISR        ;TIMER 0 INTERRUPT SERVICE ROUTINE
          .SPACE  54         ;REMAINDER OF VECTOR SECTION
          .data              ;ASSEMBLE INTO DATA SECTION
STACKS    .word   809F00h    ;INIT STACK POINTER 
S_RATE    .word   150        ;Fs=7.5 MHz/(2*S_RATE)
IE_REG    .word   100h       ;ENABLE TIMER 0 INTERRUPT
ST_REG    .word   2000h      ;SET STATUS REGISTER
          .text              ;ASSEMBLE INTO TEXT SECTION
BEGIN     LDP     STACKS     ;INIT DATA PAGE 
          LDI     @STACKS,SP ;SP -> 809F00h 
          CALL    BB_SET     ;INIT TCLK1 AND SERIAL PORT REG
          LDI     @IE_REG,IE ;ENABLR XINT1
          OR      @ST_REG,ST ;ENABLE GIE BIT IN STATUS REG
          LDI     0,R6       ;CLEAR INPUT REGISTER
LOOP      IDLE               ;WAIT FOR RECEIVE INTERRUPT
          LDI     R6,R7      ;R7 = NEW INPUT SAMPLE (OUTPUT)
          BR      LOOP       ;BRANCH BACK AGAIN
ISR       CALL    BB_IO_I    ;OUTPUT R7 AND INPUT INTO R6
          RETI               ;RETURN FROM INTERRUPT
          .END               ;END


