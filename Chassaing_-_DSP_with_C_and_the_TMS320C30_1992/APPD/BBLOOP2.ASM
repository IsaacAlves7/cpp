;BBLOOP2.ASM-TMS320C30/BURR-BROWN LOOP FOR CHANNELS A AND B
          .title "BBLOOP2.ASM" ;LOOP FOR CHANNEL A AND B
          .global RESET, BEGIN, S_RATE, BB_SET, BB_IO2P
          .sect   "vectors"  ;VECTOR SECTION 
RESET     .word   BEGIN      ;RESET VECTOR
          .data              ;ASSEMBLE INTO DATA SECTION
STACKS    .word   809F00h    ;INIT STACK POINTER 
S_RATE    .word   150        ;Fs=7.5 MHz/(2*S_RATE)
          .text              ;ASSEMBLE INTO TEXT SECTION
BEGIN     LDP     STACKS     ;INIT DATA PAGE 
          LDI     @STACKS,SP ;SP -> 809F00h 
          CALL    BB_SET     ;INIT TCLK1 AND SERIAL PORT REG
          LDI     0,R4       ;CLEAR INPUT A REGISTER    
          LDI     0,R6       ;CLEAR INPUT B REGISTER
LOOP      LDI     R4,R5      ;R5=NEW INPUT A SAMPLE (OUTPUT A)
          LDI     R6,R7      ;R7=NEW INPUT B SAMPLE (OUTPUT B)
          CALL    BB_IO2P    ;OUTPUT R5 & R7 AND INPUT R4 & R6
          BR      LOOP       ;BRANCH BACK AGAIN
          .END               ;END


