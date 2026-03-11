;BBLOOP.ASM-TMS320C30/BURR-BROWN LOOP FOR CHANNEL A
          .title "BBLOOP.ASM" ;LOOP FOR CHANNEL A
          .global RESET, BEGIN, S_RATE, BB_SET, BB_IO_P
          .sect   "vectors"  ;VECTOR SECTION 
RESET     .word   BEGIN      ;RESET VECTOR
          .data              ;ASSEMBLE INTO DATA SECTION
STACKS    .word   809F00h    ;INIT STACK POINTER 
S_RATE    .word   150        ;Fs=7.5 MHz/(2*S_RATE)
          .text              ;ASSEMBLE INTO TEXT SECTION
BEGIN     LDP     STACKS     ;INIT DATA PAGE 
          LDI     @STACKS,SP ;SP -> 809F00h 
          CALL    BB_SET     ;INIT TCLK1 AND SERIAL PORT REG
          LDI     0,R6       ;CLEAR INPUT REGISTER
LOOP      LDI     R6,R7      ;R7 = NEW INPUT SAMPLE (OUTPUT)
          CALL    BB_IO_P    ;OUTPUT R7 AND INPUT INTO R6
          BR      LOOP       ;BRANCH BACK AGAIN
          .END               ;END


