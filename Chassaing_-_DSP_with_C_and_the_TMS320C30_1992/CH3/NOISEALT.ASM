;NOISEALT.ASM-PSEUDORANDOM NOISE GENERATOR
          .TITLE "NOISEALT"  ;PSEUDORANDOM NOISE GENERATOR
          .SECT   "VECTORS"  ;ASSEMBLE INTO VECTOR SECTION
RESET     .WORD   BEGIN      ;RESET VECTOR
          .DATA              ;ASSEMBLE INTO DATA SECTION
STACKS    .WORD   809F00H    ;INIT STACK POINTER DATA
IO_ADDR   .WORD   804000H    ;I/O ADDRESS
SEED      .WORD   7E521603h  ;INITIAL SEED VALUE
MINUS     .WORD   0FFFFF000h ;NEGATIVE LEVEL
PLUS      .WORD   1000h      ;POSITIVE LEVEL
          .TEXT              ;ASSEMBLE INTO TEXT SECTION
BEGIN     LDP     STACKS     ;INIT DATA PAGE
          LDI     @STACKS,SP ;SP-> 809F00H
          LDI     @SEED,R0   ;R0 = INITIAL SEED VALUE
LOOP      LDI     R0,R4      ;PUT SEED IN R4
          LSH     -31,R4     ;MOVE BIT 31 TO LSB    =>R4
          LDI     R0,R2      ;R2=R0=SEED
          LSH     -30,R2     ;MOVE BIT 30 TO LSB    =>R2
          ADDI    R2,R4      ;ADD BITS (31+30)      =>R4
          LDI     R0,R2      ;R2=R0=SEED
          LSH     -28,R2     ;MOVE BIT 28 TO LSB    =>R2
          ADDI    R2,R4      ;ADD BITS (31+30+28)   =>R4
          LDI     R0,R2      ;R2=R0=SEED
          LSH     -17,R2     ;MOVE BIT 17 TO LSB    =>R2
          ADDI    R2,R4      ;ADD BITS (31+30+28+17)=>R4
          AND     1,R4       ;MASK LSB OF R4
          LDIZ    @MINUS,R7  ;IF R4=0, R7 = @MINUS
          LDINZ   @PLUS,R7   ;IF R4=1, R7 = @PLUS
          LSH     1,R0       ;SHIFT SEED LEFT BY 1
          OR      R4,R0      ;PUT R4 INTO LSB OF R0 
          LDI     @IO_ADDR,AR0 ;1ST ADDR OF EXP BUS=>AR0 
          STI     R7,*AR0    ;OUTPUT R7
          BR      LOOP       ;REPEAT FOR NEXT NOISE SAMPLE
          .END               ;END



