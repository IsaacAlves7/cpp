;IIRMAIN.ASM-MAIN PROGRAM. CALLS IIR FILTER SUBROUTINE 
           .TITLE  "IIRMAIN"       ;CALLS HYPERSIGNAL SUB   
           .GLOBAL MAIN,BEGIN,FILT ;REF/DEF SYMBOLS
           .SECT   "VECTORS"       ;ASSEMBLE INTO VECTOR SECTION
MAIN       .WORD   BEGIN           ;RESET VECTOR
           .DATA                   ;ASSEMBLE INTO DATA SECTION
IN_ADDR    .WORD   804000h         ;INPUT PORT ADDRESS
OUT_ADDR   .WORD   804001H         ;OUTPUT PORT ADDRESS
LENGTH     .SET    345             ;LENGTH OF IMPULSE RESPONSE
           .TEXT                   ;ASSEMBLE INTO TEXT SECTION
BEGIN      LDP     IN_ADDR         ;INIT DATA PAGE
           LDI     @IN_ADDR,AR5    ;A/D ADDR -> AR5
           LDI     @OUT_ADDR,AR6   ;D/A ADDR -> AR6           
           LDI     LENGTH,R4       ;R4 = LOOP COUNTER
LOOP       FLOAT   *AR5,R3         ;INPUT NEW SAMPLE -> R3
           CALL    FILT            ;GO TO SUBROUTINE FILT   
           FIX     R3,R1           ;R1=INTEGER(R3)
           STI     R1,*AR6         ;OUTPUT INTEGER VALUE
           SUBI    1,R4            ;DECREMENT R4
           BNZ     LOOP            ;BRANCH UNTIL R4 < 0
WAIT       BR      WAIT            ;WAIT
          .END



