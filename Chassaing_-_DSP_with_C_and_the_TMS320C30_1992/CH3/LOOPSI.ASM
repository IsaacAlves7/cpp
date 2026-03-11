;LOOPSI.ASM-LOOP PROGRAM USING INTERRUPT.CALLS AICCOMS.ASM
        .title "LOOPSI"       ;TESTS AIC
        .global RESET,BEGIN,AICSEC,AICSET,AICSET_I,AICIO_I,SPSET
        .sect   "vectors"     ;VECTOR SECTION
RESET   .word   BEGIN         ;RESET VECTOR
        .space  4             ;SKIP 4 WORDS
        .word   ISR           ;SP 0 TRANSMIT INTERRUPT SERVICE ROUTINE
        .space  58            ;REMAINDER OF VECTOR SECTION
        .data                 ;ASSEMBLE INTO DATA SECTION
STACKS  .word   809F00h       ;INIT STACK POINTER DATA
AICSEC  .word   1428h,1h,4A96h,67h    ;SP0 AIC CONFIG DATA
        .text                 ;ASSEMBLE INTO TEXT SECTION
BEGIN   LDP     STACKS        ;INIT DATA PAGE 
        LDI     @STACKS,SP    ;SP -> 809F00h 
        CALL    AICSET_I      ;INIT AIC
        LDI     0,R7          ;R7 = 0h  (OUTPUT)
LOOP    IDLE                  ;WAIT FOR TRANSMIT INTERRUPT
        LDI     R6,R7         ;R7 = NEW INPUT SAMPLE (OUTPUT)
        BR      LOOP          ;BRANCH BACK AGAIN
ISR     CALL    AICIO_I       ;OUTPUT R7 AND PLACE INPUT IN R6
        RETI                  ;RETURN FROM INTERRUPT
        .END                  ;END







