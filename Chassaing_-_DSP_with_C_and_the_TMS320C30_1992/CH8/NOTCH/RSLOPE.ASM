;SLOPE PROGRAM TO OBTAIN COSINE FROM SINE 
                .TITLE  "RSL0PE"
                .OPTION X               ;FOR SYMBOL XREF
                .GLOBAL RSLOPE          ;REF/DEF SYMBOL
                .DATA                   ;ASSEMBLE INTO DATA SECTION
SINE_ADRS       .WORD   SIN_OLD         ;OLD SINE SAMPLE ADDRESS   
WT              .FLOAT  4.74            ;AS SAMPLING DT
                .BSS    SIN_OLD,1       ;1 SPACE FOR OLD SINE SAMPLE
                .TEXT                   ;ASSEMBLE INTO TEXT SECTION
RSLOPE          PUSH    DP              ;SAVE DATA PAGE
                PUSH    R0              ;R0=PREVIOUS RESULT
                LDP     WT              ;INIT DP
                LDI     @SINE_ADRS,AR5  ;SIN_OLD => AR5
                LDF     R6,R0           ;NEW SINE SAMPLE=>R0
                SUBF    *AR5,R6         ;NEW SAMPLE MINUS OLD SAMPLE
                MPYF    @WT,R6          ;TIMES DELTA X
                STF     R0,*AR5         ;NEW SAMPLE TO BECOME OLD => AR5 
                POP     R0              ;RETRIEVE R0
                POP     DP              ;RETRIEVE DP
                RETS                    ;RETURN FROM SUBROUTINE
                .END
               



