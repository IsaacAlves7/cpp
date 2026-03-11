*AIC ON THE AIB COMMUNICATIONS ROUTINE
        .TITLE  "AICBCOM"       ;AIC ON THE AIB COMMUNICATIONS ROUTINES
        .GLOBAL AICSEC,AICSET,AICIO,INITSP,INITAIC,SECEND,SEC,TWAIT ;REF/DEF 
        .OPTION X               ;FOR XREF
        .DATA                   ;ASSEMBLE INTO DATA SECTION
PBASE   .WORD   808000h         ;PERIPHERAL BASE ADDR
SETSP   .WORD   2970300h        ;SERIAL PORT SET-UP DATA
GOSP    .WORD   0C000000h       ;SERIAL PORT START DATA
ATABLE  .WORD   AICSEC          ;AIC INIT TABLE ADDR
        .TEXT                   ;START ADDR = 809800h
AICSET  PUSH    AR0             ;SAVE AR0
        PUSH    AR1             ;SAVE AR1
        PUSH    R0              ;SAVE R0
        PUSH    R1              ;SAVE R1
        LDI     @PBASE,AR0      ;AR0 -> 808000h 
        LDI     @ATABLE,AR1     ;AR1 -> AIC INIT DATA
        XOR     R1,R1           ;CLEAR R1
        STI     R1,*+AR0(58h)   ;CLEAR SERIAL1 DATA TRANSMIT REG (DTR)
INITSP  LDI     @SETSP,R0       ;RESET -> SP: 16 BITS, EXT CLKS, STD MODE
        STI     R0,*+AR0(50h)   ;FSX = OUTPUT, & INT ENABL (SP GLOBAL REG)
        LDI     131h,R0         ;X & R PORT CONTROL REGS DATA
        STI     R0,*+AR0(52h)   ;FSX/DX/CLKX = SP OPERATIONAL PINS
        STI     R0,*+AR0(53h)   ;FSR/DR/CLKR = SP OPERATIONAL PINS
        LDI     @GOSP,R0        ;R0 = VALUE TO START THE SP 
        OR      *+AR0(50h),R0   ;R0 = DATA TO END RESET OF SP
        STI     R0,*+AR0(50h)   ;START C30 SERIAL COM 
INITAIC LDI     06h,RC          ;RC = 6
        RPTB    SECEND          ;REPEAT 7 DATA TRANSMITS OF SECONDARY COMM
SEC     CALL    TWAIT           ;WAIT FOR DATA TO TRANSMIT
        LDI     *AR1++(1),R0    ;AR1 -> NEXT AIC INIT DATA 
SECEND  STI     R0,*+AR0(58h)   ;DTR =  CURRENT AIC DATA
        POP     R1              ;RESTORE R1
        POP     R0              ;RESTORE R0
        POP     AR1             ;RESTORE AR1
        POP     AR0             ;RESTORE AR0
        RETS                    ;RETURN TO MAIN PROGRAM
;AIC TRANSFER ROUTINE. D/A DATA IS TO BE IN R7 - PLACES A/D IN R6.
AICIO   PUSH    AR0             ;SAVE AR0
        PUSH    R0              ;SAVE R0
        LDI     @PBASE,AR0      ;AR0 -> 808000h
        CALL    TWAIT           ;WAIT FOR DATA TO TRANSMIT
        LDI     R7,R6           ;LOAD R7 INTO R6
        AND     -4,R6           ;MAKE SURE TWO LSBs ARE 0s
        STI     R6,*+AR0(58h)   ;DTR=NEXT DATA FOR AIC D/A
        LDI     *+AR0(5Ch),R6   ;R0=DRR DATA FROM AIC D/A
        LSH     16,R6           ;MULT INPUT BY 2**16 FOR SIGN EXT
        ASH     -16,R6          ;AND DIVIDE BY 2**16 KEEPING SIGN
        POP     R0              ;RESTORE R0
        POP     AR0             ;RESTORE AR0
        RETS                    ;RETURN FROM SUBROUTINE
;TRANSMIT WAIT ROUTINE
TWAIT   LDI     *+AR0(50h),R0   ;R0 = CONTENT OF SP GLOBAL CONTROL REG
        AND     2h,R0           ;SEE IF TRANSMIT BUFFER IS READY
        BZ      TWAIT           ;IF NOT READY, TRY AGAIN
        RETS                    ;ELSE,RETURN TO CALLING PROGRAM
        .END                    ;END


