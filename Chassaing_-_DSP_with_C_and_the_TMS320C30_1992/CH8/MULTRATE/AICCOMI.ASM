;AICCOMI.ASM-TMS320C30/AIC COMMUNICATION PROGRAM (A DIFFERENT VERSION)

        .title  "AICCOMI.ASM"   ;TMS320C30/AIC COMMUNICATION USING INTERRUPT
        .global AICSEC,AICSET,INITT,INITSP,INITAIC,SEC,SECEND,AICIO_I,TWAIT
        .data                   ;ASSEMBLE INTO DATA SECTION
PBASE   .word   808000h         ;PERIPHERAL BASE ADDR
SETSP   .word   2970300h        ;SERIAL PORT SET-UP DATA
GOSP    .word   0C000000h       ;SERIAL PORT START DATA
ATABLE  .word   AICSEC          ;AIC INIT TABLE ADDR
        .text                   ;START ADDR = 809800h
AICSET  PUSH    AR0             ;SAVE AR0
        PUSH    AR1             ;SAVE AR1
        PUSH    R0              ;SAVE R0
        PUSH    R1              ;SAVE R1
        LDI     02h,IOF         ;AIC RESET = 0
        LDI     @PBASE,AR0      ;AR0 -> 808000h 
        LDI     @ATABLE,AR1     ;AR1 -> AIC INIT DATA
INITT   LDI     0h,R0           ;TIMER CLK 0 = H1/2 *(AIC MASTER CLK)
        STI     R0,*+AR0(28h)   ;INIT TIMER 0 PERIOD REG (TCLK0 = 7.5 or 8 MHZ)
        XOR     R1,R1           ;CLEAR R1
        STI     R1,*+AR0(24h)   ;START TIMER 0 COUNTER REG @ 0
        LDI     03C1h,R0        ;INIT TIMER 0 GLOBAL REG
        STI     R0,*+AR0(20h)   ;RESET TIMER 0
        STI     R1,*+AR0(48h)   ;CLEAR SERIAL DATA TRANSMIT REG (DTR)
INITSP  LDI     @SETSP,R0       ;RESET -> SP: 16 BITS, EXT CLKS, STD MODE
        STI     R0,*+AR0(40h)   ;FSX = OUTPUT, & INT ENABL (SP GLOBAL REG)
        LDI     131h,R0         ;X & R PORT CONTROL REGS DATA
        STI     R0,*+AR0(42h)   ;FSX/DX/CLKX = SP OPERATIONAL PINS
        STI     R0,*+AR0(43h)   ;FSR/DR/CLKR = SP OPERATIONAL PINS
        LDI     @GOSP,R0        ;R0 = VALUE TO START THE SP 
        OR      *+AR0(40h),R0   ;R0 = DATA TO END RESET OF SP
        STI     R0,*+AR0(40h)   ;START C30 SERIAL COM 
INITAIC OR      06h,IOF         ;XF0 = 1 *(AIC RESET = 1)
        LDI     06h,RC          ;RC = 6
        RPTB    SECEND          ;REPEAT 7 DATA TRANSMITS OF SECONDARY COMM
SEC     CALL    TWAIT           ;WAIT FOR DATA TRANSMIT
        LDI     *AR1++(1),R0    ;AR1 -> NEXT AIC INIT DATA 
SECEND  STI     R0,*+AR0(48h)   ;DTR =  CURRENT AIC DATA
        POP     R1              ;RESTORE R1
        POP     R0              ;RESTORE R0
        POP     AR1             ;RESTORE AR1
        POP     AR0             ;RESTORE AR0
        RETS                    ;RETURN TO MAIN PROGRAM
; ------------------- AIC TRANSFER ROUTINE -------------------------
;D/A DATA IS TO BE IN R7 - PLACES A/D IN R6. ALL OTHER REGISTERS NOT MODIFIED.
AICIO_I PUSH    AR0             ;SAVE AR0
        LDI     @PBASE,AR0      ;AR0 -> 808000h
        AND     -4,R7           ;MAKE SURE TWO LSBs ARE 0s
        STI     R7,*+AR0(48h)   ;DTR = NEXT DATA FOR AIC D/A
        LDI     *+AR0(4Ch),R6   ;R0 = DRR DATA FROM AIC A/D
        LSH     16,R6           ;MULT INPUT BY 2**16 FOR SIGN EXT
        ASH     -16,R6          ;AND DIVIDE BY 2**16 KEEPING SIGN
        POP     AR0             ;RESTORE AR0 
        RETS                    ;RETURN TO MAIN PROGRAM
; ------------------- TRANSMIT WAIT ROUTINE ------------------------
TWAIT   LDI     *+AR0(40h),R0   ;R0 = CONTENT OF SP GLOBAL CONTROL REG
        AND     02h,R0          ;SEE IF TRANSMIT BUFFER IS READY
        BZ      TWAIT           ;IF NOT READY ,TRY AGAIN
        RETS                    ;ELSE RETURN TO CALLING SECTION
        .END                    ;END






