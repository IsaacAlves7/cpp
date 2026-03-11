;BBCOM.ASM-TMS320C30/BURR-BROWN COMM ROUTINES.POLLING OR INTERRUPT 
        .title "BBCOM.ASM"    ;COMM ROUTINES FOR CHANNELS A&B
        .global S_RATE,BB_SET,INITSP,BB_IO_I,BB_IO_P,BB_IO2I,BB_IO2P
        .data                 ;ASSEMBLE INTO DATA SECTION
P_CNTL  .word   808064h       ;PRIMARY BUS CONTROL REGISTER
PBASE   .word   808000h       ;PERIPHERAL BASE ADDR
SETSP   .word   2BC0000h      ;SERIAL PORT SET-UP DATA
GOSP    .word   0C000000h     ;SERIAL PORT START DATA
        .text
BB_SET  PUSH    AR0           ;SAVE AR0
        PUSH    R0            ;SAVE R0
        PUSH    R1            ;SAVE R1
        LDI     -249,R0       ;DATA FOR 0 WAIT STATE
        LDI     @P_CNTL,AR0   ;AR0->PRIMARY BUS CONTROL REG
        AND     *AR0,R0       ;AND CONTROL REG WITH 0 WAIT DATA
        STI     R0,*AR0       ;INITIATE SRAM TO 0 WAIT STATE
        LDI     @PBASE,AR0    ;AR0 -> 808000h 
        LDI     020h,IOF      ;SET FX1 LOW
INITT1  XOR     R1,R1         ;CLEAR R1
        LDI     01h, R0       ;R0 = DATA FOR TIMER 1 PERIOD REG
        STI     R1,*+AR0(38h) ;INIT TIMER 1 PERIOD REG(3.75 MHz)
        STI     R1,*+AR0(34h) ;START TIMER 1 COUNTER REG @ 0
        LDI     03C1h,R0      ;INIT TIMER 1 GLOBAL REG
        STI     R0,*+AR0(30h) ;RESET TIMER 1
        STI     R1,*+AR0(58h) ;CLEAR DATA TRANSMIT REG (DTR)
INITSP  LDI     @SETSP,R0     ;RESET->SP:32 BITS,EXT CLKS,STD MODE
        STI     R0,*+AR0(50h) ;FSX=OUTPUT,& INT ENABL(SP GLOBAL REG)
        LDI     131h,R0       ;X & R PORT CONTROL REGS DATA
        STI     R0,*+AR0(52h) ;FSX/DX/CLKX = SP OPERATIONAL PINS
        STI     R0,*+AR0(53h) ;FSR/DR/CLKR = SP OPERATIONAL PINS
        LDI     @GOSP,R0      ;R0 = VALUE TO START THE SP 
        OR      *+AR0(50h),R0 ;R0 = DATA TO END RESET OF SP
        STI     R0,*+AR0(50h) ;START C30 SERIAL COM 
INITT0  LDI     @S_RATE,R0    ;TIMER CLK 0 = 7.5 MHZ / (2 * S_RATE)
        STI     R0,*+AR0(28h) ;INIT TIMER 0 PERIOD REG (TCLK0)
        STI     R1,*+AR0(24h) ;START TIMER 0 COUNTER REG @ 0
        LDI     3C1h,R0       ;INIT TIMER 0 GLOBAL REG
        STI     R0,*+AR0(20h) ;RESET TIMER 0
        POP     R1            ;RESTORE R1
        POP     R0            ;RESTORE R0
        POP     AR0           ;RESTORE AR0
        RETS                  ;RETURN TO MAIN PROGRAM
;D/A IN R7 AND A/D IN R6. POLLS TIMER 0 INTERRUPT FLAG
BB_IO_P PUSH    AR0           ;SAVE AR0
        PUSH    R0            ;SAVE R0
WAIT1   LDI     IF,R0         ;R0=INTERRUPT FLAG REG CONTENT
        AND     100h,R0       ;SEE IF TIMER 0 INTER FLAG SET
        BZ      WAIT1         ;BRANCH TO WAIT IF NOT SET
        AND     -257,IF       ;CLEAR TIMER 0 INTERRUPT FLAG
        OR      60h,IOF       ;SET FX1 HIGH
        LDI     @PBASE,AR0    ;AR0 -> 808000h
        CALL    IN_OUT        ;OUTPUT D/A AND INPUT A/D
        POP     R0            ;RESTORE R0
        POP     AR0           ;RESTORE AR0 
        LDI     20h,IOF       ;SET FX1 LOW TO INIT A/D-D/A CONVERSION
        RETS                  ;RETURN TO MAIN PROGRAM
;D/A IN R7 AND A/D IN R6.CALLED FROM TIMER 0 INTERRUPT SERVICE ROUTINE
BB_IO_I PUSH    AR0           ;SAVE AR0
        OR      60h,IOF       ;SET FX1 HIGH
        LDI     @PBASE,AR0    ;AR0 -> 808000h
        CALL    IN_OUT        ;OUTPUT D/A AND INPUT A/D
        POP     AR0           ;RESTORE AR0 
        LDI     20h,IOF       ;SET FX1 LOW TO INIT A/D-D/A CONVERSION
        RETS                  ;RETURN TO MAIN PROGRAM
;I/O UTILITY ROUTINE 
IN_OUT  LDI     R7,R6         ;SAVE CONTENT OF THE OUTPUT REG
        LSH     14,R6         ;MULT OUTPUT BY 2**14
        STI     R6,*+AR0(58h) ;DTR=NEXT DATA FOR D/A
        LDI     *+AR0(5Ch),R6 ;R6=DRR DATA FROM A/D
        ASH     -14,R6        ;DIVIDE INPUT BY 2**14 KEEPING SIGN
        RETS                  ;RETURN TO MAIN PROGRAM
;TWO INPUTS AND OUTPUTS ROUTINES WITH 16-BIT ACCURACY
;INPUT:R5 AND R7 (A AND B D/A).   OUTPUT:R4 AND R6 (A AND B A/D)
BB_IO2P PUSH    AR0           ;SAVE AR0. POLLS TIMER 0 INTER FLAG
        PUSH    R0            ;SAVE R0
WAIT2   LDI     IF,R0         ;R0=INTERRUPT FLAG REG CONTENT
        AND     100h,R0       ;SEE IF TIMER 0 INTER FLAG IS SET
        BZ      WAIT2         ;BRANCH TO WAIT IF NOT SET
        AND     -257,IF       ;CLEAR TIMER 0 INTERRUPT FLAG
        OR      60h,IOF       ;SET FX1 HIGH
        LDI     @PBASE,AR0    ;AR0 -> 808000h
        CALL    IN_OUT2       ;OUTPUT D/A AND INPUT A/D
        POP     R0            ;RESTORE R0
        POP     AR0           ;RESTORE AR0 
        LDI     20h,IOF       ;SET FX1 LOW TO INIT A/D-D/A CONVERSION
        RETS                  ;RETURN TO MAIN PROGRAM
BB_IO2I PUSH    AR0           ;SAVE AR0.CALLED FROM TIMER 0 INTER
        OR      60h,IOF       ;SET FX1 HIGH
        LDI     @PBASE,AR0    ;AR0 -> 808000h
        CALL    IN_OUT2       ;OUTPUT D/A AND INPUT A/D
        POP     AR0           ;RESTORE AR0 
        LDI     020h,IOF      ;SET FX1 LOW TO INIT A/D-D/A CONVERSION
        RETS                  ;RETURN TO MAIN PROGRAM
;I/O  2 CHANNELS UTILITY ROUTINE 
IN_OUT2 LDI     R5,R4         ;SAVE CONTENT OF OUTPUT A
        LDI     R7,R6         ;SAVE CONTENT OF OUTPUT B
        LSH     16,R4         ;OUTPUT A IN UPPER 16(CLEAR LOWER 16)
        AND     0FFFFh,R6     ;CLEAR OUTPUT B UPPER 16 BITS
        OR      R4,R6         ;R6 = R4 (OUTPUTS A AND B)
        STI     R6,*+AR0(58h) ;DTR = NEXT DATA FOR D/A
        LDI     *+AR0(5Ch),R6 ;R6 = DRR DATA FROM A/D
        LDI     R6,R4         ;R4=R6(A IN UPPER 16 AND B IN LOWER 16)
        ASH     -16,R4        ;DIVIDE INPUT A BY 2**16 KEEPING SIGN
        LSH     16,R6         ;PUT INPUT B IN UPPER 16 BITS
        ASH     -16,R6        ;DIVIDE INPUT B BY 2**16 KEEPING SIGN
        RETS                  ;RETURN TO MAIN PROGRAM
        .END                  ;END





