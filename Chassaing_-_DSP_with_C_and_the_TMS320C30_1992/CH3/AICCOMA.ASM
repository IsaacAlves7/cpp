*AICCOMA.ASM-AIC COMMUNICATION ROUTINES IN TMS320C30-POLLING OR INTERRUPT  
        .title  "AICCOMA"       ;TMS320C30/AIC COM WITH POLLING OR INTERRUPT
        .global AICSEC,AICSET,AICSET_I,AICIO_I,AICIO_P,TWAIT,SPSET,IOPRI,IOAUX 
        .data                   ;ASSEMBLE INTO DATA SECTION
PBASE   .word   808000h         ;PERIPHERAL BASE ADDRESS
SETSP   .word   0E970300h       ;SERIAL PORT SET-UP DATA
ATABLE  .word   AICSEC          ;SP0 AIC INIT TABLE ADDR
SPSET   .word   0h              ;SERIAL PORT 0 OFFSET
        .text                   ;ASSEMBLE INTO TEXT
AICSET  PUSH    AR0             ;SAVE AR0
        PUSH    AR1             ;SAVE AR1
        PUSH    AR7             ;SAVE AR7
        PUSH    R0              ;SAVE R0
        PUSH    R1              ;SAVE R1
        LDI     @SPSET,AR7      ;SERIAL PORT OFFSET
        LDI     @PBASE,AR0      ;AR0 -> 808000h
        ADDI    AR7,AR0         ;SERIAL PORT OFFSET + PBASE
        LDI     1h,R0           ;TIMER CLK  = H1/2 *(AIC MASTER CLK)
        STI     R0,*+AR0(28h)   ;INIT TIMER PERIOD REG(TCLK0=7.5 MHZ)
        LDI     02C1h,R0        ;INIT TIMER GLOBAL REG
        STI     R0,*+AR0(20h)   ;RESET TIMER 
        CMPI    0,AR7           ;DETERMINE REQUESTED SERIAL PORT
        BNZ     SP1A            ;BRANCH TO CONFIGURE SERIAL PORT 1
        LDI     62h,IOF         ;PRI AIC RESET = 0
        B       SP0A            ;BRANCH TO CONFIGURE SERIAL PORT 0
SP1A    LDI     26h,IOF         ;AUX AIC RESET = 0
SP0A    LDI     @ATABLE,AR1     ;AR1 -> AIC INIT DATA
        RPTS    99              ;REPEAT NEXT INSTRUCTION 100 TIMES
        NOP                     ;KEEP IOF LOW FOR A WHILE
        LDI     111h,R0         ;X & R PORT CONTROL REGS DATA
        STI     R0,*+AR0(42h)   ;FSX/DX/CLKX = SP OPERATIONAL PINS
        STI     R0,*+AR0(43h)   ;FSR/DR/CLKR = SP OPERATIONAL PINS
        LDI     @SETSP,R0       ;RESET -> SP: 16 BITS, EXT CLKS, STD MODE
        STI     R0,*+AR0(40h)   ;FSX = OUTPUT, & INT ENABL (SP GLOBAL REG)
        LDI     0,R0            ;CLEAR R0
        STI     R0,*+AR0(48h)   ;CLEAR SERIAL PORT XMIT REGISTER
        CMPI    0,AR7           ;DETERMINE REQUESTED SERIAL PORT 
        BNZ     SP1B            ;BRANCH TO CONFIGURE SERIAL PORT 1
        OR      06h,IOF         ;BRING PRI AIC OUT OF RESET
        B       SP0B            ;BRANCH TO CONFIGURE SERIAL PORT 0
SP1B    OR      60h,IOF         ;BRING AUX AIC OUT OF RESET
SP0B    LDI     03h,RC          ;RC = 3 TO TRANSMIT 4 VALUES
        RPTB    SECEND          ;REPEAT 4 DATA TRANSMITS OF SECONDARY COMM
        CALL    TWAIT           ;WAIT FOR DATA TRANSMIT
        LDI     03h,R0          ;VALUE FOR SEC XMIT REQUEST FOR AIC
        STI     R0,*+AR0(48h)   ;SECONDARY XMIT REQUEST FOR AIC
        CALL    TWAIT
        LDI     *AR1++(1),R0    ;AR1 -> NEXT AIC INIT DATA
SECEND  STI     R0,*+AR0(48h)   ;DTR =  CURRENT AIC DATA
        POP     R1              ;RESTORE R1
        POP     R0              ;RESTORE R0
        POP     AR7             ;RESTORE AR7
        POP     AR1             ;RESTORE AR1
        POP     AR0             ;RESTORE AR0
        RETS                    ;RETURN FROM SUBROUTINE
AICSET_I                        ;-----CONFIG FOR INTERRUPT-----
        PUSH    AR7             ;SAVE AR7
        CALL    AICSET          ;CALL AICSET ROUTINE
        LDI     @SPSET,AR7      ;SERIAL PORT OFFSET
        LDI     0h,IF           ;CLEAR IF REGISTER
        CMPI    0,AR7           ;DETERMINE REQUESTED SERIAL PORT
        BNZ     SP1C            ;BRANCH TO CONFIGURE SERIAL PORT 1
        OR      10h,IE          ;ENABLE EXINT0 CPU INTERRUPT
        B       SP0C            ;BRANCH TO LABEL A2
SP1C    OR      40h,IE          ;ENABLE EXINT1 CPU INTERRUPT
SP0C    OR      2000h,ST        ;GLOBAL INTERRUPT ENABLE
        POP     AR7             ;RESTORE AR7
        RETS                    ;RETURN FROM SUBROUTINE      
;--------------------TRANSMIT WAIT ROUTINE-------------------------
TWAIT   PUSH    AR0             ;SAVE AR0
        PUSH    AR7             ;SAVE AR7
        PUSH    R0              ;SAVE R0
        LDI     @PBASE,AR0      ;AR0 -> 0808000h
        LDI     @SPSET,AR7      ;SERIAL PORT OFFSET
        ADDI    AR7,AR0         ;SERIAL PORT OFFSET + PBASE
TW1     LDI     *+AR0(40h),R0   ;R0 = CONTENT OF SP GLOBAL CONTROL REG
        AND     02h,R0          ;SEE IF TRANSMIT BUFFER IS READY
        BZ      TW1             ;IF NOT READY, TRY AGAIN
        POP     R0              ;RESTORE R0
        POP     AR7             ;RESTORE AR7
        POP     AR0             ;RESTORE AR0
        RETS                    ;RETURN FROM SUBROUTINE
;--------------------AIC TRANSFER ROUTINE--------------------------
AICIO_I LSH     2,R7            ;TWO LSB MUST = 0 FOR PRIMARY AIC COM
IO      PUSH    AR0             ;SAVE AR0 
        PUSH    AR7             ;SAVE AR7
        LDI     @SPSET,AR7      ;SERIAL PORT OFFSET
        LDI     @PBASE,AR0      ;AR0 -> 0808000h
        ADDI    AR7,AR0         ;SERIAL PORT OFFSET + PBASE
        STI     R7,*+AR0(48h)   ;DTR = NEXT DATA FOR AIC D/A
        LDI     *+AR0(4Ch),R6   ;R6 = DRR DATA FROM AIC A/D
        LSH     16,R6           ;LEFT SHIFT FOR SIGN EXTENSION
        ASH     -18,R6          ;RIGHT SHIFT KEEPING SIGN
        POP     AR7             ;RESTORE AR7 
        POP     AR0             ;RESTORE AR0
        RETS
;--------------------AIC POLLING ROUTINE---------------------------
AICIO_P CALL    TWAIT           ;WAIT FOR DATA TO BE TRANSFERED
        CALL    AICIO_I         ;CALL AIC TRANSFER ROUTINE
        RETS                    ;RETURN FROM SUBROUTINE     
SW_IO   LSH     2,R7            ;PREPARE FOR SECONDARY AIC COM
        OR      03h,R7          ;SET TWO LSB FOR SECONDARY COM
        CALL    TWAIT           ;WAIT FOR DATA TO BE TRANSFERED
        CALL    IO              ;CALL AIC TRANSFER ROUTINE
        LDI     R1,R7           ;LOAD SECONDARY COM DATA INTO R7
        CALL    TWAIT           ;WAIT FOR DATA TO BE TRANSFERED
        CALL    IO              ;CALL AIC TRANSFER ROUTINE
        RETS
;SUBROUTINES FOR PRIMARY OR AUXILIARY INPUT
IOPRI   PUSH    R1              ;SAVE R1
        LDI     067h,R1         ;LOAD SECONDARY COM DATA INTO R1
        CALL    SW_IO           ;CALL IO ROUTINE TO SWITCH INPUTS
        POP     R1              ;RESTORE R1
        RETS
IOAUX   PUSH    R1              ;SAVE R1
        LDI     077h,R1         ;LOAD SECONDARY COM DATA INTO R1
        CALL    SW_IO           ;CALL IO ROUTINE TO SWITCH INPUTS
        POP     R1              ;RESTORE R1
        RETS
            







