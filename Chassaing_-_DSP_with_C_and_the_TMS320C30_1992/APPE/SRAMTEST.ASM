;This Program tests the external SRAM on the TMS320C30 module.
;The address locations are 0h - 03FFFh. First the locations are tested with
;ones, and then with zeros. If the tests are successful, R0 and R1 will end 
;up with  zeros. If these registers are not zero, then the tests are not 
;successful and the error can be found in (AR0 - 1)
                .TITLE  "SRAMTEST"      ;ADDRESSES TO TEST : 0 - 3FFFh
                .OPTION X               ;FOR XREF SYMBOLS
                .GLOBAL RESET, BEGIN    ;REF/DEF SYMBOLS
                .SECT   "VECTORS"       ;ASSEMBLE INTO VECTOR SECTION
RESET           .WORD   BEGIN           ;RESET VECTOR
                .DATA                   ;ASSEMBLE INTO DATA SECTION
FIRST_ADDR      .WORD   0h              ;ADDRESS OF FIRST SRAM LOCATION
MID_ADDR        .WORD   2000h           ;ADDRESS OF MIDDLE SRAM LOCATION
SRAM_CNTL       .WORD   808064h         ;PRIMARY BUS CONTROL REGISTER
                .TEXT                   ;ASSEMBLE INTO TEXT SECTION
BEGIN           LDP     FIRST_ADDR      ;INIT DATA PAGE 
                LDI     -225,R0         ;INTERNAL 0 WAIT STATE(-193 FOR WS=1)
                LDI     @SRAM_CNTL,AR0  ;PRIMARY BUS CNTL REG => AR0
                AND     *AR0,R0         ;SET DESIRED BITS
                STI     R0,*AR0         ;REPLACE CONTENTS PRIMARY BUS CNTL REG
                LDI     @FIRST_ADDR,AR0 ;FIRST SRAM ADDRESS => AR0
                LDI     @MID_ADDR,AR1   ;MIDDLE ADDRESS OF SRAM => AR1
                LDI     -1,R0           ;INITIAL VALUES ONES
                RPTS    1FFFh           ;INITIALIZE ALL SRAM
                STI     R0,*AR0++(1)    ;
||              STI     R0,*AR1++(1)    ;
                LDI     @FIRST_ADDR,AR0 ;FIRST SRAM ADDRESS => AR0
                LDI     4000h, R7       ;R7 = TOTAL # OF SRAM ADDRESSES
FIRST_T         LDI     *AR0++(1),R1    ;R1 = CONTENTS OF SRAM
                XOR     R0, R1          ;R1 = 0 IF R0=R1
                BNZ     END_T           ;IF NOT ZERO,END TESTING WITH ERROR
                SUBI    1, R7           ;ELSE, DECREMENT COUNTER
                BNZ     FIRST_T         ;BRANCH BACK IF COUNTER NOT ZERO
                LDI     @FIRST_ADDR,AR0 ;ELSE,FIRST SRAM ADDRESS => AR0
                LDI     @MID_ADDR,AR1   ;MIDDLE ADDRESS => AR1
                LDI     0,R0            ;INITIAL VALUES ZEROS
                RPTS    1FFFh           ;INITIALIZE ALL SRAM
                STI     R0,*AR0++(1)    ;
||              STI     R0,*AR1++(1)    ;
                LDI     @FIRST_ADDR,AR0 ;FIRST SRAM ADDRESS => AR0
                LDI     4000h, R7       ;R7 = TOTAL # OF SRAM ADDRESSES
LAST_T          LDI     *AR0++(1),R1    ;R1 = CONTENTS OF SRAM
                XOR     R0, R1          ;R1 = 0 IF R0=R1
                BNZ     END_T           ;IF NOT ZERO,END TESTING WITH ERROR
                SUBI    1, R7           ;ELSE,DECREMENT COUNTER
                BNZ     LAST_T          ;BRANCH BACK IF COUNTER NOT ZERO
END_T           SWI                     ;EMU30 SOFTWARE INTERRUPT
                .END                    ;END


