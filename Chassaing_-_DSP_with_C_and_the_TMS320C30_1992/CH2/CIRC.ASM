;CIRC.ASM - DEMONSTRATES CIRCULAR BUFFER & DATA TRANSFER
           .TITLE "CIRC.ASM"   ;PROGRAM TITLE
           .GLOBAL BEGIN       ;REF/DEF SYMBOLS 
           .DATA               ;ASSEMBLE INTO DATA SECTION 
HN         .FLOAT  1           ;1ST VALUE IN HN ARRAY (H3)
           .FLOAT  2           ;2ND HN VALUE (H2)
           .FLOAT  0           ;3RD HN VALUE (H1)
H0         .FLOAT  3           ;4TH HN VALUE (H0)
           .BSS    X3,LENGTH   ;4 SPACES FOR XN IN BSS(INTERNAL RAM) 
XN         .FLOAT  2           ;1ST XN VALUE (X3)
           .FLOAT  3           ;2ND XN VALUE (X2)
           .FLOAT  4           ;3RD XN VALUE (X1)
X0         .FLOAT  0           ;4TH XN VALUE (X0)
H_ADDR     .WORD   HN          ;STARTING ADDR OF HN ARRAY                 
XSRAM_ADDR .WORD   XN          ;STARTING ADDR OF XN IN SRAM
XRAM_ADDR  .WORD   X3          ;START ADDR RESERVED FOR XN IN RAM
X0_ADDR    .WORD   X3+LENGTH-1 ;BOTTOM SAMPLE ADDR OF X0 IN RAM
LENGTH     .SET    H0-HN+1     ;# OF HN VALUES
;           +--------+       +------------+
; LOW ADDR  | H3 = 1 |       |    X3 = 2  | 
;           +--------+       +------------+
;           | H2 = 2 |       |    X2 = 3  | 
;           |--------+       +------------+ 
;           | H1 = 0 |       |    X1 = 4  |              
;           +--------+       +------------+
; HIGH ADDR | H0 = 3 |       |    X0 = 0  |
;           +--------+       +------------+
           .TEXT                   ;ASSEMBLE INTO TEXT SECTION
BEGIN      LDP     H_ADDR          ;INIT DATA PAGE (0)
           LDI     @XSRAM_ADDR,AR0 ;XN ADDR IN SRAM -> AR0
           LDI     @XRAM_ADDR,AR1  ;BEGINNING ADDR IN RAM 0->AR1     
           LDI     LENGTH,BK       ;BK = 4 (LENGTH OF CIRC BUFFER)
           LDI     LENGTH-1,RC     ;REPEAT COUNTER REG RC = 3
           RPTB    XFER            ;BLOCK REPEAT TIL XFER
           LDF     *AR0++,R0       ;TO TRANSFER SAMPLES XN
XFER       STF     R0,*AR1++%      ;FROM SRAM TO RAM 0(START @ 809800H)
           LDI     @H_ADDR,AR0     ;START ADDR OF HN (H3) -> AR0
           LDI     @X0_ADDR,AR1    ;BOTTOM ADDR(X0)->AR1 (LAST XN ADDR)
           LDF     6,R5            ;SET R5=6
           STF     R5,*AR1++%      ;STORE R5 @ X0 ADDR THEN AR1->X3  
           LDF     0,R0            ;INIT R0=0
           LDF     0,R2            ;INIT R2=0
           RPTS    LENGTH-1        ;EXECUTE NEXT 2 INSTR 4 TIMES
           MPYF    *AR0++,*AR1++%,R0  ;(AR0)*(AR1) -> R0
||         ADDF    R0,R2           ;IN PARALLEL WITH ACC INTO R2
           ADDF    R0,R2           ;LAST ACC -> R2=H3*X3+...+H0*X0
WAIT       BR      WAIT            ;WAIT
          .END                     ;END




