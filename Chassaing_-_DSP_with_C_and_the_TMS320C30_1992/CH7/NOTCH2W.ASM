;NOTCH2W.ASM-TWO-WEIGHT ADAPTIVE NOTCH FILTER
            .TITLE "NOTCH2W"        ;2-WEIGHT ADAPTIVE NOTCH FILTER
            .GLOBAL RESET,BEGIN,FILT,ADAPT         ;REF/DEF SYMBOLS
            .SECT   "VECTORS"       ;ASSEMBLE INTO VECT SECTION       
RESET       .WORD   BEGIN           ;RESET VECTOR
            .DATA                   ;ASSEMBLE INTO DATA SECTION
SIN_ADDR    .WORD   804000H         ;SINE PORT ADDRESS
COS_ADDR    .WORD   804001H         ;COSINE PORT ADDRESS
DN_ADDR     .WORD   804002H         ;D+N PORT ADDRESS
OUT_ADDR    .WORD   804003H         ;OUTPUT PORT ADDRESS
SC_ADDR     .WORD   SC              ;RAM ADDR OF SINE+COSINE SAMPLES
COEFF_ADDR  .WORD   COEFF           ;ADDR OF COEFF H(N-1)
ERF_ADDR    .WORD   ERR_FUNC        ;ADDR OF ERROR FUNCTION
BETA        .FLOAT  5.0E-10         ;RATE OF ADAPTATION 
LENGTH      .SET    2               ;FILTER LENGTH N = 2
            .BSS    SC,LENGTH       ;N SPACE FOR SINE+COSINE SAMPLES
            .BSS    COEFF,LENGTH    ;N SPACE FOR COEFFICIENTS
            .BSS    ERR_FUNC,1      ;1 SPACE FOR ERROR FUNTION
            .TEXT                   ;ASSEMBLE INTO TEXT SECTION
BEGIN       LDI     @SIN_ADDR,AR2   ;SINE ADDR       -> AR2
            LDI     @COS_ADDR,AR5   ;COSINE ADDR     -> AR5
            LDI     @DN_ADDR,AR3    ;D+N ADDR        -> AR3
            LDI     @OUT_ADDR,AR4   ;OUTPUT ADDR     -> AR4
            LDI     @ERF_ADDR,AR6   ;ERROR FUNC ADDR -> AR6
            LDI     LENGTH,BK       ;FILTER LENGTH N -> BK
            LDI     @COEFF_ADDR,AR0 ;COEFF H(N-1) ADDRESS    -> AR0
            LDI     @SC_ADDR,AR1    ;SAMPLE SINE+COSINE ADDR -> AR1
            LDF     0,R0            ;INIT R0=0 
            RPTS    LENGTH-1        ;PERFORM NEXT 2 INST. N TIMES
            STF     R0,*AR0++       ;INIT COEFF TO 0, IN // WITH
||          STF     R0,*AR1++%      ;INIT SINE+COSINE SAMPLES TO 0
LOOP        FLOAT   *AR2,R3         ;INPUT SINE SAMPLE INTO R3
            STF     R3,*AR1++%      ;STORE SINE SAMPLE IN RAM
            FLOAT   *AR5,R3         ;INPUT COSINE SAMPLE INTO R3
            STF     R3,*AR1++%      ;STORE COSINE SAMPLE IN RAM
            FLOAT   *AR3,R4         ;INPUT SIGNAL+NOISE(D+N) INTO R4
            LDI     @COEFF_ADDR,AR0 ;H(N-1) ADDR -> AR0
            CALL    FILT            ;CALL FIR SUBROUTINE FILT
            SUBF    R0,R4,R0        ;ERROR = DN - Y --> R0
            FIX     R0,R1           ;CONVERT R0 TO INTEGER -> R1
            STI     R1,*AR4         ;STORE R1 (ERROR) INTO OUTPUT PORT
            MPYF    @BETA,R0        ;R0 = ERROR FUNCTION = BETA*ERROR
            STF     R0,*AR6         ;STORE ERROR FUNCTION
            LDI     @COEFF_ADDR,AR0 ;H(N-1) ADDR -> AR0
            CALL    ADAPT           ;CALL ADAPTING SUBROUTINE
            BR      LOOP            ;REPEAT WITH NEXT SAMPLE
;FIR FILTER SUBROUTINE
FILT        MPYF    *AR0++,*AR1++%,R0 ;H1(n)*x1(n)=y1(n) -> R0
            LDF     0,R2              ;R2 = 0
            MPYF    *AR0++,*AR1++%,R0 ;H2(n)*x2(n)=y2(n) -> R0
||          ADDF    R0,R2,R2          ;R2 = y1(n)
            ADDF    R0,R2,R0          ;y1(n)+y2(n)=y(n)  -> R0
            RETSU                     ;RETURN FROM SUBROUTINE
;ADAPTATION SUBROUTINE                
ADAPT       MPYF    *AR6,*AR1++%,R0   ;ERROR FUNCTION*x1(n) -> R0
            LDF     *AR0,R3           ;H1(n) -> R3
            MPYF    *AR6,*AR1++%,R0   ;ERROR FUNCTION*x2(n) -> R0
||          ADDF    R3,R0,R2          ;H1(n)+ERROR FUNCTION*x1(n)->R2
            LDF     *+AR0,R3          ;H2(n) -> R3
||          STF     R2,*AR0++         ;H1(n+1)=H1(n)+ERROR FUNCTION*x1(n)
            ADDF    R3,R0,R2          ;H2(n)+ERROR FUNCTION*x2(n)->R2
            STF     R2,*AR0           ;H2(n+1)=H2(n)+ERROR FUNCTION*x2(n)
            RETSU                     ;RETURN FROM SUBROUTINE
            .END                      ;END


