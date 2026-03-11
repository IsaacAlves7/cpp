;ADAPTER.ASM-ADAPTIVE STRUCTURE FOR NOISE CANCELLATION.OUTPUT AT E
           .TITLE "ADAPTER.ASM"    ;ADAPTIVE FIR USING AIC
           .GLOBAL RESET,BEGIN,AICSET,AICSEC,IOPRI,IOAUX,SPSET 
           .SECT   "vectors"       ;ASSEMBLE INTO VECT SECTION       
RESET      .WORD   BEGIN           ;RESET VECTOR
           .DATA                   ;ASSEMBLE INTO DATA SECTION
AICSEC     .WORD   1428h,1h,5EBEh,73h   ;FOR AIC, Fs = 8 KHz
N_ADDR     .WORD   N+LENGTH-1      ;LAST ADDR OF NOISE SAMPLES
COEFF_ADDR .WORD   COEFF           ;ADDR OF COEFF H(N-1)
ERF_ADDR   .WORD   ERR_FUNC        ;ADDR OF ERROR FUNCTION
BETA       .FLOAT  2.5E-12         ;RATE OF ADAPTATION CONSTANT 
LENGTH     .SET    50              ;FILTER LENGTH N 
           .BSS    COEFF,LENGTH    ;N SPACE FOR COEFFICIENTS
           .BSS    ERR_FUNC,1      ;1 SPACE FOR ERROR FUNTION
N          .USECT  "XN_BUFF",LENGTH ;NOISE SAMPLES BUFFER SIZE
           .TEXT                   ;ASSEMBLE INTO TEXT SECTION
BEGIN      LDP     SPSET           ;INIT DATA PAGE
           LDI     10h,R1          ;SERIAL PORT 1 OFFSET
           STI     R1,@SPSET       ;ENABLE SERIAL PORT 1
           CALL    AICSET          ;INIT AIC/TMS320C30
           LDI     @ERF_ADDR,AR6   ;ERROR FUNC ADDR       ->AR6
           LDI     LENGTH,BK       ;FILTER LENGTH N       ->BK
           LDI     @COEFF_ADDR,AR0 ;COEFF H(N-1) ADDRESS  ->AR0
           LDI     @N_ADDR,AR1     ;LAST NOISE SAMPLE ADDR->AR1
           LDF     0,R0            ;INIT R0=0 
           RPTS    LENGTH-1        ;PERFORM NEXT 2 INST. N TIMES
           STF     R0,*AR0++       ;INIT COEFF TO 0, IN // WITH
||         STF     R0,*AR1++%      ;INIT NOISE SAMPLES TO 0
LOOP       CALL    IOAUX           ;OUTPUT,THEN GET N FROM AUX IN
           FLOAT   R6,R3           ;INPUT NOISE SAMPLE INTO R3
           STF     R3,*AR1++%      ;STORE IN RAM
           LDI     @COEFF_ADDR,AR0 ;H(N-1) ADDR -> AR0
           CALL    FILT            ;CALL FIR SUBROUTINE FILT
           CALL    IOPRI           ;OUTPUT,THEN GET D+N FROM PRI IN
           FLOAT   R6,R4           ;R4=FLOATING POINT(D+N)                
           SUBF3   R0,R4,R0        ;ERROR => R0 = (D+N)-Y
           FIX     R0,R7           ;R7=INTEGER(R0)
           MPYF    @BETA,R0        ;R0=ERR FUNC=BETA*ERROR
           STF     R0,*AR6         ;STORE ERROR FUNCTION 
           LDI     LENGTH-2,RC     ;RESET REPEAT COUNTER
           LDI     @COEFF_ADDR,AR0 ;H(N-1) ADDR -> AR0
           CALL    ADAPT           ;CALL ADAPTING SUBROUTINE
           BR      LOOP            ;REPEAT WITH NEXT SAMPLE
FILT       LDF     0,R2               ;R2 = 0                 
           RPTS    LENGTH-1           ;NEXT 2 INST.(LENGTH-1)TIMES
           MPYF3   *AR0++,*AR1++%,R0  ;H(N-1-i)*x(n-(N-1-i)),i=0,...,N-1 
||         ADDF3   R0,R2,R2           ;ACCUM     
           ADDF3   R0,R2,R0           ;ADD LAST PRODUCT=y(n) -> R0
           RETSU                      ;RETURN FROM SUBROUTINE
ADAPT      MPYF3   *AR6,*AR1++%,R0 ;ERR FUNC * x(n-(N-1)) -> R0
           LDF     *AR0,R3            ;H(N-1) -> R3           
           RPTB    LOOP_END           ;REPEAT (N-2 TIMES)UNTIL LOOP_END
           MPYF3   *AR6,*AR1++%,R0 ;ERR FUNC * x(n-(N-1-i)) -> R0
||         ADDF3   R3,R0,R2           ;H(N-1-i)+ERR FUNC * x(n-(N-1-i))
LOOP_END   LDF     *+AR0(1),R3        ;LOAD SUBSEQUENT H(k) -> R3
||         STF     R2,*AR0++          ;STORE/UPDATE COEFF
           ADDF3   R3,R0,R2           ;H(n+1) = H(n) + ERR FUNC * x(n)
           STF     R2,*AR0            ;STORE/UPDATE COEFF
           RETSU                      ;RETURN FROM SUBROUTINE
           .END                       ;END




