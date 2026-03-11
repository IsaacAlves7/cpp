;ADAPTDNA.ASM-ADAPTIVE FILTER WITH SIGNAL+NOISE @ PRIMARY
           .TITLE "ADAPTPDNA.ASM"  ;ADAPTIVE FIR FILTER
           .GLOBAL RESET,BEGIN,FILT,ADAPT ;REF/DEF SYMBOLS
           .SECT   "VECTORS"       ;ASSEMBLE INTO VECTOR SECTION       
RESET      .WORD   BEGIN           ;RESET VECTOR
           .DATA                   ;ASSEMBLE INTO DATA SECTION
D_ADDR     .WORD   804000H         ;DESIRED SIGNAL+NOISE PORT ADDR
DN_ADDR    .WORD   804001H         ;INPUT TO ADAPT FILTER PORT ADDR
OUT_ADDR   .WORD   804002H         ;OUTPUT PORT ADDRESS
DNS_ADDR   .WORD   DNS+LENGTH-1    ;LAST ADDR OF NOISE SAMPLES
COEFF_ADDR .WORD   COEFF           ;ADDR OF COEFF H(N-1)
ERF_ADDR   .WORD   ERR_FUNC        ;ADDR OF ERROR FUNCTION
BETA       .FLOAT  2.5E-10         ;RATE OF ADAPTATION
LENGTH     .SET    50              ;FILTER LENGTH N 
           .BSS    DNS,LENGTH      ;N SPACE FOR NOISE SAMPLES
           .BSS    COEFF,LENGTH    ;N SPACE FOR COEFFICIENTS
           .BSS    ERR_FUNC,1      ;1 SPACE FOR ERROR FUNTION
           .TEXT                   ;ASSEMBLE INTO TEXT SECTION
BEGIN      LDI     @D_ADDR,AR2     ;D+N ADDR            -> AR2
           LDI     @DN_ADDR,AR3    ;NOISE ADDRESS       -> AR3
           LDI     @OUT_ADDR,AR4   ;OUTPUT ADDR         -> AR4
           LDI     @ERF_ADDR,AR6   ;ERROR FUNC ADDR     -> AR6
           LDI     LENGTH,BK       ;FILTER LENGTH N     -> BK
           LDI     @COEFF_ADDR,AR0 ;COEFF H(N-1) ADDRESS-> AR0
           LDI     @DNS_ADDR,AR1   ;LAST N SAMPLE ADDR  -> AR1
           LDF     0,R0            ;INIT R0=0 
           RPTS    LENGTH-1        ;PERFORM NEXT 2 INST. N TIMES
           STF     R0,*AR0++       ;INIT COEFF TO 0, IN // WITH
||         STF     R0,*AR1++%      ;INIT DN SAMPLES TO 0
MAIN       LDI     LENGTH,AR5      ;AR5 = FILTER LENGTH
           SUBI    1,AR5           ;DECREMENT AR5
LOOP       FLOAT   *AR3,R3         ;INPUT NOISE SAMPLE INTO R3
           STF     R3,*AR1++(1)%   ;STORE IN RAM
           FLOAT   *AR2,R4         ;INPUT D+N SAMPLE INTO R4
           LDI     @COEFF_ADDR,AR0 ;H(N-1) ADDR -> AR0
FILT       LDF     0,R2            ;R2 = 0                 
           RPTS    LENGTH-1        ;NEXT 2 INST.(LENGTH-1)TIMES
           MPYF    *AR0++,*AR1++%,R0 ;H(N-1-i)*x(n-(N-1-i))i=0..N-1 
||         ADDF    R0,R2,R2        ;ACCUMULATE
           ADDF    R0,R2,R0        ;ADD LAST PRODUCT = y(n) -> R0 
           SUBF    R0,R4,R0        ;ERROR = (D+N)-Y -> R0
           FIX     R0,R1           ;CONVERT INTO INTEGER
           STI     R1,*AR4         ;OUTPUT
           MPYF    @BETA,R0        ;R0=ERR FUNC=BETA*ERROR
           STF     R0,*AR6         ;STORE ERROR FUNCTION 
           LDI     LENGTH-2,RC     ;RESET REPEAT COUNTER
           LDI     @COEFF_ADDR,AR0 ;H(N-1) ADDR -> AR0
ADAPT      MPYF    *AR6,*AR1++%,R0 ;ERR FUNC*x(n-(N-1))->R0
           LDF     *AR0,R3         ;H(N-1) -> R3           
           RPTB    LOOP_END        ;REPEAT(N-2 TIMES) UNTIL LOOP_END
           MPYF    *AR6,*AR1++%,R0 ;ERR FUNC * x(n-(N-1-i)) -> R0
||         ADDF    R3,R0,R2        ;H(N-1-i)+ERR FUNC * x(n-(N-1-i))
LOOP_END   LDF     *+AR0(1),R3     ;LOAD SUBSEQUENT H(k) -> R3
||         STF     R2,*AR0++       ;STORE/UPDATE COEFFICIENT
           DBNZD   AR5,LOOP        ;DELAYED BRANCH FOR NEXT n
           ADDF    R3,R0,R2        ;H(n+1) = H(n) + ERR FUNC * x(n)
           STF     R2,*AR0         ;STORE/UPDATE THE LAST COEFF
           NOP                     ;ADDED TO USE DELAYED BRANCH
           BR      MAIN            ;REPEAT FOR 50 OUTPUT POINTS
           .END                    ;END


