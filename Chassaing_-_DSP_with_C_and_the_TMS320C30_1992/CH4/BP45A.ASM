;BP45A.ASM-FIR BANDPASS FILTER WITH 45 COEFF,USING TMS320C30 CODE
        .TITLE  "BP45A.ASM"      ;BP FIR,Fc=1 kHz,Fs=10 kHz
        .GLOBAL FILT,MAIN,BEGIN  ;GLOBAL REF/DEF               
        .DATA                    ;ASSEMBLE INTO DATA SECTION
XN_ADDR .WORD   XN+LENGTH-1      ;(LAST) SAMPLE ADDRESS
HN_ADDR .WORD   COEFF            ;ADDR OF COEFF H(N-1)  
IO_IN   .WORD   804000H          ;INPUT PORT ADDRESS
IO_OUT  .WORD   804001H          ;OUTPUT PORT ADDRESS
XN      .USECT  "XN_BUFF",LENGTH ;BUFFER SIZE OF SAMPLES  
        .SECT   "VECTORS"        ;VECTOR SECTION
MAIN    .WORD   BEGIN            ;BEGIN @ RESET (0H)
        .TEXT                    ;ASSEMBLE INTO TEXT
BEGIN   LDI     @IO_IN,AR5       ;POINTER TO INPUT PORT ADDR
        LDI     @IO_OUT,AR6      ;POINTER TO OUTPUT PORT ADDR
        LDI     LENGTH,BK        ;SIZE OF CIRCULAR BUFFER
        LDI     @XN_ADDR,AR1     ;LAST SAMPLE ADDR -> AR1
        LDI     @COEFF,AR0       ;COEFF H(N-1) ADDR-> AR0
        LDF     0,R0             ;INITIALIZE R0
        RPTS    LENGTH-1         ;INITIALIZE ALL XN TO 0
        STF     R0,*AR1--%       ;START/END AT LAST SAMPLE
FILT    LDI     LENGTH,AR4       ;LENGTH IN AR4
        SUBI    1,AR4            ;AR4=LOOP COUNTER 
LOOP    FLOAT   *AR5,R3          ;INPUT NEW SAMPLE
        STF     R3,*AR1++        ;STORE NEWEST SAMPLE
        LDI     @HN_ADDR,AR0     ;AR0 POINTS TO COEFF H(N-1)
        LDF     0,R0             ;INITIALIZE R0
        LDF     0.0,R2           ;INITIALIZE R2
        RPTS    LENGTH-1         ;REPEAT LENGTH-1 TIMES
        MPYF    *AR0++,*AR1++%,R0  ;R0 = HN*XN 
||      ADDF    R0,R2            ;R2 IS THE ACCUMULATOR
        DBNZD   AR4,LOOP         ;DELAYED BRANCH UNTIL AR4<0
        ADDF    R0,R2            ;LAST VALUE ACCUMULATED
        FIX     R2,R0            ;FLOAT(R2) TO INTEGER(R0)
        STI     R0,*AR6          ;OUTPUT R0 TO IO_OUT
WAIT    BR      WAIT             ;WAIT 
        .DATA                    ;ASSEMBLE INTO DATA SECTION
COEFF   .FLOAT  -1.839E-3        ;H44
        .FLOAT  -2.657E-3        ;H43
        .FLOAT  -1.437E-7        ;H42
        .FLOAT   3.154E-3        ; .
        .FLOAT   2.595E-3        ; .
        .FLOAT  -4.159E-3        
        .FLOAT  -1.540E-2
        .FLOAT  -2.507E-2
        .FLOAT  -2.547E-2 
        .FLOAT  -1.179E-2
        .FLOAT   1.392E-2
        .FLOAT   4.206E-2 
        .FLOAT   5.888E-2
        .FLOAT   5.307E-2
        .FLOAT   2.225E-2
        .FLOAT  -2.410E-2
        .FLOAT  -6.754E-2
        .FLOAT  -8.831E-2
        .FLOAT  -7.475E-2
        .FLOAT  -2.956E-2
        .FLOAT   3.030E-2
        .FLOAT   8.050E-2
        .FLOAT   1.000E-1
        .FLOAT   8.050E-2
        .FLOAT   3.030E-2
        .FLOAT  -2.956E-2
        .FLOAT  -7.475E-2
        .FLOAT  -8.831E-2
        .FLOAT  -6.754E-2
        .FLOAT  -2.410E-2
        .FLOAT   2.225E-2
        .FLOAT   5.307E-2
        .FLOAT   5.888E-2
        .FLOAT   4.206E-2
        .FLOAT   1.392E-2
        .FLOAT  -1.179E-2
        .FLOAT  -2.547E-2
        .FLOAT  -2.507E-2
        .FLOAT  -1.540E-2
        .FLOAT  -4.159E-3
        .FLOAT   2.595E-3
        .FLOAT   3.154E-3
        .FLOAT  -1.437E-7
        .FLOAT  -2.657E-3        ;H1
H0      .FLOAT  -1.839E-3        ;H0
LENGTH  .SET    H0-COEFF+1       ;LENGTH = 45
        .END                     ;END

