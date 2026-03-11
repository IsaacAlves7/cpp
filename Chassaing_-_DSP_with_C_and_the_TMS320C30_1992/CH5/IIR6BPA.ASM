*IIR6BPA.ASM-SIXTH-ORDER IIR BANDPASS, Fc=1250 Hz,IMPULSE LENGTH 345 
        .GLOBAL IIR               ;REF/DEF
        .TEXT                     ;ASSEMBLE INTO TEXT
BEGIN   LDI     @IO_IN,AR5        ;POINTER TO IO_INPORT
        LDI     @IO_OUT,AR6       ;POINTER TO IO_OUTPORT
        LDI     LENGTH-1,AR3      ;LENGTH OF RESPONSE
        LDI     @C_ADDR,AR0       ;AR0 POINTS TO COEFF
        LDI     @D_ADDR,AR1       ;AR1 POINTS TO DELAY
        CALL    IIR               ;CALL SUBROUTINE IIR
END     BR      END               ;WAIT
IIR     FLOAT   *AR5,R3           ;STAGE INPUT
        MPYF    *AR0++,*AR1++,R0  ;b[i][0]*dly[i][0]
        LDI     STAGES-1, RC      ;INITIALIZE STAGE COUNTER
        RPTB    LOOP              ;REPEAT LOOP RC TIMES
        MPYF    *AR0++,*AR1--,R1  ;b[i][1]*dly[i][1]
||      SUBF    R0,R3             ;input-b[i][0]*dly[i][0]
        MPYF    *AR0++,*AR1++,R0  ;a[i][1]*dly[i][0]
||      SUBF    R1,R3,R2 ;dly=input-b[i][0]*dly[i][0]-b[i][1]*dly[i][1]
        MPYF    *AR0++,*AR1--,R1  ;a[i][2]*dly[i][1]
        ADDF    R0,R1,R3          ;a[i][2]*dly[i][1]+a[i][1]*dly[i][0]
        LDF     *AR1,R4           ;dly[i][2]
||      STF     R2,*AR1++         ;dly[i][0] = dly
        MPYF    R2,*AR0++,R2      ;dly*a[i][0]
||      STF     R4,*AR1++         ;dly[i][1] = dly[i][0]
LOOP    MPYF    *AR0++,*AR1++,R0  ;b[i+1][0]*dly[i+1][0]
||      ADDF    R2,R3             ;STAGE OUTPUT=NEXT STAGE INPUT
        LDI     @C_ADDR,AR0       ;AR0 POINTS TO COEFF
        DBNZD   AR3,IIR           ;DELAYED BRANCH
        FIX     R3,R1             ;FLOAT OUTPUT TO INTEGER
        STI     R1,*AR6           ;OUTPUT TO IO_PORT
        LDI     @D_ADDR,AR1       ;AR1 POINTS TO DELAY
        RETS                      ;RETURN FROM SUBROUTINE
        .DATA          ;coeff:b[i][0],b[i][1],a[i][1],a[i][2],a[i][0]
COEFF .FLOAT  -1.4435E+00,9.4880E-01,0.0000E+00,-5.3324E-02,5.3324E-02
      .FLOAT  -1.3427E+00,8.9515E-01,0.0000E+00,-5.3324E-02,5.3324E-02
      .FLOAT  -1.3082E+00,9.4378E-01,0.0000E+00,-5.3324E-02,5.3324E-02 
DLY    .FLOAT  0.0,0.0,0.0,0.0,0.0,0.0                  
LENGTH .SET     345               ;LENGTH OF IMPULSE
STAGES .SET     3                 ;NUMBER OF STAGES
C_ADDR .WORD    COEFF             ;ADDRESS OF COEFF
D_ADDR .WORD    DLY               ;ADDRESS OF DELAY
IO_IN  .WORD    804000H           ;POINTER TO IO_INPORT
IO_OUT .WORD    804001H           ;POINTER TO IO_OUTPORT
       .SECT    "VECTORS"         ;VECTOR SECTION
MAIN   .WORD    BEGIN             ;BEGIN IN VECTOR 0H
       .END                       ;END

