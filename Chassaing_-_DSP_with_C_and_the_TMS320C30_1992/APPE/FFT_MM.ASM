**************************************************************************
*
*   COMPLEX, RADIX-2 FFT IN THE TMS320C30.  (MEMORY-TO-MEMORY)
*
* THE PROGRAM IS BASED ON THE FORTRAN CODE IN THE BOOK "DFT/FFT AND
* CONVOLUTION ALGORITHMS" BY BURRUS AND PARKS, P. 111.  THIS PROGRAM IS
* APPLICABLE FOR ALL FFT SIZES;  THE DIFFERENCE IS THE ACCOMPANYING FILE
* THAT CONTAINS THE SINE/COSINE TABLES AND THE VALUES OF THE VARIABLES
* N AND M.  SINCE THE TMS320C30 HAS BIT-REVERSED ADDRESSING CAPABILITY,
* IT IS ASSUMED THAT, DURING OUTPUT, THE DATA ARE READ IN BIT-REVERSED
* ORDER.  THEN, THE FFT OUTPUT WILL BE IN NORMAL ORDER.
*
* WRITTEN BY: PANOS PAPAMICHALIS
* REVISION: 1.0
* COPYRIGHT TEXAS INSTRUMENTS INC. 1988
*
**************************************************************************
*
        .GLOBL  FFT               ; ENTRY POINT FOR EXECUTION
        .GLOBL  N                 ; FFT SIZE
        .GLOBL  M                 ; LOG2(N)
        .GLOBL  SINE              ; ADDRESS OF SINE TABLE
        .GLOBL  INPUT
 
        .BSS    INP,512          ; MEMORY WITH INPUT DATA
 
        .TEXT
 
;   INITIALIZE
 
FFTSIZ  .WORD     N
LOGFFT  .WORD     M
SINTAB  .WORD     SINE
INPUT   .WORD     INP
 
FFT:    LDP       FFTSIZ          ; COMMAND TO LOAD DATA PAGE POINTER
 
        LDI     @FFTSIZ,IR1
        LSH     -2,IR1            ; IR1=N/4, POINTER FOR SIN/COS TABLE
        LDI     0,AR6             ; AR6 HOLDS THE CURRENT STAGE NUMBER
        LDI     @FFTSIZ,IR0
        LSH     1,IR0             ; IR0=2*N1 (BECAUSE OF REAL/IMAG)
        LDI     @FFTSIZ,R7        ; R7=N2
        LDI     1,AR7             ; INITIALIZE REPEAT COUNTER OF FIRST LOOP
        LDI     1,AR5             ; INITIALIZE IE INDEX (AR5=IE)
 
;   OUTER LOOP
 
LOOP:   NOP     *++AR6(1)         ; CURRENT FFT STAGE
        LDI     @INPUT,AR0        ; AR0 POINTS TO X(I)
        ADDI    R7,AR0,AR2        ; AR2 POINTS TO X(L)
        LDI     AR7,RC
        SUBI    1,RC              ; RC SHOULD BE ONE LESS THAN DESIRED #
 
; FIRST LOOP
 
        RPTB    BLK1
        ADDF    *AR0,*AR2,R0      ; R0=X(I)+X(L)
        SUBF    *AR2++,*AR0++,R1  ; R1=X(I)-X(L)
        ADDF    *AR2,*AR0,R2      ; R2=Y(I)+Y(L)
        SUBF    *AR2,*AR0,R3      ; R3=Y(I)-Y(L)
        STF     R2,*AR0--         ; Y(I)=R2  AND...
||      STF     R3,*AR2--         ; Y(L)=R3
BLK1    STF     R0,*AR0++(IR0)    ; X(I)=R0  AND...
||      STF     R1,*AR2++(IR0)    ; X(L)=R1 AND AR0,2 = AR0,2 + 2*N1
 
; IF THIS IS THE LAST STAGE, YOU ARE DONE
 
        CMPI    @LOGFFT,AR6
        BZD     END
 
; MAIN INNER LOOP
 
        LDI     4,AR1             ; INIT LOOP COUNTER FOR INNER LOOP
        LDI     2,AR0
        ADDI    @INPUT,AR0        ; INITIAL (X(I),Y(I)) POINTER
        LDI     @SINTAB,AR4       ; INITIALIZE IA INDEX (AR4=IA)
INLOP:  ADDI    AR5,AR4           ; IA=IA+IE; AR4 POINTS TO SINE
        ADDI    R7,AR0,AR2        ; (X(L),Y(L)) POINTER
        LDI     AR7,RC
        SUBI    1,RC              ; RC SHOULD BE ONE LESS THAN DESIRED #
        LDF     *AR4,R6           ; R6=SIN
 
; SECOND LOOP
 
        RPTB    BLK2
        SUBF    *AR2,*AR0,R2      ; R2=X(I)-X(L)
        SUBF    *+AR2,*+AR0,R1    ; R1=Y(I)-Y(L)
        MPYF    R2,R6,R0          ; R0=R2*SIN AND...
||      ADDF    *+AR2,*+AR0,R3    ; R3=Y(I)+Y(L)
        MPYF    R1,*+AR4(IR1),R3  ; R3=R1*COS AND...
||      STF     R3,*+AR0          ; Y(I)=Y(I)+Y(L)
        SUBF    R0,R3,R4          ; R4=R1*COS-R2*SIN
        MPYF    R1,R6,R0          ; R0=R1*SIN AND...
||      ADDF    *AR2,*AR0,R3      ; R3=X(I)+X(L)
        MPYF    R2,*+AR4(IR1),R3  ; R3=R2*COS AND...
||      STF     R3,*AR0++(IR0)    ; X(I)=X(I)+X(L) AND AR0=AR0+2*N1
        ADDF    R0,R3,R5          ; R5=R2*COS+R1*SIN
BLK2    STF     R5,*AR2++(IR0)    ; X(L)=R2*COS+R1*SIN, INCR AR2 AND..     .
||      STF     R4,*+AR2          ; Y(L)=R1*COS-R2*SIN
 
        CMPI    R7,AR1
        BNED    INLOP             ; LOOP BACK TO THE INNER LOOP (DELAYED)
        LDI     AR1,AR0
        ADDI    2,AR1             ; INCREMENT INNER LOOP COUNTER
        ADDI    @INPUT,AR0        ; NEW (X(I),Y(I)) POINTER
 
        LSH     1,AR7             ; INCREMENT LOOP COUNTER FOR NEXT TIME
 
        BRD     LOOP              ; NEXT FFT STAGE (DELAYED)
        LSH     1,AR5             ; IE=2*IE
        LDI     R7,IR0            ; N1=N2
        LSH     -1,R7             ; N2=N2/2
END     RETS  
        .END




