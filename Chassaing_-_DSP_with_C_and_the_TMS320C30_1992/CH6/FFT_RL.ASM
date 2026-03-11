;FFT_RL-FROM DSP APPLICATIONS WITH THE TMS320C30, VOL.3 WITH BIT REVERSAL
;ROUTINE MODIFIED
* Name:
*    fft_rl --- radix-2 real FFT to be called as a C function.
*
* Synopsis:
*    int  fft_rl(N, M, data)
*    int  N	    FFT size: N=2**M
*    int  M	    Number of stages = log2(N)
*    float *data    Array with input and output data
*
* Description:
*    Generic function to do a radix-2 FFT computation on the 320C30.
*    The data array is N-long, with only real data.  The output is stored in
*    the same locations with real and imaginary points R and I as follows:
*	  R(0), R(1),..., R(N/2), I(N/2-1),..., I(1)
*    The program is based on the FORTRAN program in the paper by Sorensen et
*    al., June 1987 issue of Trans. on ASSP.
*    The computation is done in place, and the original data is destroyed.
*    Bit reversal is implemented at the beginning  of the function.  If this
*    is not necessary, this part can be commented out.
*    The sine/cosine table for the twiddle factors is expected to be supplied
*    during link time, and it should have the following format:
*
*	       .global	 _sine
*	       .data
*    _sine     .float	 value1 	= sin(0*2*pi/N)
*	       .float	 value2 	= sin(1*2*pi/N)
*	       ......
*	       .float	 value(N/2)    = cos((N/4-1)*2*pi/N)
*
*    The values value1 to value(N/4) are the first quarter of the sine
*    period, and value(N/4+1) to value(N/2) are the first quarter of the
*    cosine period.
*
* Stack structure upon the call:
*	       +--------------+
*    -FP(4)    |    data      |
*    -FP(3)    |    M	      |
*    -FP(2)    |    N	      |
*    -FP(1)    | return addr  |
*    -FP(0)    |   old FP     |
*	       +--------------+
*
* Registers used: R0, R1, R2, R3, R4, R5, AR0, AR1, AR2, AR4, AR5
*    IR0, IR1, RS, RE, RC
*
* AUTHOR: PANOS E. PAPAMICHALIS
*	  TEXAS INSTRUMENTS			  OCTOBER 13, 1987
*
*****************************************************************************
*
FP   .set      AR3

        .GLOBL  _fft_rl         ; ENTRY POINT FOR EXECUTION
	.GLOBL	_sine		; ADDRESS OF SINE TABLE

	.BSS   FFTSIZ,1
	.BSS   LOGFFT,1
	.BSS   INPUT,1

        .TEXT

SINTAB	  .word     _sine

;   INITIALIZE C FUNCTION

_fft_rl:  PUSH      FP             ; SAVE DEDICATED REGISTERS
	  LDI	    SP,FP
	  PUSH	    R4
	  PUSH	    R5
	  PUSH	    AR4
	  PUSH	    AR5

	  LDI	    *-FP(2),R0	   ; MOVE ARGUMENTS TO LOCATIONS MATCHING
	  STI	    R0,@FFTSIZ	   ;	THE NAMES IN THE PROGRAM
	  LDI	    *-FP(3),R0
	  STI	    R0,@LOGFFT
	  LDI	    *-FP(4),R0
	  STI	    R0,@INPUT

; DO THE BIT-REVERSING AT THE BEGINNING
;;MODIFICATIONS FROM D.G.CHANDLER/N.L.CHANG 
;;(DETAILS ON SIGNAL PROCESSING,FALL 90)
	LDI	@FFTSIZ,RC		; RC=N
	SUBI	1,RC		; RC SHOULD BE ONE LESS THAN DESIRED #
	LDI	@FFTSIZ,IR0
	LSH	-1,IR0		   ; IR0=HALF THE SIZE OF FFT=N/2
        SUBI    R0,R0   ;;(ADDED) R0=0
	LDI	R0,AR0     ;;(ADDED) AR0=0
	LDI	R0,AR1     ;;(ADDED) AR1=0
        LDI     @INPUT,IR1  ;;(ADDED),IR1 FOR PROPER OFFSET
        RPTB    BITRV
	CMPI	AR1,AR0 	   ; XCHANGE LOCATIONS ONLY
	BGE	CONT		   ;	IF AR0<AR1
	LDF	*+AR0(IR1),R0    ;;(ADDED) PRE-INCREMENT INDEX TO AR0
||	LDF	*+AR1(IR1),R1  ;;(ADDED) PRE-INCREMENT INDEX TO AR1
	STF	R0,*+AR1(IR1)    ;;(ADDED) PRE-INCREMENT INDEX TO AR1
||	STF	R1,*+AR0(IR1)  ;;(ADDED) PRE-INCREMENT INDEX TO AR0
CONT	NOP	*AR0++
BITRV	NOP	*AR1++(IR0)B

;   LENGTH-TWO BUTTERFLIES

	LDI	@INPUT,AR0	; AR0 POINTS TO X(I)
	LDI	IR0,RC		; REPEAT N/2 TIMES
	SUBI	1,RC		; RC SHOULD BE ONE LESS THAN DESIRED #

        RPTB    BLK1
	ADDF	*+AR0,*AR0++,R0 	; R0=X(I)+X(I+1)
	SUBF	*AR0,*-AR0,R1		; R1=X(I)-X(I+1)
BLK1	STF	R0,*-AR0		; X(I)=X(I)+X(I+1)
||	STF	R1,*AR0++		; X(I+1)=X(I)-X(I+1)

; FIRST PASS OF THE DO-20 LOOP (STAGE K=2 IN DO-10 LOOP)

	LDI	@INPUT,AR0	; AR0 POINTS TO X(I)
	LDI	2,IR0		; IR0=2=N2
	LDI	@FFTSIZ,RC
	LSH	-2,RC		; REPEAT N/4 TIMES
	SUBI	1,RC		; RC SHOULD BE ONE LESS THAN DESIRED #

	RPTB	BLK2
	ADDF	*+AR0(IR0),*AR0++(IR0),R0    ; R0=X(I)+X(I+2)
	SUBF	*AR0,*-AR0(IR0),R1	     ; R1=X(I)-X(I+2)
	NEGF	*+AR0,R0		     ; R0=-X(I+3)
||	STF	R0,*-AR0(IR0)		     ; X(I)=X(I)+X(I+2)
BLK2	STF	R1,*AR0++(IR0)		     ; X(I+2)=X(I)-X(I+2)
||	STF	R0,*+AR0		     ; X(I+3)=-X(I+3)

; MAIN LOOP (FFT STAGES)

	LDI	@FFTSIZ,IR0
	LSH	-2,IR0	      ; IR0=INDEX FOR E
	LDI	3,R5	      ; R5 HOLDS THE CURRENT STAGE NUMBER
	LDI	1,R4	      ; R4=N4
	LDI	2,R3	      ; R3=N2
LOOP	LSH	-1,IR0	      ; E=E/2
	LSH	1,R4	      ; N4=2*N4
	LSH	1,R3	      ; N2=2*N2

;    INNER LOOP (DO-20 LOOP IN THE PROGRAM)

	LDI	@INPUT,AR5    ; AR5 POINTS TO X(I)
INLOP	LDI	IR0,AR0
	ADDI	@SINTAB,AR0   ; AR0 POINTS TO SIN/COS TABLE
	LDI	R4,IR1	      ; IR1=N4

	LDI	AR5,AR1
	ADDI	1,AR1	      ; AR1 POINTS TO X(I1)=X(I+J)
	LDI	AR1,AR3
	ADDI	R3,AR3	      ; AR3 POINTS TO X(I3)=X(I+J+N2)
	LDI	AR3,AR2
	SUBI	2,AR2	      ; AR2 POINTS TO X(I2)=X(I-J+N2)
	ADDI	R3,AR2,AR4    ; AR4 POINTS TO X(I4)=X(I-J+N1)

	LDF	*AR5++(IR1),R0	   ; R0=X(I)
	ADDF	*+AR5(IR1),R0,R1   ; R1=X(I)+X(I+N2)
	SUBF	R0,*++AR5(IR1),R0  ; R0=-X(I)+X(I+N2)
||	STF	R1,*-AR5(IR1)	   ; X(I)=X(I)+X(I+N2)
	NEGF	R0		   ; R0=X(I)-X(I+N2)
	NEGF	*++AR5(IR1),R1	   ; R1=-X(I+N4+N2)
||	STF	R0,*AR5 	   ; X(I+N2)=X(I)-X(I+N2)
	STF	R1,*AR5 	   ; X(I+N4+N2)=-X(I+N4+N2)


; INNERMOST LOOP

	LDI	@FFTSIZ,IR1
	LSH	-2,IR1		   ; IR1=SEPARATION BETWEEN SIN/COS TBLS
	LDI	R4,RC
	SUBI	2,RC		   ; REPEAT N4-1 TIMES

	RPTB	BLK3
	MPYF	*AR3,*+AR0(IR1),R0	; R0=X(I3)*COS
	MPYF	*AR4,*AR0,R1		; R1=X(I4)*SIN
	MPYF	*AR4,*+AR0(IR1),R1	; R1=X(I4)*COS
||	ADDF	R0,R1,R2		; R2=X(I3)*COS+X(I4)*SIN
	MPYF	*AR3,*AR0++(IR0),R0	; R0=X(I3)*SIN
	SUBF	R0,R1,R0		; R0=-X(I3)*SIN+X(I4)*COS !!!
	SUBF	*AR2,R0,R1		; R1=-X(I2)+R0 !!!
	ADDF	*AR2,R0,R1		; R1=X(I2)+R0 !!!
||	STF	R1,*AR3++		; X(I3)=-X(I2)+R0 !!!
	ADDF	*AR1,R2,R1		; R1=X(I1)+R2
||	STF	R1,*AR4--		; X(I4)=X(I2)+R0 !!!
	SUBF	R2,*AR1,R1		; R1=X(I1)-R2
||	STF	R1,*AR1++		; X(I1)=X(I1)+R2
BLK3	STF	R1,*AR2--		; X(I2)=X(I1)-R2

	SUBI	@INPUT,AR5
	ADDI	R4,AR5			; AR5=I+N1
	CMPI	@FFTSIZ,AR5
	BLTD	INLOP			; LOOP BACK TO THE INNER LOOP
	ADDI	@INPUT,AR5
	NOP
	NOP

	ADDI	1,R5
	CMPI	@LOGFFT,R5
	BLE	LOOP

; RESTORE THE REGISTER VALUES AND RETURN

	  POP	    AR5
	  POP	    AR4
	  POP	    R5
	  POP	    R4
	  POP	    FP
	  RETS






