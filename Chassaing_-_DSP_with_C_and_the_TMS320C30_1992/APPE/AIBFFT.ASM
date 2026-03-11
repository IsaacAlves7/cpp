; AIB/TMS320C30 A/D and D/A communications program, modified for FFT_256.ASM 
          .title "AIBFFT.ASM" 
          .option X
          .global TIM_INT, INITAIB, AIBIN, AIBOUT, S_RATE
          .sect   "vectors"     ; named section in FFT_256.CMD
          .space  8             ; skip 8 spaces
TIMER0    .word   TIM_INT       ; Timer 0 interrupt vector
          .space  22+32         ; remainder of vector section
          .data                 ; assemble into data section
IO_ADDR   .word   0804002h      ; Start address for AIB I/O
IE_REG    .word   100h          ; IE register data to enable Timer 0
PER_ADDR  .word   808028h       ; Timer 0 Period Register
PERIOD    .word   S_RATE        ; Sample Rate data for Timer 0 Period Register
TCNTL     .word   03C3h         ; Timer 0  control data
ST_REG    .word   02000h        ; Status Register data
          .text                 ; assemble into text section
INITAIB PUSH    AR0             ; Save AR0
        PUSH    R0              ; Save R0
        PUSH    DP              ; Save Data Page
        LDP     PER_ADDR        ; Init Data Page
        LDI     @PER_ADDR,AR0   ; AR0 -> Timer 0 Period register
        LDI     @IE_REG,IE      ; Enable INT0
        LDI     @PERIOD,R0      ; R0 = Sample Rate
        STI     R0,*AR0--(8)    ; Sample Rate to Timer 0
        LDI     @TCNTL,R0       ; R0 = Timer 0 control data
        STI     R0,*AR0         ; Load Timer 0 Control Register
        OR      @ST_REG,ST      ; Enable interrupts
        POP     DP              ; Restore Data Page
        POP     R0              ; Restore R0
        POP     AR0             ; Restore AR0
        RETS                    ; Return to calling routine
; AIB A/D input routine.   output : R0 = A/D data
AIBIN   PUSH    AR0             ; Save AR0
        PUSH    DP              ; Save Data Page
        LDP     IO_ADDR         ; Init Data Page
        LDI     @IO_ADDR,AR0    ; AR0 -> FIRST ADDRESS OF EXPANSION BUS
        LDI     *AR0,R0         ; Input A/D sample from AIB
        ASH     -16,R0          ; Put A/D data in lower 16 bits (signed)
        POP     DP              ; Restore Data Page
        POP     AR0             ; Restore AR0
        RETS                    ; Return to calling routine
; AIB II D/A output routine.  input : R0 = D/A data
AIBOUT  PUSH    AR0             ; Save AR0
        PUSH    DP              ; Save Data Page
        LDP     IO_ADDR         ; Init Data Page
        LDI     @IO_ADDR,AR0    ; AR0 -> FIRST ADDRESS OF EXPANSION BUS
        LSH     16,R0           ; Put D/A data in upper 16 bits 
        STI     R0,*AR0         ; Output D/A sample to AIB
        POP     DP              ; Restore Data Page
        POP     AR0             ; Restore AR0
        RETS                    ; Return to calling routine
        .END


