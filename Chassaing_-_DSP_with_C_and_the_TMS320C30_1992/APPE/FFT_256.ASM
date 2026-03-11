;256-POINT FFT USING THE AIB AND TMS320C30 MODULE
          .tITLE "FFT_256.ASM" 
          .option X
          .global RESET,S_RATE,BEGIN,WAIT,INITAIB,AIBIN,AIBOUT,SQRT
          .global INPUT,FFT,TIM_INT
          .sect   "vectors"     ; named section in FFT_256.CMD
RESET     .word   BEGIN         ; reset vector
          .data                 ; assemble into data section
STACKS    .word   809F00h       ; Init Stack Pointer data
SRAM_CNTL .word   0808064h      ; Primary Bus Cntl. Reg. address
S_RATE    .SET    1F4h          ; Set-up conversion rate to 8KHZ(Timer 0)
IN_ADDR   .word   in_buff       ; Starting address of the input buffer
IN_PTR    .word   in_buff       ; Input buffer pointer
OUT_ADDR  .word   out_buff      ; Starting address of the output buffer
OUT_PTR   .word   out_buff      ; Output buffer pointer
SYNC_OUT  .word   08000h        ; FFT negative output sync pulse
N         .SET    256           ; 256 total number of points
SCALE     .float  0.010         ; 1 / (2 to the N - 1) or rounded to 1 / 100
SC_ADDR   .word   SCALE         ; Address of the input scaling factor
          .bss    in_buff,N     ; Make room for the input sample buffer
          .bss    out_buff,N    ; Make room for the output sample buffer
; Constants for the SQRT routine
CNST1     .SET    0.5           
CNST2     .SET    1.5
CNST3     .float  1.103553391   ; Adjusted 1.0
CNST4     .float  0.780330086   ; Adjusted SQRT(1/2)
SMSK      .word   0FF7FFFFFh
          .text
BEGIN   LDP     STACKS          ; Init Data Page 
        LDI     16,R0           ; R0 = internal 0 wait state data
        LDI     @SRAM_CNTL,AR0  ; AR0 -> Primary Bus Cntl. Reg.
        OR      *AR0,R0         ; Set the desired bits
        STI     R0,*AR0         ; Primary Bus Cntl. Reg. = R0
        LDI     @STACKS,SP      ; SP -> 0809F00h 
        LDI     @IN_ADDR,R0     ; R0 = start address of input buffer 
        STI     R0,@IN_PTR      ; Init input buffer pointer address 
        LDI     @OUT_ADDR,AR0   ; AR0 -> start address of output buffer
        STI     AR0,@OUT_PTR    ; Init output buffer pointer address
        LDI     @SYNC_OUT,R0    ; R0 = Negative sync pulse
        ADDI    N-1,AR0         ; AR0 -> last output buffer sample
        STI     R0,*AR0         ; Last sample = sync for the FFT
        CALL    INITAIB         ; Init AIB interface
        IDLE                    ; Enable and wait for Timer 0 interrupt 
N_FRAME LDI     @OUT_ADDR,R0    ; R0 = start address of output buffer
        ADDI    N,R0            ; R0 = last address of output buffer + 1
WAIT    LDI     @OUT_PTR,R1     ; R1 = Current output buffer pointer
        CMPI    R1,R0           ; See if the last FFT point is output
        BGT     WAIT            ; If not, then go back and check again
        LDI     @IN_ADDR,R0     ; else, R0 = start address of input buffer
        STI     R0,@IN_PTR      ; Init input buffer pointer address
        LDI     @OUT_ADDR,R0    ; R0 = start address of output buffer 
        STI     R0,@OUT_PTR     ; Init output buffer pointer address
        LDI     @INPUT,AR0      ; AR0 -> start addr. FFT input working buffer   
        LDI     @OUT_ADDR,AR1   ; AR1 -> start address of output buffer      
        LDI     N,IR0           ; IR0 = address offset for Bit Reversed
        LDI     *AR0++(IR0)B,R0 ; R0 = first output sample
        RPTS    N-3             ; Repeat next instruction N-2 ( 254 ) times 
        LDI     *AR0++(IR0)B,R0 ; R0 = next output sample
||      STI     R0,*AR1++       ; Out Buff. = Prev. R0, AR1 -> Out Buff. + 1 
        STI     R0,*AR1         ; Out Buff. = R0 ( output sample 255 )
        LDI     @IN_ADDR,AR0    ; AR0 -> start address of input buffer     
        LDI     @INPUT,AR1      ; AR1 -> start addr. FFT input working buffer   
        LDI     @SC_ADDR,AR2    ; AR2 -> address of scale mult. factor
        LDF     0.0,R1          ; R0 = 0.0
        LDI     2,IR0           ; IR0 = offset for FFT real values 
        LDI     N-1,RC          ; RC = N-1  ( 255 )
        FLOAT   *AR0++,R0       ; R0 = flt. pt. of first input sample 
        RPTB    W_FILL          ; Repeat next two instructions RC+1 times
        MPYF3   *AR2,R0,R0      ; R0 = scaled input sample
||      STF     R1,*+AR1        ; FFT buff. + 1 = 0.0  ( imaginary value ) 
W_FILL  FLOAT   *AR0++,R0       ; R0 = flt. pt. of next input sample 
||      STF     R0,*AR1++(IR0)  ; FFT buff. = R0, AR0 -> next real location 
        CALL    FFT             ; Compute FFT
        LDI     @INPUT,AR0      ; AR0 -> start addr. FFT input working buffer
        LDI     @INPUT,AR1      ; AR1 -> start addr. FFT input working buffer
        CALL    SQ_SUM          ; Calculate all the sum of the squares
        LDI     @INPUT,AR0      ; AR0 -> start addr. FFT input working buffer
        LDI     @INPUT,R5       ; R5 -> start addr. FFT input working buffer
        ADDI    N+N-2,R5        ; R5 -> end addr. FFT input working buffer - 1
        LDF     *AR0,R0         ; R0 = first sum of the squares
N_SQRT  CALL    SQRT            ; SQRT( sum of the squares )
        CMPI    AR0,R5          ; Does AR0 = R5
        BGTD    N_SQRT          ; If not branch after three instructions
        FIX     R0,R1           ; R1 = flt. pt of R0
        STI     R1,*AR0++(2)    ; FFT buff. = R1, AR0 = FFT buff. + 2
        LDF     *AR0,R0         ; R0 = next sum of the squares
;  Branch made here to "N_SQRT" if necessary 
        BR      N_FRAME         ; Branch for next FFT frame
; Sum of the squares routine 
; input: AR0 -> Input buff., AR1 -> Output buff. 
SQ_SUM  LDF     *AR0++,R2       ; R2 = first value of input buff.
        LDF     *AR0++,R1       ; R1 = next value of input buff.
        MPYF3   R2,R2,R0        ; R0 = first value of input buff. ** 2
        MPYF3   R1,R1,R1        ; R1 = next value of input buff. ** 2
        LDI     2,IR0           ; IR0 = 2
        LDI     N-2,RC          ; RC = N-2 ( 254 )
        RPTB    SQ_END          ; Repeat next three instructions RC+1 times
        MPYF3   *AR0,*AR0++,R0  ; R0 = input buff. ** 2
||      ADDF3   R0,R1,R2        ; R2 = prev. R0 + R1  ( sum of squares )
        LDF     *AR0,R1         ; R1 = next value of input buff.
SQ_END  MPYF    *AR0++,R1,R1    ; R1 = input buff. ** 2, AR0 -> next input buff.
||      STF     R2,*AR1++(IR0)  ; out buff. = sum of sq., AR1 -> out buff. + 2
        ADDF3   R0,R1,R2        ; R2 = last sum of squares
        STF     R2,*AR1         ; out buff. = last sum of squares
        RETS                    ; return to calling procedure
; Start of Square Root from page 297 of "Digital Signal Processing Applications
; with the TMS320 Family" Theory, Algorithms, and Implementation,  Volume 3
SQRT    LDF     R0,R3           ; Save V
; Get approximation to 1/V. For V = (1+M)*2**E
; And 0 <= M < 1, for E even : X[0] = (1-M/2)*2**-E/2
; And for E odd : X[0] = SQRT(1/2)*(1-M/2)*2**-E/2
        PUSHF   R0              ; Save V as flt. pt. V = (1+M)*2**E
        POP     R2              ; R2 <= V as integer
        XOR     @SMSK,R2        ; R2 <= complement all but sign
        LDI     R2,R1           ; R1 <= (1-M/2)*2**-E
        LDI     R2,R4           ; R4 <= R2
        LSH     8,R1            ; R1 <= R1 exp. removed
        ASH     -1,R2           ; R2 <= R2 with -E/2 exp.
        PUSH    R2              ; Save R2 as integer
        POPF    R2              ; R2 <= flt. pt.
        LDE     R2,R1           ; R1 <= (1-M/2)*2**-E/2
        LDF     @CNST3,R2       ; R2 <= 1.1... for odd E
        LSH     7,R4            ; Test LSB of E (as sign)
        LDFNN   @CNST4,R2       ; If E even R2 <= 0.78...
        MPYF    R2,R1           ; R1 <= corrected estimate
; Generate V/2 (uses MPYF)
        MPYF    CNST1,R0        ; R0 <= V/2 trunc.
        RND     R0              ; R0 <= RND V/2
; Newton iteration for Y(X) = X - V**2 = 0 ...
        MPYF    R1,R1,R2        ; R2 <= X[0]**2
        MPYF    R0,R2           ; R2 <= (V/2) * X[0]**2
        SUBRF   CNST2,R2        ; R2 <= 1.5 - (V/2) * X[0]**2
        MPYF    R2,R1           ; R1 <= X[1] = X[0] * (1.5 - (V/2)*X[0]**2)
        MPYF    R1,R1,R2        ; R2 <= X[1]**2
        MPYF    R0,R2           ; R2 <= (V/2) * X[1]**2
        SUBRF   CNST2,R2        ; R2 <= 1.5 - (V/2) * X[1]**2
        MPYF    R2,R1           ; R1 <= X[2] = X[1] * (1.5 - (V/2)*X[1]**2)
        MPYF    R1,R1,R2        ; R2 <= X[2]**2
        MPYF    R0,R2           ; R2 <= (V/2) * X[2]**2
        SUBRF   CNST2,R2        ; R2 <= 1.5 - (V/2) * X[2]**2
        MPYF    R2,R1           ; R1 <= X[3] = X[2] * (1.5 - (V/2)*X[2]**2)
        RND     R1              ; Round before *
        MPYF    R1,R1,R2        ; R2 <= X[3]**2
        RND     R2              ; Round before *
        MPYF    R0,R2           ; R2 <= (V/2) * X[3]**2
        SUBRF   CNST2,R2        ; R2 <= 1.5 - (V/2) * X[3]**2
        RND     R2              ; Round before *
        MPYF    R2,R1           ; R1 <= X[4] = X[3] * (1.5 - (V/2)*X[3]**2)
; Invert Final result and return
        POP     R2              ; R2 <= return address
        BUD     R2              ; Return (delayed)
        RND     R3              ; Round before *
        RND     R1              ; Round before *
        MPYF    R1,R3,R0        ; R0 = SQRT(V) = V*SQRT(1/V)
; Timer 0 Interrupt service routine
TIM_INT PUSH    R0              ; Save R0
        PUSH    AR0             ; Save AR0
        LDI     @OUT_PTR,AR0    ; AR0 -> current output sample pointer
        LDI     *AR0++,R0       ; R0 = curr. out sample, AR0 -> next out sample
        CALL    AIBOUT          ; Call AIB D/A output routine
        STI     AR0,@OUT_PTR    ; Save output sample pointer
        LDI     @IN_PTR,AR0     ; AR0 -> current input sample pointer 
        CALL    AIBIN           ; Call AIB A/D input routine
        STI     R0,*AR0++       ; curr. in sample = R0, AR0 -> next in sample
        STI     AR0,@IN_PTR     ; Save input sample pointer
        POP     AR0             ; Restore AR0
        POP     R0              ; Restore R0
        RETI                    ; Return from interrupt
        .END                                        


