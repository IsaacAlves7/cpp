;VECS.ASM-INTERRUPT AND RESET VECTORS                         
       .asect "vecs",0h ;interrupt and reset vectors
       .ref _c_int00    ;initialization reset
       .ref _c_int05    ;serial port transmit interrupt routine
       .ref _c_int11    ;DMA interrupt
       .ref _c_int99    ;unexpected interrupt handler
reset  .word   _c_int00
int0   .word   _c_int99
int1   .word   _c_int99
int2   .word   _c_int99
int3   .word   _c_int99
xint0  .word   _c_int05
rint0  .word   _c_int99
xint1  .word   _c_int99
rint1  .word   _c_int99
tint0  .word   _c_int99
tint1  .word   _c_int99
dint   .word   _c_int11





