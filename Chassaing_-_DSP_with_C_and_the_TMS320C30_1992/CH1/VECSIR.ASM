;VECSIR.ASM-INTERRUPT AND RESET VECTORS 
        .asect  "vecs",0h ;interrupt and reset vectors
        .ref    _c_int00  ;compiler defined C init reset
reset:  .word   _c_int00





