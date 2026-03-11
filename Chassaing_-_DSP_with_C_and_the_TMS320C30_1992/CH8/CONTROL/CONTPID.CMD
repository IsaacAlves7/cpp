-c        
vecs
contpid
-O contpid.out
-l rts.lib
-m contpid.map
                               

/* SPECIFY THE SYSTEM MEMORY MAP */

MEMORY
{
   VECS:   org = 0        len = 0x40
   SRAM:   org = 0x40     len = 0x3FC0
   RAM:    org = 0x809800 len = 0x800   
}

/* SPECIFY THE SECTIONS ALLOCATION INTO MEMORY */

SECTIONS
{
   .text:   {} > SRAM                 /* CODE                                */
   .cinit:  {} > RAM                  /* INITIALIZATION TABLES               */
   .stack:  {} > RAM                  /* SYSTEM STACK                        */
   .bss:    {} > RAM
   vecs:    {} > VECS
 }

