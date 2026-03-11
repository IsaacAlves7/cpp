-c                                         /* LINK USING C CONVENTIONS      */
vecs.obj
evmsweep.obj
introuts.obj
-l rts.lib 
-o evmsweep.out
-m evmsweep.map


/* SPECIFY THE SYSTEM MEMORY MAP */

MEMORY
{
   VECS:   org = 0        len = 0x40
   SRAM:   org = 0x40     len = 0x3CF0
   RAM:    org = 0x809800 len = 0x800
}
                                     

/* SPECIFY THE SECTIONS ALLOCATION INTO MEMORY */

SECTIONS
{                                     
   .text:   {} > SRAM                /* CODE                                */
   .cinit:  {} > RAM                 /* INITIALIZATION TABLES               */
   .stack:  {} > RAM                 /* SYSTEM STACK                        */
   .bss:    {} > RAM                 /* GLOBAL & STATIC VARS (SEE NOTE 2)   */
   .sysmem: {} > RAM                 /* DYNAMIC MEMORY - DELETE IF NOT USED */
   .data:   {} > RAM
   .vecs:   {} > VECS
}




