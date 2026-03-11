/*FFT8.CMD         COMMAND FILE FOR LINKING  */
-c                 /*USING C CONVENTION      */
vecsir             /*INTERR/RESET DEF        */
fft8.obj           /*MAIN PROGRAM            */
fft.obj            /*FFT function            */
-O fft8.out        /*LINKED COFF OUTPUT FILE */
-l rts30.lib       /*RUN-TIME LIBRARY SUPPORT*/
MEMORY
{
 VECS:  org = 0        len = 0x40   /*INTERRUPT VECTORS*/   
 SRAM:  org = 0x40     len = 0x3FC0 /*USER STATIC RAM  */
 RAM:   org = 0x809800 len = 0x800  /*INTERNAL RAM     */ 
 IO:    org = 0x804000 len = 0x200  /*I/O              */
}
SECTIONS
{
 .text:   {} > SRAM      /*CODE                        */
 .cinit:  {} > RAM       /*INITIALIZATION TABLES       */
 .stack:  {} > RAM       /*SYSTEM STACK                */
 vecs:    {} > VECS      /*VECTOR SECTION              */ 
 OUT_ADDR 0x804001: {} > IO     /*OUTPUT PORT ADDRESS  */         
 }


