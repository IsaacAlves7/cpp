/*LOOPALL.CMD-COMMAND FILE FOR LINKING       */
-c                 /*USING C CONVENTION      */
vecsir.obj         /*INTERR/RESET DEF        */
loopall.obj        /*MAIN PROGRAM            */
-O loopall.out     /*LINKED COFF OUTPUT FILE */
-l rts30.lib       /*RUN-TIME LIBRARY SUPPORT*/
MEMORY
{
 VECS:  org = 0        len = 0x40   /*INTERRUPT VECTORS*/   
 SRAM:  org = 0x40     len = 0x3FC0 /*USER STATIC RAM  */
 RAM:   org = 0x809800 len = 0x800  /*INTERNAL RAM     */ 
}
SECTIONS
{
 .text:   {} > SRAM      /*CODE                        */
 .cinit:  {} > RAM       /*INITIALIZATION TABLES       */
 .stack:  {} > RAM       /*SYSTEM STACK                */
 .bss:    {} > RAM       /*BSS SECTION                 */
 vecs:    {} > VECS      /*VECTOR SECTION              */
 }


