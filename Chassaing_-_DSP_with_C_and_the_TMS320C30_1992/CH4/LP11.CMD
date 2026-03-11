/*LP11.CMD     COMMAND FILE                     */ 
LP11.OBJ      /*FIR PROGRAM                     */
-E BEGIN      /*SPECIFIES ENTRY POINT FOR OUTPUT*/
-O LP11.OUT   /*LINKED COFF OUTPUT FILE         */
MEMORY
{
 VECS: org = 0        len = 0x40   /*VECTOR LOCATIONS*/
 SRAM: org = 0x40     len = 0x3FC0 /*16K EXTERNAL RAM*/
 RAM : org = 0x809800 len = 0x800  /*2K INTERNAL RAM */
}
SECTIONS
{
 .data:  {} > SRAM     /*DATA SECTION                */
 .text:  {} > SRAM     /*CODE                        */
 .cinit: {} > SRAM     /*INITIALIZATION TABLES       */
 .stack: {} > RAM      /*SYSTEM STACK                */
 .bss:   {} > RAM      /*BSS SECTION IN RAM          */
 VECTORS:{} > VECS     /*RESET & INTERRUPT VECTORS   */
 XN_BUFF  ALIGN(64): {} > RAM  /*CIRCULAR BUFFER     */
}
