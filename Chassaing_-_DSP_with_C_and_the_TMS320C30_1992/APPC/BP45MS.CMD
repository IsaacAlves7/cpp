/*BP45MS.CMD    COMMAND FILE                    */ 
BP45S.OBJ     /*FIR SUB FROM HYPERSIGNAL        */
BP45M.OBJ     /*FIR MAIN PROGRAM                */
-E BEGIN      /*SPECIFIES ENTRY POINT FOR OUTPUT*/
-O BP45MS.OUT /*LINKED COFF OUTPUT FILE         */
MEMORY
{
 VECS: org = 0        len = 0x40   /*VECTOR LOCATIONS*/
 SRAM: org = 0x40     len = 0x3FC0 /*16K EXTERNAL RAM*/
 RAM : org = 0x809800 len = 0x800  /*2K INTERNAL RAM */
 IO  : org = 0x804000 len = 0x2000 /*-IOSTRB I/O     */
}
SECTIONS
{
 .data:  {} > RAM      /*DATA SECTION IN INTERNAL RAM*/
 .text:  {} > SRAM     /*CODE                        */
 .cinit: {} > SRAM     /*INITIALIZATION TABLES       */
 .stack: {} > RAM      /*SYSTEM STACK                */
 .bss:   {} > RAM      /*BSS SECTION IN RAM          */
 vecs:   {} > VECS     /*RESET & INTERRUPT VECTORS   */
 XN_BUFF  ALIGN(64): {} > RAM  /*CIRCULAR BUFFER     */
 IN_PORT  804000H  : {} > IO   /*INPUT PORT ADDRESS  */
 OUT_PORT 804001H  : {} > IO   /*OUTPUT PORT ADDRESS */
}
