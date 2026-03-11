/*FIRGEN.CMD        COMMAND FILE                     */ 
FIRGEN.OBJ          /*GENERIC FIR PROGRAM            */
COEFGEN.OBJ         /*COEFFICIENTS FILE              */
AICCOMA.OBJ         /*AIC COMM ROUTINES IN TMS320C30 */        
-E BEGIN            /*SPECIFIES ENTRY POINT          */
-O FIRGEN.OUT       /*LINKED COFF OUTPUT FILE        */
MEMORY
{
 VECS: org = 0        len = 0x40   /*VECTOR LOCATIONS*/
 SRAM: org = 0x40     len = 0x3FC0 /*16K EXTERNAL RAM*/
 RAM : org = 0x809800 len = 0x800  /*2K INTERNAL RAM */
}
SECTIONS
{
 .data:  {} > RAM      /*DATA SECTION                */
 .text:  {} > SRAM     /*CODE                        */
 VECTORS:{} > VECS     /*RESET & INTERRUPT VECTORS   */
 XN_BUFF  ALIGN(64): {} > RAM      /*CIRCULAR BUFFER */
}
