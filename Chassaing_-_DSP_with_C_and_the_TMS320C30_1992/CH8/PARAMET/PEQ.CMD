/*PEQ.CMD          COMMAND FILE                    */ 
PEQ.OBJ          /*PEQ PROGRAM                     */
AICCOMA.OBJ      /*AIC COMM ROUTINES IN TMS320C30  */        
-E BEGIN         /*SPECIFIES ENTRY POINT FOR OUTPUT*/
-O PEQ.OUT       /*LINKED COFF OUTPUT FILE         */
MEMORY
{
 VECS: org = 0        len = 0x40   /*VECTOR LOCATIONS*/
 SRAM: org = 0x40     len = 0x3FC0 /*16K EXTERNAL RAM*/
 RAM : org = 0x809800 len = 0x800  /*2K INTERNAL RAM */
 IO  : org = 0x804000 len = 0x10   /*INPUT OFFSET    */
}
SECTIONS
{
 .data:  {} > RAM      /*DATA SECTION                */
 .text:  {} > SRAM     /*CODE                        */
 VECTORS:{} > VECS     /*RESET & INTERRUPT VECTORS   */
 XN_BUFF   ALIGN(512): {} > RAM  /*ALIGN CIRC BUFFER */
 COF_TABLE ALIGN(64) : {} > RAM  /*ALIGN COEF TABLE  */ 
}
