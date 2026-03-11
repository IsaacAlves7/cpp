/*PREFIR.CMD      COMMAND FILE           */ 
-E BEGIN        /*SPECIFIES ENTRY POINT  */
PREFIR.OBJ      /*OBJECT CODE            */
-O PREFIR.OUT   /*LINKED COFF OUTPUT FILE*/
MEMORY
{
 VECS: org = 0        len = 0x40   /*VECTOR LOCATIONS */
 SRAM: org = 0x40     len = 0x3FC0 /*PRIMARY BUS (16K)*/
 RAM : org = 0x809800 len = 0x800  /*2K INTERNAL RAM  */
 IO  : org = 0x804000 len = 0x2000 /*8K -IOSTRB       */ 
}
SECTIONS
{
 .data:  {} > SRAM               /*DATA SECTION IN SRAM        */
 .text:  {} > SRAM               /*CODE                        */
 .cinit: {} > SRAM               /*INITIALIZATION TABLES       */
 .stack: {} > RAM                /*SYSTEM STACK                */
 .bss:   {} > RAM                /*BSS SECTION IN RAM          */
 vecs:   {} > VECS               /*RESET & INTERRUPT VECTORS   */
 IN_PORT    804000h  : {} > IO   /*INPUT PORT ADDRESS          */
 OUT_PORT   804002h  : {} > IO   /*OUTPUT PORT ADDRESS         */
 XN_BUFFER  ALIGN(64): {} > RAM  /*ALIGN CIRCULAR BUFFER IN RAM*/     
}



