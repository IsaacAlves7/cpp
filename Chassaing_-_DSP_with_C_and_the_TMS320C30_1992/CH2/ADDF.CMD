/*ADDF.CMD                COMMAND FILE           */ 
ADDF.OBJ                /*OBJECT CODE            */
-O ADDF.OUT             /*LINKED COFF OUTPUT FILE*/
MEMORY
{
 VECS: org = 0        len = 0x40   /*VECTOR LOCATIONS */
 SRAM: org = 0x40     len = 0x3FC0 /*PRIMARY BUS (16K)*/
 RAM : org = 0x809800 len = 0x800  /*2K INTERNAL RAM  */
}
SECTIONS
{
 .data:  {} > RAM       /*DATA SECTION IN INTERNAL RAM*/
 .text:  {} > SRAM      /*CODE                        */
 .cinit: {} > SRAM      /*INITIALIZATION TABLES       */
 .stack: {} > RAM       /*SYSTEM STACK                */
 .bss:   {} > RAM       /*BSS SECTION IN RAM          */
 vecs:   {} > VECS      /*RESET & INTERRUPT VECTORS   */
}
