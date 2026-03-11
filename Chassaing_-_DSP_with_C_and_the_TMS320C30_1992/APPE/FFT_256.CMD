/* FFT_256.CMD */
FFT_256.OBJ       /*TMS320C30 256-POINT FFT PROGRAM                */
AIBFFT.OBJ        /*AIB/C30 Module A/D-D/A communication program   */
FFT_MM.OBJ        /*Texas Instruments FFT routine                  */
SIN256.OBJ        /*FFT sine values                                */
-O FFT_256.OUT    /*linked COFF output file                        */
-M FFT_256.MAP    /*memory map listing                             */
MEMORY
{
   VECS:   org = 0h       len = 0C0h    /*INTERRUPT VECTORS AND RESERVED*/
   SRAM:   org = 0C0h     len = 03F40h  /*USER STATIC RAM               */
   RAM:    org = 0809800h len = 0800h   /*INTERNAL RAM                  */
}
SECTIONS
{
   vectors:  {} > VECS                  /*INTERRUPT VECTORS             */
   .text:    {} > SRAM                  /*CODE                          */
   .data:    {} > SRAM                  /*SRAM-RESIDENT TABLES&CONSTANTS*/
   .bss:     {} > RAM                   /*BSS DATA                      */
}




