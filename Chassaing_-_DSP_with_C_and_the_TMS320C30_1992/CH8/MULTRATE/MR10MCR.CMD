/*MR10MCR.CMD              REAL-TIME MULTIRATE FILTER COMMAND FILE           */
mr10mcr.obj                             /* MAIN PROGRAM                      */
aiccoma.obj                             /* AIC COMMUNICATION WITH INTERRUPT  */
-o mr10mcr.out                          /* OUTPUT FILE                       */
MEMORY
{
  VECS:   org = 0H       len = 0C0H     /* INTERRUPT VECTORS                 */
  SRAM:   org = 0C0H     len = 03F40H   /* USER STATIC RAM                   */
  RAM:    org = 0809800H len = 0800H    /* INTERNAL ZERO WAIT STATE RAM      */
}
SECTIONS
{
  VECTORS:  {} > VECS               /* INTERRUPT VECTORS                     */
  .text:    {} > SRAM               /* CODE                                  */
  .data:    {} > SRAM               /* SRAM-RESIDENT TABLES/CONSTANTS        */
  .bss:     {} > RAM                /* INTERNAL RAM                          */
  F_DATA    align(64): {} > RAM     /* ALIGN CRCLR BUFFER FOR FILTER INPUT   */
  I_DATA    align(16): {} > RAM     /* ALIGN CRCLR BUFFER FOR INTRPLTR INPUT */
}
