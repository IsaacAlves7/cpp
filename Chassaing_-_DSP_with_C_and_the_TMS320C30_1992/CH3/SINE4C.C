/*SINE4C.C SINE GENERATION WITH 4 POINTS.USE OF INTERRUPT               */
#define period 375                           /*timer 0 period reg value */
#define VEC_ADDR (volatile int *) 0x00       /*addr of interrupt vectors*/
#define TIMER_ADDR (volatile int *) 0x808020 /*addr of timer 0          */
int sin_table[4] = {0, 16384, 0, -16384};    /*values in sine table     */
int loop = 0;                                /*declare global variable  */
volatile int *IO_OUTPUT = (volatile int *) 0x804001; /*output port addr */

void c_int09()                   /*TINT0 Interrupt handler     */
{
 *IO_OUTPUT = sin_table[loop];  /*output value from sine table */             
 if (loop < 3) ++loop;          /*increment to next value      */
 else loop = 0;                 /*go to beginning of sine table*/
}

main()
{
  volatile int *INTVEC = VEC_ADDR;    /*pointer to interrupt vector*/
  volatile int *TIMER = TIMER_ADDR;   /*pointer to timer 0 address */
  INTVEC[9] = (volatile int) c_int09; /*Install interrupt 9 handler*/
  TIMER[8] = period;                  /*Set period register        */
  TIMER[0] = 0x2C1;                   /*Set Global Cntrl Register  */
  asm ("        LDI    100h,IE");     /*Enable TINT0 interrupt     */
  asm ("        OR     2000h,ST");    /*Enable GIE Bit             */
  for (;;);                           /*Wait for interrupt         */
}
                                        










