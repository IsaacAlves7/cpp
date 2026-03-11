/*VOLATILE.C - USE OF VOLATILE VARIABLE*/
#include <dos.h>
#include <stdio.h>
void interrupt (*old_int_routine)(void);/*declare interrupt function pointer*/
volatile int clock = 0;                 /*declare a volatile variable*/

void interrupt timer()
{
  clock++;
  if (clock % 18 == 0)               /*approximately 18.2 clock ticks/second*/
    printf("The following seconds have elapsed %i\n", clock/18);
  old_int_routine();                 /*call old interrupt routine*/
}

main()
{
  old_int_routine = getvect(0x8);    /*save original interrupt vector*/
  setvect(0x8, timer);               /*install new interrupt handler*/
  while (clock < 200);               /*idle for 200 clock ticks */
  setvect(0x8, old_int_routine);     /*install orginal interrupt vector*/
}