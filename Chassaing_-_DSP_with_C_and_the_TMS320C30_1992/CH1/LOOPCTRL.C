/*LOOPCTRL.C-LOOP PROGRAM WITH CONTROL, WITH GRAPHICS*/
#include "aiccom.c"
#include "loopcmd.h"
int AICSEC[4] = {0x1428,0x1,0x4A96,0x67};
int ampt = 4, loop = 1;
volatile int *host = (volatile int *) 0x804000;

void command_process()
{
  int i;
  i = *host & 0xFF;
  if (i && i < 128)    /* valid keystroke */
  {
    *host = i;
    switch(i)
    {
      case HIGHER: if (ampt < 8) {++ampt; break;}
                   else break;
      case LOWER:  if (ampt > 1) {--ampt; break;}
                   else break;
    }
    *host = 0x00;
  }
}

void main(void)
{
  int data_IN, data_OUT;
  AICSET();
  do
  {
    command_process(); 
    data_IN = UPDATE_SAMPLE(data_OUT);   
    data_OUT = data_IN / ampt;
  }
  while (loop != 0);
}  
  

