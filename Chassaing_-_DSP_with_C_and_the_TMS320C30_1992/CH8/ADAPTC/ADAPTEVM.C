/*ADAPTEVM.C-REAL-TIME ADAPTIVE USING TABLE LOOKUP*/
#include "caiccom.c"
#include "scdat"
#include "adaptcmd.h"
#define TRANSFER 8
#define NONE 0x00
#define GLOBAL 0
#define EXPANSION 0
#define PRIMARY 4
#define N 20
#define beta 1e-6
int AICSEC[7] = {0x3,0x0e1c,0x3,0x4286,0x3,0x63,0x0};

void command_process();
int ampt = 4, loop = 1;
volatile int *dma = (volatile int *) 0x808000;
volatile int *host = (volatile int *) 0x804000;
volatile int *bus = (volatile int *) 0x808060;

void main(void)
{
  double yo, E, IN_data, OUT_data, Y_high=0, Y_low=0, Y_cen=0,Y_shift;
  int I, J, K, Y1, IN_int, OUT_int;
  double W[N];
  double Delay[N];
  init_evm();
  asm("  LDI  400h,IE");
  asm("  OR  2000h,ST"); 
  Y_high=-10000;
  Y_low=10000;
  for (K=0; K < N; K++)
     { W[K]=0;
       Delay[K]=0; }
  AICSET();
  do
  {
         command_process(); 
 	 OUT_data = (((Y_high-Y_low)*(OUT_data+1)/2)+Y_low)/4;
 	 TWAIT(); 
         PBASE[0x48] = (int)OUT_data << 2;    
     	 IN_data = PBASE[0x4C] << 16 >> 18;                   
	 if (IN_data > Y_high)
	     Y_high = IN_data;
	 if (IN_data < Y_low )
	     Y_low = IN_data;
         IN_data = (2*((IN_data)-Y_low)/(Y_high-Y_low))-1;
         Y_shift=(IN_data/2)+1; 
	 Y1= Y_shift * 100;
	 if (yo > IN_data)         
	     Delay[0] = yc[Y1];
	 if (yo < IN_data)	    
 	     Delay[0] = ys[Y1];
         OUT_data = 0;
	 yo=IN_data;
	 for (I = 0; I < N; I++)
		{
		OUT_data += (W[I] * Delay[I]);
		}
	 E = IN_data - OUT_data;
	 for (J=N; J > 0; J--)
		{
		W[J] = W[J] + (2*beta*E*Delay[J]);
		if (J != 0)
			Delay[J] = Delay[J-1];
		}
 	}
  while (loop != 0);
}  

void command_process()
{
  int i;
  while(dma[TRANSFER]);
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
    while(*host & 0xFF);
    *host = NONE;
  }
}

void c_int11()
{
  while(*host & 0xFF);
  *host = NONE;
  dma[GLOBAL] = 0;
  asm("  AND  0FFFFh,IE");
}        

void init_evm()
{ 
  bus[EXPANSION] = 0x0;
  bus[PRIMARY] = 0x1000;
  *host = NONE;
  asm("  OR    800h,ST");
}     

void c_int05()
{
  for (;;);
}

void c_int99()
{
  for (;;);
}





