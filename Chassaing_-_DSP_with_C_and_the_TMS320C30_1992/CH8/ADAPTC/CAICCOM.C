
/*AIC COMMUNICATION ROUTINE   */
#define SETSP 0x2970300       /*SERIAL PORT setup data*/
#define GOSP 0x0c000000
extern int AICSEC[7];
volatile int *PBASE = (volatile int *) 0x808000; 
void AICSET()
{
  int loop;
  asm("   LDI 02h,IOF");  /*AIC RESET */
  PBASE[0x28] = 0;        /*init TIMER 0 PERIOD REG(TCLK0=7.5MHZ*/ 
  PBASE[0x24] = 0;        /*start TIMER 0 counter REG @ 0  */ 
  PBASE[0x20] = 0x03C1;   /*reset TIMER 0*/
  PBASE[0x48] = 0;        /*clear serial data transmit REG */
  PBASE[0x40] = SETSP;    /*SP: 16 bits, ext clks, std mode*/  
  PBASE[0x42] = 0x131;                       
  PBASE[0x43] = 0x131;
  PBASE[0x40] = GOSP | PBASE[0x40];
  asm("   OR 06h,IOF");
  for (loop = 0; loop <= 6; loop++)
  {
    TWAIT();
    PBASE[0x48] = AICSEC[loop];
  }
}
void TWAIT()
{
  int temp;  
  do
  {
    temp = PBASE[0x40];
    temp = temp & 0x2;
  }
  while (temp == 0);
}
  
    
      

