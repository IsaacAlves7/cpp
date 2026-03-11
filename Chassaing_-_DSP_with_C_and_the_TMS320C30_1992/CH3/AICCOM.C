/*AICCOM.C-COMUNICATION ROUTINES FOR SERIAL PORT 0 OF AIC */
#define TWAIT while (!(PBASE[0x40] & 0x2)) /*wait till XMIT buffer clear*/
extern int AICSEC[4];                      /*array defined in main prog */
volatile int *PBASE = (volatile int *) 0x808000; /*peripherals base addr*/ 

void AICSET()                        /*function to initialize AIC */ 
{
  volatile int loop;                 /*declare local variables    */
  PBASE[0x28] = 0x00000001;          /*set timer period           */ 
  PBASE[0x20] = 0x000002C1;          /*set timer control register */
  asm("    LDI      00000002h,IOF"); /*set IOF low to reset AIC   */ 
  for (loop = 0; loop < 50; loop++); /*keep IOF low for a while   */
  PBASE[0x42] = 0x00000111;          /*set xmit port control      */         
  PBASE[0x43] = 0x00000111;          /*set receive port control   */
  PBASE[0x40] = 0x0E970300;          /*set serial port global reg */
  PBASE[0x48] = 0x00000000;          /*clear xmit register        */
  asm("    OR       00000006h,IOF"); /*set IOF high to enable AIC */    
  for (loop = 0; loop < 4; loop++)   /*loop to configure AIC      */ 
  {
    TWAIT;                           /*wait till XMIT buffer clear*/
    PBASE[0x48] = 0x3;               /*enable secondary comm      */
    TWAIT;                           /*wait till XMIT buffer clear*/
    PBASE[0x48] = AICSEC[loop];      /*secondary command for SP0  */
  }
}

void AICSET_I()                      /*configure AIC, enable TINTO*/
{
  AICSET();                          /*function to configure AIC  */
  asm("    LDI      00000000h,IF");  /*clear IF Register          */
  asm("    OR       00000010h,IE");  /*enable EXINT0 CPU interrupt*/
  asm("    OR       00002000h,ST");  /*global interrupt enable    */
}
 
int UPDATE_SAMPLE(int output)       /*function to update sample   */
{
  int input;                        /*declare local variables     */
  TWAIT;                            /*wait till XMIT buffer clear */
  PBASE[0x48] = output << 2;        /*left shift and output sample*/
  input = PBASE[0x4C] << 16 >> 18;  /*input sample and sign extend*/
  return(input);                    /*return new sample           */
}  

 



