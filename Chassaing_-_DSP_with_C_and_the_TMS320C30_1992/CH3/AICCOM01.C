/*AICCOM01.C-COMMUNICATION ROUTINES FOR BOTH SERIAL PORTS 0 AND 1 OF AIC   */
#define TWAIT while (!(PBASE[0x40+SP] & 0x2)) /*wait till XMIT buffer clear*/
#define SP0 0x00                              /*serial port 0 offset       */
#define SP1 0x10                              /*serial port 1 offset       */
extern int AICSEC[4], AICSEC1[4];             /*array defined in main prog */
volatile int *PBASE = (volatile int *) 0x808000; /*peripherals base addr   */ 

void AICSET(int SP)                        /*function to initialize AIC  */ 
{
  volatile int loop;                       /*declare local variables     */
  PBASE[0x28+SP] = 0x00000001;             /*set timer period            */ 
  PBASE[0x20+SP] = 0x000002C1;             /*set timer control register  */
  if (!SP) asm("   LDI    00000062h,IOF"); /*set IOF low to reset PRI AIC*/ 
  else asm("   LDI    00000026h,IOF");     /*set IOF low to reset AUX AIC*/  
  for (loop = 0; loop < 50; loop++);       /*keep IOF low for a while    */
  PBASE[0x42+SP] = 0x00000111;             /*set xmit port control       */
  PBASE[0x43+SP] = 0x00000111;             /*set receive port control    */
  PBASE[0x40+SP] = 0x0E970300;             /*set serial port global reg  */
  PBASE[0x48+SP] = 0x00000000;             /*clear xmit register         */
  if (!SP) asm("   OR    00000006h,IOF");  /*set IOF high enable AIC PRI */    
  else asm("   OR    00000060h,IOF");      /*set IOF high enable AIC AUX */
  for (loop = 0; loop < 4; loop++)         /*loop to configure AIC       */ 
  {
    TWAIT(SP);                             /*wait till XMIT buffer clear */
    PBASE[0x48+SP] = 0x3;                  /*enable secondary comm       */
    TWAIT(SP);                             /*wait till XMIT buffer clear */
    if (!SP)  PBASE[0x48] = AICSEC[loop];  /*secondary command for SP0   */
    else PBASE[0x58] = AICSEC1[loop];      /*secondary command for SP1   */
  }
}

void AICSET_I(int SP)                      /*configure AIC, enable TINTO */
{
  AICSET(SP);                              /*function to configure AIC   */
  asm("   LDI     00000000h,IF");          /*clear IF Register           */
  if (!SP) asm("   OR      00000010h,IE"); /*enable EXINT0 CPU interrupt */
  else asm("   OR      00000040h,IE");     /*enable EXINT1 CPU interrupt */ 
  asm("   OR      00002000h,ST");          /*global interrupt enable     */
}
 
int UPDATE_SAMPLE(int SP, int output)      /*function to update sample   */
{
  int input;                               /*declare local variable      */
  TWAIT(SP);                               /*wait till XMIT buffer clear */
  PBASE[0x48+SP] = output << 2;            /*left shift and output sample*/
  input = PBASE[0x4C+SP] << 16 >> 18;      /*input sample and sign extend*/
  return(input);                           /*return new sample           */
}  
int GET_SAMPLE(int SP)                     /*function to input sample    */
{
  int input;                               /*declare local variable      */
  input = PBASE[0x4C+SP] << 16 >> 18;      /*input and sign extend sample*/
  return(input);                           /*return new sample           */
}
void PUT_SAMPLE(int SP, int output)        /*function to output sample   */
{
  TWAIT(SP);                               /*wait till XMIT buffer clear */
  PBASE[0x48+SP] = output << 2;            /*left shift and output sample*/
} 







