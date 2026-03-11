/*LOOPALL.C-LOOP AND COMUNICATIONS ROUTINES FOR SERIAL 0 PORT AIC       */
#define TWAIT while (!(PBASE[0x40] & 0x2)) /*wait till XMIT buffer clear*/
int AICSEC[4]= {0x1428,0x1,0x4A96,0x67};   /*config data for SP0 AIC    */
volatile int *PBASE = (volatile int *) 0x808000; /*peripherals base addr*/ 

void AICSET()                        /*function to initialize AIC    */ 
{
 volatile int loop;                  /*declare local variables       */
 PBASE[0x28] = 0x00000001;           /*set timer period              */ 
 PBASE[0x20] = 0x000002C1;           /*set timer control register    */
 asm("     LDI  00000002h,IOF");     /*set IOF low to reset AIC      */ 
 for (loop = 0; loop < 50; loop++);  /*keep IOF low for a while      */
 PBASE[0x42] = 0x00000111;           /*set xmit port control         */        
 PBASE[0x43] = 0x00000111;           /*set receive port control      */
 PBASE[0x40] = 0x0E970300;           /*set serial port global reg    */
 PBASE[0x48] = 0x00000000;           /*clear xmit register           */
 asm("     OR   00000006h,IOF");     /*set IOF high to enable AIC    */    
 for (loop = 0; loop < 4; loop++)    /*loop to configure AIC         */ 
  {                                                                  
   TWAIT;                            /*wait till XMIT buffer clear   */
   PBASE[0x48] = 0x3;                /*enable secondary comm         */
   TWAIT;                            /*wait till XMIT buffer clear   */ 
   PBASE[0x48] = AICSEC[loop];       /*secondary command for SP0     */
  }
}
                                                                    
int UPDATE_SAMPLE(int output)        /*function to update sample     */
{
  int input;                         /*declare local variables       */
  TWAIT;                             /*wait till XMIT buffer clear   */ 
  PBASE[0x48] = output << 2;         /*left shift and output sample  */
  input = PBASE[0x4C] << 16 >> 18;   /*input sample and sign extend  */
  return(input);                     /*return new sample             */
}  

main()
{
 int data_in, data_out;              /*initialize variables          */
 AICSET();                           /*call function to congig AIC   */
 while (1)                           /*create endless loop           */
 {
  data_in = UPDATE_SAMPLE(data_out); /*call function to update sample*/
  data_out = data_in;                /*loop input to output          */
 }
}  



