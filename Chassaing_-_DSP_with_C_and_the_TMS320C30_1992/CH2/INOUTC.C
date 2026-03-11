/*INOUTC.C-DEMONSTRATES INPUT AND OUTPUT        */
main()
{
 volatile int *IO_INPUT=(volatile int *) 0x804000;
 volatile int *IO_OUTPUT=(volatile int *)0x804001;
 int count,value;
 for (count=0; count <=4; ++count)
  {             
   value = *IO_INPUT;           /*input -> value*/
   *IO_OUTPUT = 2 * value;      /*output=2*input*/
  }
}  

