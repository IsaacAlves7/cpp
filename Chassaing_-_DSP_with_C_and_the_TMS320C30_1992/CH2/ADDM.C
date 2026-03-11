/*ADDM.C-PROGRAM IN C CALLING A FUNCTION IN ASSEMBLY*/
extern int addmfunc();  /*external assembly function*/
int temp = 10;          /*global C variable         */ 
main()
{
 volatile int *IO_OUTPUT=(volatile int *) 0x804001; /*output port addr*/
 int count;
 for (count = 0; count < 5; ++count)
  {
   *IO_OUTPUT=addmfunc(count); /*calls assembly function*/
  }
} 

