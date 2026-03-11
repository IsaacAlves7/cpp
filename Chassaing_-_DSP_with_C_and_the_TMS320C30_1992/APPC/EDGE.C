/*EDGE.C-EDGE ENHANCEMENT PROGRAM*/
main()
{
 volatile int *IO_INPUT = (volatile int *) 0x804000;
 volatile int *IO_OUTPUT = (volatile int *) 0x804001;
 int A[3] = {0,0,0};
 int i, result;
  for (i = 0; i < 75; ++i)
   {
    A[2] = *IO_INPUT;
    *IO_OUTPUT = 3*A[1] - A[2] - A[0];
    A[0] = A[1];
    A[1] = A[2];
   }
}

