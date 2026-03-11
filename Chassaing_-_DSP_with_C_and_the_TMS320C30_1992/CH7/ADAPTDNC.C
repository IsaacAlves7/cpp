/*ADAPTDNC.C-50-COEFFICIENT ADAPTIVE FILTER IN C       */
#define beta 2.5E-10       /*rate of convergence       */   
#define N 50               /*# of coefficients         */
#define NS 128             /*# of output sample points */

main()
{
  int I,T;
  float Y, E, D;
  float W[N+1];
  float Delay[N+1];
  volatile int *IO_INPUT = (volatile int*) 0x804000;
  volatile int *IO_INPUT1= (volatile int*) 0x804001;
  volatile int *IO_OUTPUT= (volatile int*) 0x804002; 
  for (T=0; T < N; T++)
  {
    W[T] = 0.0;
    Delay[T] = 0.0;
  } 
  for (T=0; T < NS; T++)   /*NS is # of output samples */
  {  
    Delay[0] = *IO_INPUT1; /*noise signal n            */
    D = *IO_INPUT;         /*desired signal+noise d+n  */
    Y = 0;                           
    for (I = 0; I < N; I++)
      Y += (W[I] * Delay[I]);          /*filter output */             
    E = D - Y;                         /*error signal  */               
    for (I = N; I > 0; I--)         
    {
      W[I] = W[I] + (beta*E*Delay[I]); /*update coeffs */ 
      if (I != 0)                     
      Delay[I] = Delay[I-1];           /*update samples*/   
    } 
    *IO_OUTPUT = E;
  }   
}









