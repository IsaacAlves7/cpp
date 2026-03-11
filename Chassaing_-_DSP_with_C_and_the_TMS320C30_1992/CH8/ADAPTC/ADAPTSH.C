/*ADAPTSH.C-ADAPTIVE FILTER WITH SHIFTED INPUT*/
#include <math.h>
#define beta 0.004                       
#define N  20                           
#define NS  256                          
#define Pi  3.1415926                     
#define shift 90
#define max_input 8192

main()
{
  int I,T;
  float IN_data;
  double Y,Y1,X, E, D, yo; 
  double W[N+1];
  double Delay[N+1];
  volatile int *IO_INPUT = (volatile int*) 0x804000;
  volatile int *IO_OUTPUT = (volatile int*) 0x804001;
  yo=0;
  for (T=0; T <= N; T++)
  {
    W[T] = 0.0;
    Delay[T] = 0.0;
  }
  for (T=0; T < NS; T++)        
  {  
    IN_data = *IO_INPUT;
    Y=IN_data/max_input;
    if (yo>=Y)
      X=acos(Y);
    else
      X=asin(Y)+(3*Pi/2);
    X=X-(shift);
    Delay[0]=cos(X);
    D=IN_data;
    Y1 = 0;                           
    yo=Y; 
    for (I = 0; I <= N; I++)
      Y1 += (W[I] * Delay[I]);              
    E = D - Y1;                       
    for (I = N; I >= 0; I--)         
    {
      W[I] = W[I] + (beta*E*Delay[I]);  
      if (I != 0)                     
      Delay[I] = Delay[I-1];              
    } 
    *IO_OUTPUT = Y1;
  }   
}





