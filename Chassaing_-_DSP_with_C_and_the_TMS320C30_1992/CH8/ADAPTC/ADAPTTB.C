/*ADAPTTB.C-ADAPTIVE FILTER USING ASIN,ACOS       */
#define beta 0.004     /*rate of adaptation       */
#define N  20          /*order of filter (N+1)    */
#define NS  256        /*number of samples        */
#define Max 8192       /*max value of input signal*/
#include "scdat"       /*table for asin,acos      */

main()
{
  int I, J, T, Y;
  double E, yo, IN_data, OUT_data;
  double W[N+1];
  double Delay[N+1];
  volatile int *IO_INPUT  = (volatile int*) 0x804000;
  volatile int *IO_OUTPUT = (volatile int*) 0x804001;
  yo=0;
  for (T=0; T <= N; T++)
  {
    W[T] = 0.0;
    Delay[T] = 0.0;
  }
  for (T=0; T < NS; T++)
  {
    IN_data = *IO_INPUT;     /*input signal                    */
    IN_data = IN_data/Max;   /*scale for range between 1 and -1*/
    Y=((IN_data)+1)*100;     /*step up array between 0 and 200 */
    if (yo > IN_data)        /*is signal falling or rising     */
      Delay[0] = yc[Y];      /*signal is falling, acos domain  */
    else
      Delay[0] = ys[Y];      /*signal is rising, asin domain   */
    OUT_data = 0;            /*init filter output to zero      */
    yo=IN_data;              /*store input                     */
    for (I = 0; I <= N; I++)
      OUT_data += (W[I] * Delay[I]);     /*filter output       */
    E = IN_data - OUT_data;              /*error signal        */
    for (J=N; J >= 0; J--)
    {
      W[J] = W[J] + (beta*E*Delay[J]);   /*update coefficients */
      if (J != 0)
      Delay[J] = Delay[J-1];             /*update data samples */
    }
    *IO_OUTPUT = OUT_data*Max;           /* output signal      */
  }
}





