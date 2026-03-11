/*IIR6BPC.C-SIXTH-ORDER IIR BANDPASS FILTER,Fc=1250 Hz  */
#define STAGES 3            /*number of 2nd-order stages*/
float A[STAGES][3]=        {/*numerator coefficients    */
{5.3324E-02, 0.0000E+00, -5.3324E-02},
{5.3324E-02, 0.0000E+00, -5.3324E-02},
{5.3324E-02, 0.0000E+00, -5.3324E-02} };
float B[STAGES][2]=        {/*denominator coefficients  */
{-1.4435E+00, 9.4879E-01},
{-1.3427E+00, 8.9514E-01},
{-1.3082E+00, 9.4377E-01} };
float DLY[STAGES][2]= {0};  /*delay samples             */
void IIR(int *IO_in, int *IO_out, int n, int len)
{
 int i, loop = 0;
 float dly, yn, input;
 while (loop < len)
 {
  ++loop;
  input = *IO_in;   
  for (i = 0; i < n; i++)
  {
  dly = input-B[i][0]*DLY[i][0]-B[i][1]*DLY[i][1];
  yn  = A[i][2]*DLY[i][1]+A[i][1]*DLY[i][0]+A[i][0]*dly;
  DLY[i][1] = DLY[i][0];
  DLY[i][0] = dly;
  input = yn;
  }
  *IO_out = yn;
 }
}
main()
{
 #define length 345
 volatile int *IO_IN  = (volatile int *) 0x804000;
 volatile int *IO_OUT = (volatile int *) 0x804001;
 IIR((int *)IO_IN, (int *)IO_OUT, STAGES, length);
} 
   


