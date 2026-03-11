/*BP45MC.C-FIR WITH C AND C-CALLED ASSEMBLY FUNCTION BP45F.ASM*/
#define N 45                      /*length of impulse response*/
float DLY[2*N];                   /*delay samples             */
const float H[N] = {/* filter coefficients*/
-1.839E-03,-2.657E-03,-4.312E-10, 3.154E-03, 2.595E-03,-4.159E-03, 
-1.540E-02,-2.507E-02,-2.547E-02,-1.179E-02, 1.392E-02, 4.206E-02, 
 5.888E-02, 5.307E-02, 2.225E-02,-2.410E-02,-6.754E-02,-8.831E-02, 
-7.475E-02,-2.956E-02, 3.030E-02, 8.050E-02, 1.000E-01, 8.050E-02, 
 3.030E-02,-2.956E-02,-7.475E-02,-8.831E-02,-6.754E-02,-2.410E-02, 
 2.225E-02, 5.307E-02, 5.888E-02, 4.206E-02, 1.392E-02,-1.179E-02, 
-2.547E-02,-2.507E-02,-1.540E-02,-4.159E-03, 2.595E-03, 3.154E-03, 
-4.312E-10,-2.657E-03,-1.839E-03};
extern void filt(float *, float *, int *, int *, int);

main ()
{
 int loop;
 volatile int *IO_INPUT = (volatile int *) 0x804000; /*in port addr */
 volatile int *IO_OUTPUT= (volatile int *) 0x804001; /*out port addr*/ 
 for (loop = 0; loop < 2*N; loop++) DLY[loop] = 0.0; /*init samples */
 filt((float *)H, (float *)DLY, (int *)IO_INPUT, (int *)IO_OUTPUT, N);
}











