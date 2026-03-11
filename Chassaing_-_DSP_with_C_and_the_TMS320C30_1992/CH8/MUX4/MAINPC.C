/*MAINPC.C HOST PROGRAM FOR 4-CHANNEL MULTIPLEXER*/
/*COMPILED USING QUICK C*/
#include "incfile.h"
#define SEND 128
#define N 128
void transfer();
void send_cmd();
void set_chan();
int index1 = 0;
char chan_no[1];
int channel[1];
int buffer1 = 1;
int fft_data[256];
float fft_real[256];
extern void plot(float, int,);
void main()
{
  int i;
  system("evmreset");
  system("evmload mainevm.out");
	_setvideomode(_MAXRESMODE);
	 _getvideoconfig(&myscreen);
	maxx = myscreen.numxpixels - 1;
	maxy = myscreen.numypixels - 1;
  screen_setup();
  init_evm();
 while(1)
 {
   send_cmd();
   do
   {
   set_chan();
   _moveto(5,5);
   _outtext("     Chan # ");
   _outtext(chan_no);
   transfer();
   }
   while (!kbhit());
 }
 getch();
 _setvideomode(_DEFAULTMODE);
}

void transfer()
{
  int i,value;
  while(READ_CMD);
  WRITE_CMD(128);
  while (READ_CMD != 128);
  CLR_READ_ACK;
  READ_DATA;
  for(i = 0; i < N; i++)
  {
    do
      UPDATE_STATUS0;
	while(!IS_READ_ACK);
    CLR_READ_ACK;
    fft_real[i] = (float)READ_DATA ;
  }

   plot_wave((float *)fft_real, N);
   WRITE_CMD(NONE);
}

void send_cmd()
{
  int command;
  command = get_command();
  switch(command)
  {
    case CHAN1:
	command = 0x11;
	*chan_no = '1';
	break;
    case CHAN2:
	*chan_no = '2';
	command = 0x23;
	break;
    case CHAN3:
	command = 0x45;
	*chan_no = '3';
	break;
    case CHAN4:
	command = 0x87;
	*chan_no = '4';
	break;
    default:
	 command = 0x11;
	 *chan_no = '1';
	 break;
  }
   channel[0] = command;
}
void set_chan()
{
/*  Transfer 1 word to EVM   */
    while(READ_CMD);
    WRITE_CMD(64);
    while(READ_CMD != 64);
    CLR_WRITE_ACK;
    WRITE_DATA(*channel);
    do
	 UPDATE_STATUS0;
    while (!IS_WRITE_ACK);
    CLR_WRITE_ACK;
   WRITE_CMD(NONE);
}






