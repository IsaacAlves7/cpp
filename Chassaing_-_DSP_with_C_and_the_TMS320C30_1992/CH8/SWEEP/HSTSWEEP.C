#include <stdlib.h>  /* INCLUDE HEADERS FOR STANDARD LIBRARY FILES USED */
#include <stdio.h>
#include <math.h>
#include <conio.h>
#include <ctype.h>
#include <bios.h>    /* INCLUDE TURBO C SPECIFIC HEADER FILES           */
#include <graphics.h>
#include <dos.h>
#include <time.h>
#include "pc_1.h"             /* EVM SUPPORT STRUCTURES AND MACROS      */
#include "pc_2.h"             /* EVM SUPPORT GLOBAL VARIABLES           */

#define DRIVER 9
#define MODE 2

#define FFTCMD 1

void init_evm(void)
{
    time_t start,finish;
    double delta;
    get_iobase();             /* GET I/O BASE ADDRESS OF EVM            */
    time(&start); delta = 0.0;
    UPDATE_STATUS0;
    while(!IS_READ_ACK)
    {
	time(&finish);
	delta = difftime(finish,start);
	if(delta > 10.0) exit(-1);
	UPDATE_STATUS0;
    }
    CLR_READ_ACK;
    WRITE_DATA(NONE);         /* CLEAR HOST REGISTER VALUE TO EVM       */
}

void get_iobase(void)
{
    char *d_options   = getenv("D_OPTIONS");  /* GET D_OPTION ENV VARIABLE  */
    int  found_base   = OFF;                  /* IOBASE (-p) OPTION FOUND   */
    if (!d_options) return;                   /* IF NO ENVVAR RETURN DEFAULT*/
    while (*++d_options)                      /* SEARCH ENTIRE STRING       */
    {
        if (d_options[-1] == '-' && toupper(d_options[0]) == 'P')
        {
            found_base = ON;
            break;
        }
    }
    if (!found_base) return;           /* IF -p NOT FOUND RETURN DEFAULT*/
    sscanf(++d_options, "%x", &iobase);/* OTHERWISE READ IN IOBASE      */
    switch(iobase)
    {
        case 0x240:
        case 0x280:
        case 0x320:
        case 0x340:
             break;
        default:
             puts("Invalid I/O base location found in");
             puts("D_OPTIONS environment variable.");
             puts("Exiting...");
             exit(EXIT_FAILURE);
    }
}

void init_graph()
    {
    int drivr,mod;
    int i;
    int freq_val;

    char *freq_str;

    drivr=DRIVER;
    mod=MODE;

    initgraph(&drivr,&mod,"c:\\turbo\\c");
    setlinestyle(0,0,3);
    setcolor(2);
    rectangle(97,67,503,253);
    setcolor(14);
    setlinestyle(0,0,1);
    for(i=70;i<=250;i+=30)
	line(86,i,95,i);
    for(i=100;i<=500;i+=20)
	line(i,255,i,264);
    settextjustify(0,1);
    setcolor(12);
    settextstyle(1,0,1);
    moveto(220,40);
    outtext("Sweep Generator");
    setcolor(14);
    settextstyle(2,1,4);
    moveto(30,155);
    outtext("Magnitude");
    settextstyle(2,0,4);
    moveto(48,70);
    outtext("+20 dB");
    moveto(48,130);
    outtext("  0 dB");
    moveto(48,190);
    outtext("-20 dB");
    moveto(48,250);
    outtext("-40 dB");
    freq_val=0;
    for(i=1;i<20;i+=2)
	{
	moveto(90+i*20,270);
	freq_val+=200;
	freq_str=itoa(freq_val,freq_str,10);
	outtext(freq_str);
	}
    moveto(250,290);
    outtext("Frequency (Hz)");
    setcolor(12);
    settextstyle(2,0,5);
    moveto(220,320);
    outtext("Press any key to exit");
    setviewport(100,70,500,250,1);
    setfillstyle(0,0);
    setcolor(11);
    return;
    }

void get_dBs(int *volt_vals,int *dB_vals,int np)
    {
    int i;
    double voltage;

    for (i=0;i<np;i++)
	{
	voltage=(double)volt_vals[i];
	dB_vals[i]=(int)(20 * log10((voltage/3100)));
	}
    return;
    }

void graph_sweep(int *values,int np)
    {
    int i;
    int x,y;

    floodfill(200,100,2);
    x=0;
    y=60-values[0]*3;
    if (y>179)
	y=179;
    moveto(x,y);
    for (i=1;i<np;i++)
	{
	x+=21;
	y=60-values[i]*3;
	if (y>179)
	    y=179;
	lineto(x,y);
	}
    return;
    }


main()
    {
    int *sweep_values;
    int *dB_values;
    int num_points;
    int key_pressed;
    int i,j;
    int ch;

    init_graph();
    init_evm();
    rewind(stdin);
    num_points=20;
    while(1)
	{
	for(i=0;i<num_points;i++)
	    {
	    WRITE_DATA(0x00);
	    while(!READ_DATA);
	    sweep_values[i]=READ_DATA;
	    WRITE_DATA(sweep_values[i]);
	    while(READ_DATA);
	    }
	get_dBs(sweep_values,dB_values,num_points);
	graph_sweep(dB_values,num_points);
	if(kbhit())
	    {
	    restorecrtmode();
	    clrscr();
	    exit(1);
	    }
	}
    }