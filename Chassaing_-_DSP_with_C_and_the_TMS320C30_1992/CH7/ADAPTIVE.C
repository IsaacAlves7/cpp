/* Program for adaptive filtering using the least means square adaptation */

#include <math.h>
#include <dos.h>
#include <stdio.h>
#include <conio.h>
#include <graphics.h>

void initial_graphics (void);
float graph_outline (void);
int format ( double position);
void LJ_Graphic();
void draw_line (double D, double Y, long T, double D_1, double Y_1);
char key_status (void);

#define BACKGROUND BLACK
#define COLORTEXT LIGHTGRAY
#define TRACE1 RED
#define TRACE2 BLUE


#define N  30
#define NS  40
#define FS  8000
#define Pi  3.1415926
#define YSCALE  2
#define XSCALE  NS/FS
#define DESIRED cos(2*Pi*T*1000/FS)
#define NOISE sin(2*Pi*T*1000/FS)

main()
{
  char ch;
  float beta;
  long loop, I, T;
  double D, Y, E, D_1, Y_1; 
  double W[N+1] = {0.0};
  double X[N+1] = {0.0};

  initial_graphics();                   //initialize graphics
  do                                    //repeat graph for new beta value
  {
    beta = graph_outline();             //draw graph and read in value of beta
    for (loop = 0; loop <= N+1; loop++) //initialize arrays to zero for new beta
    {  
      W[loop] = 0.0;                    //coefficient array
      X[loop] = 0.0;                    //array for noise samples
    }
    for (T = 0; T <= NS; T++)           //start adaptive alogrithm
    {    
      X[0] = NOISE;                     //new noise sample
      D = DESIRED;                      //desired signal
      Y = 0;                            //set output of filter to zero
      for (I = 0; I <= N; I++)
	Y = Y + (W[I] * X[I]);          //calculate filter output
      E = D - Y;                        //calculate error signal
      for (I = N; I >= 0; I--)
      {
	W[I] = W[I] + (2*beta*E*X[I]);  //update filter coefficients
	if (I != 0)                     
	  X[I] = X[I-1];                //move data sample 
      }      
      draw_line (D,Y,T,D_1,Y_1);        //plot output to display
      D_1 = D;                          //D_1 is previous desired sample
      Y_1 = Y;                          //Y_1 is previous output sample
    }   
    ch = key_status();                  //get instruction from user
    setviewport (0,0,getmaxx(),getmaxy(),0);   
    clearviewport();                    //clear screen for plot with new beta
  }
  while (ch == 0);    
  closegraph();                         //return to text mode
  return 0;
}

char key_status (void)
{
  char ch;

  while (!kbhit());
  ch = getch();
  if (ch == 0)
  {
    ch = getch();
    if (ch == 59)
      LJ_Graphic();
    if (ch == 60)
      return(0);
  }
  return(1);
}


void draw_line(double D, double Y, long T, double D_1, double Y_1)
{
  if (T != 0)
  {
    setcolor(TRACE2);
    setlinestyle(DASHED_LINE,1,1);
    line((T-1)*(400/NS), (150 - (D_1 * 150/YSCALE)), (T * 400/NS),
	 (150 - (D * (150/YSCALE))));
    setcolor(TRACE1);
    setlinestyle(SOLID_LINE,1,1);
    line((T-1)*(400/NS), (150 - (Y_1 * 150/YSCALE)), (T * 400/NS),
	 (150 - (Y * (150/YSCALE))));
    setcolor(COLORTEXT);
  }
}

void initial_graphics (void)
{

  int gdriver = DETECT, gmode, errorcode;    // request auto detection
  initgraph(&gdriver, &gmode, "");           // initialize graphics mode
  errorcode = graphresult();            // read result of initialization 
  if (errorcode != grOk)                     // an error occurred
    {
       printf("Graphics error: %s\n", grapherrormsg(errorcode));
       printf("Press any key to halt:");
       getch();
       exit(1);                                // return with error code
    }
}


float graph_outline (void)
  {
    int left, top, right, bottom, loop;
    char str[25];
    float beta;

    printf ("Enter a real number for beta (Example 0.05).  ");
    scanf ("%f", &beta);
    setbkcolor(BACKGROUND);
    cleardevice();
    setcolor(COLORTEXT);
    left = getmaxx() / 2 - 200;
    top = getmaxy() / 2 - 200;
    right = getmaxx() / 2 + 200;
    bottom = getmaxy() / 2 + 100;
    settextjustify(CENTER_TEXT,CENTER_TEXT);
    outtextxy (getmaxx()/2,top-20, "ADAPTIVE FILTER");
    rectangle (left, top, right, bottom);
    setviewport (left+1,top+1,right-1,bottom-1,1);
    setlinestyle(DASHED_LINE, 1,1);
    line(0,150,400,150);
    setlinestyle(SOLID_LINE,1,1);
    for (loop = 1; loop <= 7; loop++)
      {
	line(0,(loop*75)/2,5,(loop*75)/2);        //draws horz hash marks
	line(loop*80,300,loop*80,292);    //draws vert hash marks
      }
    setviewport (left-50,top-10,right,bottom+10,1);
    for (loop = 4; loop >= -4; loop--)
      {
	gcvt((-loop/(4.0)*YSCALE),2,str);
	outtextxy(25,(loop*75)/2+160, str);
      }
    setviewport (left-20,top,right+50,bottom+140,0);
    for (loop = 0; loop <= 5; loop++)
      {
	gcvt(loop/(5.0)*XSCALE, 10, str);
	outtextxy((loop*80) + 20, 315, str);
      }
    outtextxy(220,340, "TIME");
    settextstyle(DEFAULT_FONT, VERT_DIR, 1);
    outtextxy (-50,150, "AMPLITUDE");
    settextstyle(DEFAULT_FONT, HORIZ_DIR, 1);
    settextjustify(LEFT_TEXT,CENTER_TEXT);
    setcolor(TRACE1);
    line(50,365,75,365);
    setlinestyle(DASHED_LINE, 1,1);
    setcolor(COLORTEXT);
    outtextxy(80,365, "OUTPUT SIGNAL");
    setcolor(TRACE2);
    line(250,365,275,365);
    setcolor(COLORTEXT);
    outtextxy(280,365, "DESIRED SIGNAL");
    settextjustify(RIGHT_TEXT,CENTER_TEXT);
    outtextxy(150,385, "SAMPLING FREQ = ");
    outtextxy(350,385, "SAMPLES = ");
    outtextxy(150,400, "ORDER OF FILTER = ");
    outtextxy(350,400, "BETA = ");
    settextjustify(CENTER_TEXT,CENTER_TEXT);
    outtextxy(75,425, "F1 for Printout");
    outtextxy(225,425, "F2 for New BETA");
    outtextxy(375,425, "ENTER to Quit");
    settextjustify(LEFT_TEXT,CENTER_TEXT);
    gcvt(FS*1.0, 10, str);
    outtextxy (150, 385, str);
    gcvt(NS*1.0, 10, str);
    outtextxy (350, 385, str);
    gcvt(N+1.0, 10, str);
    outtextxy (150, 400, str);
    gcvt(beta, 6, str);
    outtextxy (350, 400, str);
    settextjustify(CENTER_TEXT,CENTER_TEXT);    
    setviewport (left+1,top+1,right-1,bottom-1,1);
    return(beta); 
}


void LJ_Graphic()
{
  int xaspect, yaspect, maxX, maxY, line, xword, pixel, xwidth, ywidth;
  double xpos, ypos, prnstep, ratio;
  char chr;
				   
  maxX = getmaxx();                         //get number of horizontal pixel
  maxY = getmaxy();                         //get number of vertical pixel
  getaspectratio(&xaspect, &yaspect);       //get the screen aspect ratio
  ratio = (double) xaspect/ (double) yaspect;      
  setviewport(0,0,maxX,maxY,0);             //set viewport for full screen

  xpos = 690;                               //initial position of prn cursor
  ypos = 500;                               //initial position of prn cursor
  prnstep = 7.2/ratio;                      //match printer aspr fo screen aspr
  fprintf (stdprn, "\x1B&E\x1B&11H\x1B&1O\x1B*p0X\x1B*p0Y\x1B*t100R");
  for (line = 0; line <= maxY; line++)
  {
    ywidth = 6;
    if (ypos < 1000.0) ywidth--;
    if (ypos < 100.0) ywidth--;
    if (ypos < 10.0) ywidth--;
    fprintf (stdprn, "\x1B&a%-*.1fh%-*.1fV", 5, xpos, ywidth, ypos);
    ypos += prnstep;
    fprintf (stdprn, "\x1B*r1A\x1B*b%dW", maxX/8);
    for (xword = 0; xword < maxX/8; xword++)
    {
      chr = 0;
      for (pixel = 0; pixel < 8; pixel++)  
      {                     //reads series of 8 pixels to create graphics char
	chr <<= 1;
	if (getpixel (xword*8+pixel, line)) chr++;
      }
      fprintf (stdprn, "%c", chr);             //sends graphic char to printer
    }
    fprintf (stdprn, "\x1B*rB");               //ends graphics line
  }
  fprintf (stdprn, "\x0C\x1B&10O\x1B&11H\x1B&E");
}
		    
