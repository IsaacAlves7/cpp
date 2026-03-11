/*EVMCTRL.C-TO PLACE TMS320C30 IN OR OUT OF RESET*/
#include "stdio.h"
main(argc,argv)
int argc;
char *argv[];
{
  if (argc!=2)
  {
    printf("\n\n\tUsage: evmctrl [go] or [stop]\n\n");
    exit(1);
  }
  if (strncmp(argv[1],"go")==0)
  {
    printf("\n\n\tStarting C30\n\n");
    outport(0x0240 + 0x0A, 0x800);
    exit(0);
  }
  if (strncmp(argv[1],"stop")==0)
  {
    printf("\n\n\tStopping C30\n\n");
    outport(0x0240 + 0x0A, 0x808);
    exit(0);
  }
  printf("\n\n\tUsage: evmctrl [go] or [stop]\n\n");
  exit(1);
}




