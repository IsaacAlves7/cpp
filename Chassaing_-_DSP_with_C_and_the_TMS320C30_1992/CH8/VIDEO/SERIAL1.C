#define SETSP1 0x02900304       
#define GOSP1 0x0C000000

volatile int *PBASE = (volatile int *) 0x808000; 

void SERIAL1SET()
{
  PBASE[0x38] = 0x04;
  PBASE[0x30] = 0x2C1;
  PBASE[0x50] = SETSP1;
  PBASE[0x52] = 0x111;
  PBASE[0x53] = 0x111;
  PBASE[0x50] = GOSP1 | PBASE[0x50];
}

void RWAIT()
{
  int temp;
  do
  {
    temp = PBASE[0x50];
    temp = temp & 0x1;
  }
  while (temp == 0);
}  