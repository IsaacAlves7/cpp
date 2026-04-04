/*
	Single Density Disk Copy Utility.
	Reads source disk a track at a time & writes track image
to destination disk a track at a time. Implements Skew 6 algorithm
for faster disk access. Allows repeats & pause for insertion of
system disk at termination.

*/
#define	numsecs	26		/* number of sectors to copy */
#define	numtrks	77		/* number of tracks to copy */
#define	version	"1.0"
#define	sect_size	128	/* bytes per sector */
#define space	0x20

char	*buffer_start,		/* ptr to start of buffer */
	*bufr_ptr;		/* ptr to current sector buffer */
unsigned	track,		/* current track # */
		sector,		/* current sector # */
		temp;		/* temporary storage */
char	xlate[] = {0,1,7,13,19,25,5,11,17,23,3,9,15,21,2,8,14,20,
			26,6,12,18,24,4,10,16,22};

main(){

   puts("\nSingle Density Disk Copy Utility   Ver",version);
   puts("\n(c) 1982   GRH Electronics  Cupertino, CA");
   puts("\n\nPlace SOURCE disk in Drive A, DESTINATION Disk in drive B &");
   puts("\n  press space key when ready. (RETURN to abort)");
   puts("\n\nCAUTION, destination disk data will be destroyed by this program!");
   buffer_start = alloc(numsecs * sect_size);	/* reserve size of track */
   if (buffer_start == 0){  puts("\n\nNot Enough Memory!"); exit(0);}
   while (getchar() == space){		/* wait for key press */
      if (bios(8,0,0) == 0){  puts("\nDrive A not Ready!"); goto abort;}
      if (bios(8,1,0) == 0){  puts("\nDrive B not ready!"); goto abort;}
      for (track = 0; track < numtrks; track++){
         bios(8,0,1);			/* re-select drive A */
         bios(9,track,0);		/* select track */
         bufr_ptr = buffer_start;
         for (sector = 1; sector <= numsecs; sector++){
            bios(10,xlate[sector],0);	/* select sector */
            bios(11,bufr_ptr,0);	/* set read ptr */
            if (bios(12,0,0)){  puts("\nRead Error!"); goto abort;}
            bufr_ptr += sect_size;
            }
         bios(8,1,1);			/* select drive B */
         bios(9,track,0);
         bufr_ptr = buffer_start;
         for (sector = 1; sector <= numsecs; sector++){
            bios(10,xlate[sector],0);
            bios(11,bufr_ptr,0);
            bios(13,0,0); /*){  puts("\nWrite Error!"); goto abort;}*/
            bufr_ptr += sect_size;
            }
         }
abort:
      puts("\n\nTo Repeat operation, Set up Source & Destination disks &");
      puts("\npress space to continue, RETURN to return to CP/M");
      }
   puts("\nPlace System disk in drive A & press any key.");
   temp = getchar();
   bios(8,0,0);			/* log on system disk */
   exit(0);
   }
