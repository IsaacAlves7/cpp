/* xref filename options
copyright 1991 Louis Baker. All rights reserved.

options:
	p ignore line after //
*/
#include <stdlib.h>
#include <stdio.h>
#include <dir.h>
#include <string.h>
#include <ctype.h>
#include <io.h>
#define bufflen 10
#define manykt 10
/*  hex: D= cr a= line feed C=form feed*/
#define CTLL 0x0c
#define CR 0x0d
#define LF 0x0a
#define iswhite(x) (((x)==' ')||((x)=='\t')||((x)== CR )||((x)=='\r'||((x)== LF )))
#define isd(x) (((x)==';')||((x)=='(')||((x)== '{')||((x)=='}'||((x)== ')'||iswhite(x))))
/*#define kwct 55*/
int kwct;
#define maxkw 200
char keywd[maxkw][bufflen+1];

void help(void)
{printf(" cross reference c code enter filename(wildcards ok)\n");
printf(" copyright 1991 Louis Baker. All rights reserved.\n");
}

#define MAXLEN 50

struct taben;

struct ref
	{
	struct taben *caller;/* caller or called, depending upon list*/
	struct ref *next;
	};

struct taben
	{char name[MAXLEN];
	struct taben *next;
	struct ref *callers,*calleds;
	}
	*root,*caller,*dummy;
// callers is list of pgms which are called by this one
// calleds is list of pgms which have called this one


char cursym[MAXLEN];

struct taben *entersym( char name[],int called,int len)
{
struct taben *el,*pres,*prev;struct ref *newen, *presr,*prevr;
int i,j,freenew=1;/* free newen unless told otherwise later*/
if(len<=0)
	{
/*	printf(" zero length entry\n");*/
	return NULL;
	}
/* reject keywords such as if(...*/
for(i=0;i< kwct ;i++)
	{
	j=strcmp(keywd[i],name);
/*if(!j)printf(" debug: rejecting keyword %s\n",name);*/
	if(!j)return NULL;
	}

/* reject pure number- the .-+ will be lost already
lets thru 6e8 of 1.6e8, but nobody's perfect. */
for(i=0;i<len;i++)
	{if(name[i]<'0'||name[i]>'9')goto create;
	}
return NULL;

/*create entry*/
create:
el= (struct taben *) malloc(sizeof(struct taben));
if(el==NULL) {fprintf(stderr," bad malloc\n");exit(0);}
el->next=NULL;
el->calleds=NULL;
el->callers=NULL;
for(i=0;i<len;i++) el->name[i]=name[i];el->name[len]='\000';
newen=NULL;
if(called)
	{
	/* prepare to add caller to list of callers*/
	if(caller==NULL  ){printf(" call %s out of {}\n",name);
	newen=NULL;
	goto table;
	}
	/* make entry */
	newen=(struct ref *) malloc(sizeof(struct ref));
	newen->next=NULL;
	newen->caller= caller;
	el->callers=newen;
	}
else el->callers=NULL;

table:

/*printf(" %s called=%d\n",name,called);
if(called && caller!=NULL)printf(" caller=%s\n",caller->name);
*/
if(root==NULL)
	{/* first entry (caller)*/
	if(called){printf(" first entry called,not caller\n");exit(0);}
/*printf(" root is %s\n",el->name);*/
	root=el;
	return el;
	}
/* else, search*/
/* insert before,or after if next==NULL.
if ==, add to list of callers if called. otherwise error
then free(el) and return existing
*/
pres=root;prev=NULL;
/* main loop to enter el into list of functions, and
to add newen to list of callers as required if function already tabulated
and if called rather than caller*/
for(;;)
	{
	if(pres==NULL)
		{/* add to tail, i.e. after prev*/
		if(prev==NULL){printf(" prev==NULL");exit(0);}
		prev->next=el;
/*printf(" adding as tail after %s\n",prev->name);*/
		break;
		}
	i=strcmp(name, pres->name);
	/* not strncmp(,,len) as sin part of sine but different !*/
	if(i<0)
		{/*insert before*/
/*printf(" adding before %s\n",pres->name);*/
		if(prev!=NULL)
			{el->next=prev->next;prev->next=el;}
		else /* new root*/
			{el->next=root;root=el;}
		break;
		}
	else if(!i)
		{/* same name- add entry to caller list if new*/
/*printf(" already in main table: %s\n",pres->name);*/
		if(called)
		{
		if(pres->callers==NULL)
			{
/*printf(" %s in table no previous callers. new caller=%s\n"
	,name,newen->caller->name);
*/			pres->callers=newen;
			freenew=0;
			break;
			/* do not free newen but free el*/
			}

		prevr=NULL;presr=pres->callers;
		for(;;)
			{
			if(presr==NULL)
				{/* add to tail, i.e. after prev*/
				prevr->next=newen;
/*printf(" entering at tail of called list of callers after %s\n"
	,prevr->caller->name);
*/				freenew=0;
				break;
				}
			i=strcmp(caller->name, presr->caller->name);
			/* not strncmp(,,len) as sin part of sine but different !*/
			if(i<0)
				{/*insert before*/
/*printf(" entering callers list before %s as %s\n"
	,presr->caller->name,caller->name);
*/				newen->next=prevr->next;prevr->next=newen;
				freenew=0;
				break;
				}
			else if(!i)
				{/* called more than once by same routine*/
/*printf(" called more than once by same routine. free newen\n");*/
				break;
				}
			prevr=presr;
			presr=presr->next;
			}
		}/* end if(called) entry into its caller list*/
/*printf(" freeing newly generated taben for %s\n",name);*/
		free(el);
		el=pres;/* el points to previous entry for function*/
		break;
		}
	/* else continue to look*/
	prev=pres;
	pres=pres->next;
	}
/*printf(" out of main table loop. root=%s\n",root->name);*/

if(!called)return el;/* if caller, just place in table.
do not place in caller's list of calleds*/
if(caller==NULL)
	{
/*	printf(" called, but caller==NULL\n");*/
	goto end;
	}
/*printf(" putting in caller's list of called\n");*/
/* called- put in caller list*/
		pres=caller;
		/* make entry */
		newen=(struct ref *) malloc(sizeof(struct ref));
		newen->next=NULL;
		newen->caller= el;
		if(pres->calleds==NULL)
			{pres->calleds=newen;freenew=0;
/*printf(" entering as first of caller list\n");*/
			goto end;
			}
		prevr=NULL;presr=pres->calleds;
		for(;;)
			{
			if(presr==NULL)
				{/* add to tail, i.e. after prev*/
/*printf(" entering at tail of caller list of called after %s\n"
	,prevr->caller->name);
*/				prevr->next=newen;
				freenew=0;
				break;
				}
			i=strcmp(name, presr->caller->name);
			/* not strncmp(,,len) as sin part of sine but different !*/
			if(i<0)
				{/*insert before*/
/*printf(" entering in caller list before %s\n",presr->caller->name);*/
				if(prevr==NULL)
					{
					newen->next=pres->calleds;
					pres->calleds=newen;
					}
				else
					{
					newen->next=prevr->next;
					prevr->next=newen;
					}
				freenew=0;
				break;
				}
			else if(!i)
				{/* called more than once by same routine*/
/*printf(" ignoring later call of %s  by %s\n",name,caller->name);*/
				break;
				}
			prevr=presr;
			presr=presr->next;
			}
if(freenew)free(newen);
end:
/*printf(" entersym returning, root=%s\n",root->name);*/
return el;
}

char *filein,*d,fileout[50];
FILE *in,*out,*keywdf, *o1,*o2;
struct ffblk ffblk;

main(argc,argv)
int argc;
char *argv[];
{
int hitas,nwhit,hitbr,index,done=0,a,cpp,igc,hitc,hitcp,i;
char *nm, kwd[bufflen+1];
struct taben *pres,*prev;struct ref *presr,*prevr;
filein= argv[1];
caller=root=NULL;
if(argc==1){help(); exit(0);}
/*printf(" argc=%d\n",argc);*/
igc=0;cpp=0;
if (argc<2)
	{
	printf(" enter filename pattern\n");
	scanf("%s",filein);
	}
else if(argc>2)
	{
	d=argv[2];
	printf(" Echo Options=%s...\n",d);
	while((*d)!=0)
		{
		  if(*d=='p'||*d=='P')cpp=1;
		}
	}
keywdf=fopen("keywd.dat","r");
if(keywd==NULL){printf(" cant find keyword list\n");return 1;}
for(kwct=0;kwct<maxkw;kwct++)
	{
	for(i=0;i<bufflen+1;i++)kwd[i]='\000';
	hitc=fscanf(keywdf,"%s",kwd);
	if(hitc==EOF)break;
	strncpy(keywd[kwct],kwd,bufflen);
	/*printf(" %s\n",keywd[kwct]);*/
	}
if(kwct+1== maxkw)printf(" hit limit %d\n",maxkw);
else printf(" %d keywords to ignore\n ",kwct);

done=findfirst( filein,&ffblk,0);
while (!done)
	{
/*process file */
printf(" processing file:%s\n",ffblk.ff_name);
	in=fopen(ffblk.ff_name,"rb");
	if (in==NULL){printf(" bad open input\n");exit(0);}
/* open it */
/*	printf(" creating %s\n",fileout);
	out=fopen(fileout,"wb");
	if(out==NULL) printf(" trouble opening %s\n",fileout);
*/
/* do it*/
	hitbr=0; index=0;nwhit=0;

	for(;;)
		{
		a=fgetc(in);
		if( a==EOF)
			{
			if(hitbr){printf(" eof with hitbr=%d\n",hitbr);break;}
			break;
			}
		else if(iswhite(a)||a==';'||a=='*')
			{/* close off symbol but don't throw away*/
			nwhit=0;
			}
		else if(a=='\'')
			{nwhit=0;index=0;
			/*gobble up but be careful about'*/
			for(;;)
				{a=fgetc(in);
				if(a=='\'')break;
				if(a=='\\')a=fgetc(in);
				}			
			}
		else if(a=='\"')
			{nwhit=0;index=0;
			/*gobble up but be careful about\"*/
			for(;;)
				{a=fgetc(in);
				if(a=='\"')break;
				if(a=='\\')a=fgetc(in);
				}
			}
		else if(a=='/')
			{nwhit=0;index=0;
			a=fgetc(in);
			if( a=='/' && cpp)
				{
				for(;;)
					{
					a=fgetc(in);
					if(a=='\n') break;
					}
				}
			else if(a!='*')ungetc(a,in);
			else
				{
				hitas=0;
				for(;;)
					{a=fgetc(in);
					if(a=='/'&&hitas)break;
					if(a=='*')hitas=1;
					else hitas=0;
					}			
				}
			}

		else if(	a=='='||
			a=='='||a=='-'||a=='%'||a=='+'
			||a=='&'||a=='&'||a=='~'||a=='?'
			||a=='<'||a=='>'||a=='!'
			||a=='^'||a==']'||a=='['
			||a=='|'||a=='\\'
			||a==','||a=='.'||a==':')
			{/* throw away*/
			nwhit=0;index=0;
			}
		else if(a=='(')
			{
			if(hitbr)
				{/* called*/
				cursym[index]='\000';
				dummy=entersym(cursym,1,index);
/*printf("------------------- CALLED:%s\n",cursym);*/
				index=0;nwhit=0;/* reset for next*/
				}
			else
				{/* caller*/
				cursym[index]='\000';
				dummy=entersym(cursym,0,index);
				if(dummy!=NULL) 
					{caller=dummy;
/*printf("------------------- caller: %s\n",caller->name);*/
					}
				index=0;nwhit=0;
				}
			}
		else if(a=='{')hitbr++;
		else if(a=='}'){hitbr--;
				/*if(!hitbr)caller=NULL;*/
				}
		else if(a==')'){index=0;nwhit=0;}
		else/* build symbol*/
			{
			if(!nwhit)index=0;
			cursym[index++]=a;nwhit=1;
			if(index>= MAXLEN)
				{
/*printf("debug warn:symbol overflow cleared\n");*/
				index=0;
				}
			}
		}   /* infinite loop over file data*/
fclose(in);
/* on to next file*/
	nextf: done=findnext(&ffblk);
	}     /* while !done*/
/* write out statistics*/
/* double loop outer: select each caller/called print caller list*/
/* list each function and all it calleds?*/
pres=root;
o1=fopen("callers","w");
/*fprintf(o1,"\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n");*/
for(;pres!=NULL;)
	{
	fprintf(o1,"\n function: %s called by:\n",pres->name);
	presr=pres->callers;
	for(;presr!=NULL;)
		{
		nm=(presr->caller)->name;
		if(nm!=NULL)fprintf(o1,"        %s\n",nm);
		presr=presr->next;
		}
	pres=pres->next;
	}

o2=fopen("callees","w");
/*fprintf(o2,"\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n");*/
pres=root;
for(;pres!=NULL;)
	{
	fprintf(o2,"\n function: %s calls:\n",pres->name);
	presr=pres->calleds;
	for(;presr!=NULL;)
		{
		nm=(presr->caller)->name;
		if(nm!=NULL)fprintf(o2,"        %s\n",nm);
		presr=presr->next;
		}
	pres=pres->next;
	}
exit(0);
return 0;
}/* end main*/

