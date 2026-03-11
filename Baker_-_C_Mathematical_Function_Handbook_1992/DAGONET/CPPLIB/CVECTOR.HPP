// Cvector.hpp
#include <stdlib.h>
#include "Complex.hpp"
static complex *present;

class Cvector {
	protected:
	public:
		int size,base;
		complex *head;
Cvector::Cvector( int s=1 , int b=0 , complex initvalue= (0.,0.) )
			//constructor
			{
			head= new complex[s];
			base=b;
			  for(size=0,present=head;size<s;
				size++,present++) *present = initvalue; 
			  size=s;		 
//			 printf(" Cvector built\n");
			}
		Cvector( Cvector&);	//copy
		~Cvector()//destructor
			{
			delete head; //printf(" Cvector axed\n");
			}
		void operator=( Cvector& rhs);
		complex operator*(Cvector& rvalue);//dot product
		inline int length()
			{return size;}
		inline complex& elemnt( int i)
			{return head[i-base];}
		inline void setelemnt( int i, complex& value)
			{head[i-base]=value;}
		void check(int);
		complex& element(int i);
		void setelement(int i, complex& value);
	};
