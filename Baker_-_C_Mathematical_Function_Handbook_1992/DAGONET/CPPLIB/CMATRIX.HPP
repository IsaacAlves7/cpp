//matrix package

#include "cvector.hpp"

class Cmatrix  {

private:
	void init(int row=1,int col=1,int b=0);
public:
	complex** m;
	int rowkt,colkt,base;

	Cmatrix(int rowkt,int colkt,int b) 
		{ init( rowkt, colkt, b);} //constructor
	Cmatrix() {init();}	   //  "
	Cmatrix( Cmatrix&);	   // init
        ~Cmatrix();		   // destructor
	void operator=(Cmatrix& );
	inline complex& elemnt(int i,int j)
		{return m[i-base][j-base];}
	void check(int,int);
	complex& element(int,int);
	inline void setelemnt(int i,int j,complex value)
		{m[i][j]=value;	}
	void setelement(int,int, complex);
	friend Cmatrix operator*(Cmatrix&,Cmatrix&); //Mat*Mat	
	friend Cvector operator*(Cmatrix&,Cvector&); //Mat*Vector
};
