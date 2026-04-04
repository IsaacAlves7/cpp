// -------- rectangl.h

#ifndef RECTANGL_H
#define RECTANGL_H

#include "dflatdef.h"

class Rect	{
    int left,top,right,bottom;
public:
	Rect() {};
	Rect(int l, int t, int r, int b)
		{ left = l; top = t; right = r; bottom = b; }
	virtual ~Rect() { }
	int Top()    const	{ return top; }
	int Bottom() const	{ return bottom; }
	int Left()	 const	{ return left; }
	int Right()  const	{ return right; }
	int Height() const	{ return bottom-top+1; }
	int Width()	 const	{ return right-left+1; }
	void SetTop(int tp)    { top = tp; }
	void SetBottom(int bt) { bottom = bt; }
	void SetLeft(int lf)   { left = lf; }
	void SetRight(int rt)  { right = rt; }
	
	Bool Inside(int x, int y) const
		{ return (Bool)(x >= left && x <= right && y >= top && y <= bottom); }
	Rect subRectangle(const Rect& rc) const;
	Bool ValidRectangle() const
		{ return (Bool) (right != -1 && top != -1); }
	Bool operator==(const Rect& rc) const
		{ return (Bool) (right == rc.right && top == rc.top &&
						 left == rc.left && bottom == rc.bottom); }
};

#endif



