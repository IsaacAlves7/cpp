// ------------- rectangl.cpp

#include "rectangl.h"

// ---- is a point within a vector?
static inline Bool within(int p, int v1, int v2)
{
	return (Bool) (p >= v1 && p <= v2);
}

// --- Produce vector end points of the overlap of two other vectors
static void subVector(int& v1, int& v2,
                        int t1, int t2, int o1, int o2)
{
    v1 = v2 = -1;
    if (within(o1, t1, t2))    {
        v1 = o1;
        if (within(o2, t1, t2))
            v2 = o2;
        else
            v2 = t2;
    }
    else if (within(o2, t1, t2))    {
        v2 = o2;
        if (within(o1, t1, t2))
            v1 = o1;
        else
            v1 = t1;
    }
    else if (within(t1, o1, o2))    {
        v1 = t1;
        if (within(t2, o1, o2))
            v2 = t2;
        else
            v2 = o2;
    }
    else if (within(t2, o1, o2))    {
        v2 = t2;
        if (within(t1, o1, o2))
            v1 = t1;
        else
            v1 = o1;
    }
}

// --- Return rectangle produced by the overlap of two other rectangles
Rect Rect::subRectangle(const Rect& rc) const
{
    Rect r(0,0,0,0);
    subVector(r.left, r.right, left, right, rc.left, rc.right);
    subVector(r.top, r.bottom, top, bottom, rc.top, rc.bottom);
    return r;
}


