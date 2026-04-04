// ------------- family.h

// ----- template of single-parent family tree
//
//  parent
//     child1, child2, ... childn
//

#ifndef FAMILY_H
#define FAMILY_H

template <class T>
class Family {
    T *parent;
    // -------- children
    T *firstchild;
    T *lastchild;
    // -------- siblings
    T *nextsibling;
    T *prevsibling;
public:
	Family();
	void SetParent(T *par)
		{ parent = par; }
    T *Parent()        { return parent; }
    T *FirstChild()    { return firstchild; }
	T *LastChild()     { return lastchild; }
	T *NextSibling()   { return nextsibling; }
	T *PrevSibling()   { return prevsibling; }
	void Append(T *child);  // append to parent's child list
	void Remove(T *child);  // remove from  "       "    "
};

template <class T>
Family<T>::Family()
{
    parent = 0;
    firstchild = 0;
    lastchild = 0;
    nextsibling = 0;
    prevsibling = 0;
}

template <class T>
void Family<T>::Append(T *child)
{
	nextsibling = 0;
	prevsibling = 0;
	if (parent != 0)	{
    	if (parent->family.firstchild == 0)
        	parent->family.firstchild = child;
		prevsibling = parent->family.lastchild;
    	if (prevsibling != 0)
			prevsibling->family.nextsibling = child;
		parent->family.lastchild = child;
	}
}

template <class T>
void Family<T>::Remove(T *child)
{
	if (parent != 0)	{
    	if (prevsibling != 0)
			prevsibling->family.nextsibling = nextsibling;
    	if (nextsibling != 0)
			nextsibling->family.prevsibling = prevsibling;
    	if (parent->family.firstchild == child)
			parent->family.firstchild = nextsibling;
    	if (parent->family.lastchild == child)
			parent->family.lastchild = prevsibling;
	}
	nextsibling = 0;
	prevsibling = 0;
}

#endif
