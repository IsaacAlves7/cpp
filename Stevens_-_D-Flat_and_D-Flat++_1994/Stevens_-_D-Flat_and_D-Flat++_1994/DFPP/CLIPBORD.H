// --------- clipbord.h

#ifndef CLIPBORD_H
#define CLIPBORD_H

#include "strings.h"

class Clipboard    {
	String clipboard;
public:
	const String& GetText() const
		{ return clipboard; }
	void SetText(const String& txt)
		{ clipboard = txt; }
};

#endif



