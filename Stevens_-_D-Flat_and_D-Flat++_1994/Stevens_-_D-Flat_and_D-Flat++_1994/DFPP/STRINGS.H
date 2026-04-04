// -------- strings.h

#ifndef STRINGS_H
#define STRINGS_H

#include <iostream.h>
#include <string.h>

// ============================
// BASIC-like String Class
// ============================
class String    {
	char *sptr;
	void putstr(const char *s);
public:
	// -------- construct a null string
	String() { sptr = 0; }
	// --- construct with char * initializer
	String(const char *s);
	// ------- copy constructor
	String(String& s);
	// -------- construct with a size and fill character
	String(int len, char fill);
	// ------- destructor
	~String() { delete sptr; }
	// ------ return the length of a string
	int Strlen() const { return strlen(sptr); }
	// ------ change the buffer length of a string
	void ChangeLength(unsigned int newlen);
	// ---- substring: right len chars
	String right(int len) const;
	// ---- substring: left len chars
	String left(int len) const;
	// ---- substring: middle len chars starting from where
	String mid(int len, int where) const;
	int FindChar(unsigned char ch) const;
	// ---------- assignment
	String& operator=(const String& s);
	// ---------- conversion to char *
	operator const char *() const { return sptr; }
	// --- concatenation operator (str1 + str2;)
	String operator+(const String& s) const;
	// --- concatenation operator (str1 += str2;)
	void operator+=(const String& s);
    // --- concatenation operator (str1 + "str")
    String operator+(const char *s) const;
	// --- concatenation operator (str1 += "str")
	void operator+=(const char *s);

	// ------- relational operators
	int operator==(const String& s) const
		{ return (stricmp(sptr,s.sptr) == 0); }
	int operator!=(const String& s) const
		{ return (stricmp(sptr,s.sptr) != 0); }
	int operator>(const String& s) const
		{ return (stricmp(sptr,s.sptr) > 0); }
	int operator<(const String& s) const
		{ return (stricmp(sptr,s.sptr) < 0); }
	int operator<=(const String& s) const
		{ return (!(*this > s)); }
	int operator>=(const String& s) const
		{ return (!(*this < s)); }
	// ------- subscript
	char& operator[](int n) const
		{ return sptr[n]; }
	// ------- stream I/O
	friend ostream& operator<< (ostream& os, const String& str);
	friend istream& operator>> (istream& is, const String& str);
};

#endif

