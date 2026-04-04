// -------- strings.cpp

#include "strings.h"

// ------- put a string into a String
void String::putstr(const char *s)
{
	delete sptr;
	sptr = 0;
	if (s)	{
		sptr = new char[strlen(s)+1];
		strcpy(sptr, s);
	}
}

// --- convert from char *
String::String(const char *s)
{
	sptr = 0;
	putstr(s);
}

// ------- copy constructor
String::String(String& s)
{
	sptr = 0;
	putstr(s.sptr);
}

// -------- construct with a size and fill character
String::String(int len, char fill)
{
	sptr = new char[len+1];
	memset(sptr, fill, len);
	*(sptr+len) = '\0';
}

// -------- assignment
String& String::operator=(const String& s)
{
	if (this != &s)
		putstr(s.sptr);
	return *this;
}

// ------- concatenation operator (str1 + str2)
String String::operator+(const String& s) const
{
	int ln = 0;
	if (sptr)
		ln += strlen(sptr);
	if (s.sptr)
		ln += strlen(s.sptr);
	String temp(ln, 0);
	if (sptr)
		strcpy(temp.sptr, sptr);
	if (s.sptr)
		strcat(temp.sptr, s.sptr);
	return temp;
}

// --- concatenation operator (str1 += str2;)
void String::operator+=(const String& s)
{
	String temp = *this;
	temp = temp + s;
	*this = temp;
}

// --- concatenation operator (str1 + "str")
String String::operator+(const char *s) const
{
	String s1(s);
	String temp = *this;
	return temp + s1;
}

// --- concatenation operator (str1 += "str")
void String::operator+=(const char *s)
{
	String s1(s);
	String temp = *this + s1;
	*this = temp;
}

// ------ substring: right len chars
String String::right(int len) const
{
	String tmp;
	if (sptr && len)	{
		int sln = strlen(sptr);
		if (len > sln)
			len = sln;
		if (len)
			tmp = String(sptr + sln - len);
	}
	return tmp;
}

// ------ substring: left len chars
String String::left(int len) const
{
	String tmp(len, 0);
	if (sptr && len)	{
		int sln = strlen(sptr);
		if (len > sln)
			len = sln;
		if (len)
			strncpy(tmp.sptr, sptr, len);
	}
	return tmp;
}

// ------ substring: middle len chars starting from where
String String::mid(int len, int where) const
{
	String tmp(len, 0);
	if (len && sptr && where < (int) strlen(sptr))
		strncpy(tmp.sptr,sptr+where,len);
	return tmp;
}

// ---- find offset to first instance of specified char
int String::FindChar(unsigned char ch) const
{
	if (sptr)	{
		char *cp = strchr(sptr, ch);
		if (cp)
			return (int) (cp - sptr);
	}
	return -1;
}

// ------ change the buffer length of a string
void String::ChangeLength(unsigned int newlen)
{
	if (sptr == 0 || *sptr == '\0')
		*this = String(newlen, 0);
	else	{
		char *cp = new char[newlen+1];
		memset(cp, '\0', newlen+1);
		if (sptr != 0)	{
			strncpy(cp, sptr, newlen);
			delete sptr;
		}
		sptr = cp;
	}
}

// ------- stream I/O
ostream& operator<< (ostream& os, const String& str)
{
	if (str.sptr)
		os << str.sptr;
	return os;
}

istream& operator>> (istream& is, const String& str)
{
	if (str.sptr)	{
		*str.sptr = '\0';
		while (*str.sptr == '\0' && !is.eof())
			is.getline(str.sptr, strlen(str.sptr)+1);
	}
	return is;
}



