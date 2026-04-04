// --------- notice.h

#ifndef NOTICE_H
#define NOTICE_H

#include "dialog.h"
#include "pbutton.h"
#include "label.h"
#include "editbox.h"

// ----------- Notice message window
class NoticeBox : public Dialog	{
	Label notice;
	// ----- command button
public:
	NoticeBox(const String& text, const String& ttl = "", int ht = 5);
};

// ----------- Notice message window
class NoticeOKBox : public NoticeBox	{
	// ----- command button
	OKButton ok;
public:
	NoticeOKBox(const String& text, const String& ttl);
};

inline void NoticeOK(const String& txt, const String& ttl = "Notice")
{
	NoticeOKBox(txt, ttl).Execute();
}

// --------- Error message window
class ErrorMessageBox : public NoticeBox	{
public:
	ErrorMessageBox(const String& text, const String& ttl);
	virtual void Execute();
};

inline void ErrorMessage(const String& txt, const String& ttl = "Error")
{
	ErrorMessageBox(txt, ttl).Execute();
}

// ------- Yes/No query window
class YesNoBox : public Dialog	{
	Label notice;
	// ----- command buttons
	YesButton yes;
	NoButton  no;
    virtual void Keyboard(int key);
public:
	YesNoBox(const String& txt, const String& ttl);
};

inline Bool YesNo(const String& txt, const String& ttl = "Question")
{
	YesNoBox yn(txt, ttl);
	yn.Execute();
	return yn.OKExit();
}

// ----- Input Box
class InputMessageBox : public Dialog	{
	Label notice;
	EditBox editor;
	// ----- command buttons
	OKButton ok;
	CancelButton cancel;
	String& inputtext;
    virtual void OKFunction();
public:
	InputMessageBox(int len, const String& txt, String& msg, const String& ttl);
	EditBox& Editor() { return editor; }
};

inline Bool InputBox(int len, const String& txt, String& msg,
					int maxlen = 0, const String& ttl = "Input")
{
	InputMessageBox inbx(len, txt, msg, ttl);
	inbx.Editor().SetTextLength(maxlen);
	inbx.Execute();
	return inbx.OKExit();
}


#endif
