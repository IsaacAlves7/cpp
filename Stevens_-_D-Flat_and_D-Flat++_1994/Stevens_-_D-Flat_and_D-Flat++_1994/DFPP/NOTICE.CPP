// ------------- notice.cpp

#include <ctype.h>
#include "notice.h"

static Color col = {
	BLACK,			// fg
	GREEN,			// bg
	BLACK,			// selected fg
	GREEN,			// selected bg
	BLACK,			// frame fg
	GREEN,			// frame bg
	BLACK,			// highlighted fg
	GREEN			// highlighted bg
};

NoticeBox::NoticeBox(const String& txt, const String& ttl, int ht) :
		Dialog(ttl, ht, 20),
		notice(txt, 1, 2, (DFWindow *)this)
{
	unsigned int wd = max(10,max(txt.Strlen()+2,ttl.Strlen()+2));
	Size(Left()+wd-1, Bottom());
	SetColor(col);
}

NoticeOKBox::NoticeOKBox(const String& txt, const String& ttl) :
		NoticeBox(txt, ttl, 7),
		ok(0, 4, (DFWindow *)this)
{
	ok.Move(Left()+((Width()-8)/2), ok.Top());
}

static Color ecol = {
	YELLOW, RED,
	YELLOW, RED,
	YELLOW, RED,
	YELLOW, RED
};

ErrorMessageBox::ErrorMessageBox(const String& txt, const String& ttl) :
		NoticeBox(txt, ttl)
{
	SetColor(ecol);
}

void ErrorMessageBox::Execute()
{
	desktop.speaker().Beep();
	NoticeBox::Execute();
}

static Color yncol = {
	BLACK,GREEN,
	BLACK,GREEN,
	BLACK,GREEN,
	BLACK,GREEN
};

YesNoBox::YesNoBox(const String& txt, const String& ttl) :
		Dialog(ttl, 7, 20),
		notice(txt, 1, 2, (DFWindow *)this),
		yes(0, 4, (DFWindow *)this),
		no(0, 4, (DFWindow *)this)
{
	unsigned int wd = max(21,max(txt.Strlen()+2,ttl.Strlen()+2));
	Size(Left()+wd-1, Bottom());
	yes.Move(Left()+((wd-17)/2), yes.Top());
	no.Move(Left()+((wd-17)/2)+9, no.Top());
	SetColor(yncol);
}

void YesNoBox::Keyboard(int key)
{
	switch (tolower(key))	{
		case 'n':
			CancelFunction();
			break;
		case 'y':
			OKFunction();
			break;
		default:
			Dialog::Keyboard(key);
			break;
	}
}

InputMessageBox::InputMessageBox(int len, const String& txt,
							 String& msg, const String& ttl) :
		Dialog(ttl, 7, 20),
		notice(txt, 1, 2, (DFWindow *)this),
		editor(1, 2, (DFWindow *)this),
		ok(0, 4, (DFWindow *)this),
		cancel(0, 4, (DFWindow *)this),
		inputtext(msg)
{
	int wd = len;
	wd += txt.Strlen();
	unsigned int tlen = wd;
	wd = max(wd,ttl.Strlen());
	wd = max(21, wd+2);
	Size(Left()+wd-1, Bottom());
	notice.Move(Left()+(wd-tlen)/2,notice.Top());
	editor.Move(notice.Right()+2, notice.Top());
	editor.Size(editor.Right()+len-2, editor.Top());
	editor.SetText(msg);
	ok.Move(Left()+((wd-17)/2), ok.Top());
	cancel.Move(Left()+((wd-17)/2)+9, cancel.Top());
	SetColor(col);
}

void InputMessageBox::OKFunction()
{
	inputtext = editor.GetText();
	Dialog::OKFunction();
}


