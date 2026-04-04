// --------- fileopen.h

#ifndef FILEOPEN_H
#define FILEOPEN_H

#include "dialog.h"
#include "directry.h"
#include "editbox.h"
#include "label.h"
#include "pbutton.h"

// ------------ File Open dialog box
class FileOpen : public Dialog    {
    // -----File Open Dialog Box Controls:
protected:
    // ----- file name editbox
    Label filelabel;
    EditBox filename;
    // ----- drive:path display
    PathNameLabel dirlabel;
    // ----- files list box
    Label fileslabel;
    FileNameListBox files;
    // ----- directories list box
    Label dirslabel;
    DirectoryListBox dirs;
    // ----- drives list box
    Label disklabel;
    DriveListBox disks;
    // ----- command buttons
    OKButton ok;
    CancelButton cancel;
    HelpButton help;
    // ------ file open data members
    String filespec;
    // ------ private member functions
    void SelectFileName();
    void ShowLists();
	void PrepFilename(String& fname); // (called from ctor: cannot be virtual)
    // --- functions inherited from DFWindow
    virtual void ControlSelected(DFWindow *Wnd);
    virtual void ControlChosen(DFWindow *Wnd);
    virtual void EnterFocus(DFWindow *Wnd);
    virtual void OKFunction();
	Bool GetFileName(String& fname);
public:
    FileOpen(const String& spec = "*.*", const String& ttl = "File Open");
    const String& FileName() { return filespec; }
};

class SaveAs : public FileOpen    {
    virtual void OKFunction();
public:
    SaveAs(const String& spec = "",
		   const String& fname = "",
		   const String& ttl = "Save As");
};

#endif

