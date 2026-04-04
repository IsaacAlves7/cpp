// ---------- fileopen.cpp

#include <io.h>
#include <dir.h>
#include "fileopen.h"
#include "notice.h"

// ----------------- File Open Dialog Box
FileOpen::FileOpen(const String& spec, const String& ttl) :
        Dialog(ttl, 19, 57),
        // ----- file name editbox
        filelabel ("~Filename:", 3, 2, 10, this),
        filename  (13, 2, 1, 40, this),
        // ----- drive:path display
        dirlabel  (3, 4, 50, this),
        // ----- files list box
        fileslabel("F~iles:", 3, 6, 7, this),
        files     (spec, 3, 7, this),
        // ----- directories list box
        dirslabel ("~Directories:", 19, 6, 13, this),
        dirs      (19, 7, this),
        // ----- drives list box
        disklabel ("Dri~ves:", 34, 6, 8, this),
        disks     (34, 7, this),
        // ----- command buttons
        ok        (46, 8, this),
        cancel    (46,11, this),
        help      (46,14, this),
        // ------ file open data members
        filespec(spec)
{
	PrepFilename(filespec);
}

void FileOpen::PrepFilename(String& fname)
{
    filename.SetText(fname);
	int x = fname.Strlen();
	if (x > 0)	{
		filename.SelectText(0, 0, x, 0);
		filename.SetCursor(x, 0);
		int wd = filename.ClientWidth();
		if (x >= wd)
			filename.SetLeftMargin(x-wd+1);
	}
    filename.Paint();
}

// --- Get selected filename: files listbox->filename editbox
void FileOpen::SelectFileName()
{
    int sel = files.Selection();
    if (sel != -1)    {
        String fname;
		files.ExtractTextLine(fname, sel);
		PrepFilename(fname);
    }
}

// ---- called when user "selects" a control
//      e.g. changes the selection on a listbox
void FileOpen::ControlSelected(DFWindow *Wnd)
{
    if (Wnd == (DFWindow *) &files)
        // --- user selected a filename from list
        SelectFileName();
    else if (Wnd == (DFWindow *) &dirs ||
            Wnd == (DFWindow *) &disks)    {
        // --- user is selecting a different drive or directory
		PrepFilename(filespec);
    }
}

// ---- called when user "chooses" a control
//      e.g. chooses the current selection on a listbox
void FileOpen::ControlChosen(DFWindow *Wnd)
{
    if (Wnd == (DFWindow *) &files)
        // --- user chose a filename from filename list
        OKFunction();
    else if (Wnd == (DFWindow *) &dirs)    {
        // --- user chose a directory from directory list
		String dr;
		dirs.ExtractTextLine(dr, dirs.Selection());
        chdir(dr);
        ShowLists();
    }
    else if (Wnd == (DFWindow *) &disks)    {
        // --- user chose a drive from drive list
		String dr;
		disks.ExtractTextLine(dr, disks.Selection());
        if (desktop.SetDisk(dr[0]))
	        ShowLists();
    }
}

Bool FileOpen::GetFileName(String& fname)
{
    filename.ExtractTextLine(fname, 0);
	char drive[3] = "";
	char dir[65] = "";
	char name[9] = "";
	char ext[5] = "";
	int pt = fnsplit(fname, drive, dir, name, ext);

	if (pt & DRIVE)
        if (!desktop.SetDisk(drive[0]))
			return False;

	if (pt & DIRECTORY)	{
		int ln = strlen(dir);
		if (ln > 1)	{
			char *cp = dir+ln-1;
			if (*cp == '\\')
				*cp = '\0';
		}
		if (chdir(dir) == -1)	{
			ErrorMessage("Invalid path");
			return False;
		}
	}

	if (pt & FILENAME)	{
		fname = name;
		if (pt & EXTENSION)
			fname += ext;
		if (fname.FindChar('*') == -1 && fname.FindChar('?') == -1)	{
			if (chdir(fname) == -1)
				return True;
			// --- this is a subdirectory name
			fname = filespec;
		}
	}
    filespec = fname;
	PrepFilename(filespec);
    ShowLists();
	return False;
}

// ---- called when user chooses OK command
void FileOpen::OKFunction()
{
	String fname;
	if (GetFileName(fname) == True)	{
    	if (access(fname, 0) == 0)    {
        	filespec = fname;
        	Dialog::OKFunction();
    	}
    	else
        	// ---- No file as specified
        	ErrorMessage("File does not exist");
	}
}

// ------ refresh the current directory display and 
//        the directories and files list after
//        user changes filespec, drive, or directory
void FileOpen::ShowLists()
{
    dirlabel.FillLabel();
    dirlabel.Show();
    dirs.FillList();
    dirs.Show();
    files.FillList(filespec);
    files.Show();
}

// ------- called just before a control gets the focus
void FileOpen::EnterFocus(DFWindow *Wnd)
{
    if (Wnd == (DFWindow *) &files)
        // --- The file name list box is getting the focus
        SelectFileName();
	Dialog::EnterFocus(Wnd);
}

SaveAs::SaveAs(const String& spec, const String& fname, const String& ttl)
        : FileOpen(spec, ttl)
{
	PrepFilename(String(fname));
}

// ---- called when user chooses OK command
void SaveAs::OKFunction()
{
	String fname;
	if (GetFileName(fname) == True)	{
    	if (access(fname, 0) != 0)    {
        	// ---- chosen file does not exist
           	filespec = fname;
           	Dialog::OKFunction();
    	}
    	else    {
        	// ---- file exists
        	String msg = fname + " already exists. Replace?";
        	if (YesNo(msg))    {
            	filespec = fname;
            	Dialog::OKFunction();
        	}
    	}
	}
}

