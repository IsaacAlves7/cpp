// ---------- directry.cpp

#include <direct.h>
#include "directry.h"
#include "desktop.h"

void PathNameLabel::FillLabel()
{
    SetText(desktop.CurrentPath());
}


DriveListBox::DriveListBox(int lf, int tp, DFWindow *par)
                        : ListBox(lf, tp, 10, 10, par)
{
    SetAttribute(BORDER);
    currdrive = getdisk();
    union REGS regs;
    for (unsigned int dr = 0; dr < 26; dr++)    {
        setdisk(dr);
        if (getdisk() == dr)    {
            // ----- test for remapped B drive
            if (dr == 1)    {
                regs.x.ax = 0x440e; // IOCTL func 14
                regs.h.bl = dr+1;
                int86(0x21, &regs, &regs);
                if (regs.h.al != 0)
                    continue;
            }
            
            String drname(" :");
            drname[0] = dr+'A';

            // ---- test for network or RAM disk
            regs.x.ax = 0x4409;     // IOCTL func 9
            regs.h.bl = dr+1;
            int86(0x21, &regs, &regs);
            if (!regs.x.cflag)    {
                if (regs.x.dx & 0x1000)
                    drname += " (Net)";
                else if (regs.x.dx == 0x0800)
                    drname += " (RAM)";
            }
            AddText(drname);
        }
    }
    setdisk(currdrive);
    SetScrollBars();
}


DirectoryListBox::DirectoryListBox(int lf, int tp, DFWindow *par)
                        : ListBox(lf, tp, 10, 13, par)
{
    SetAttribute(BORDER | SORTED);
    FillList();
}

void DirectoryListBox::FillList()
{
    ClearText();
    int ax;
    struct ffblk ff;
	String directory;
    ax = desktop.FindFirst("*.*", &ff, FA_DIREC);
    while (ax == 0)    {
        if ((ff.ff_attrib & FA_DIREC) != 0)
            if (!(*ff.ff_name == '.' && *(ff.ff_name+1) == '\0'))
                directory += String(ff.ff_name) + "\n";
        ax = findnext(&ff);
    }
	SetText(directory);
    SetScrollBars();
}

FileNameListBox::FileNameListBox(const String& filespec,int lf,int tp, DFWindow *par)
                        : ListBox(lf, tp, 10, 14, par)
{
    SetAttribute(BORDER | SORTED);
    FillList(filespec);
}

void FileNameListBox::FillList(const String& filespec)
{
    ClearText();
    int ax;
    struct ffblk ff;
	String list;
	String fs = filespec;
	if (fs.Strlen() == 0)
		fs = "*.*";
    ax = desktop.FindFirst(fs, &ff, 0);
    while (ax == 0)    {
        list += String(ff.ff_name) + "\n";
        ax = findnext(&ff);
    }
	SetText(list);
    SetScrollBars();
}



