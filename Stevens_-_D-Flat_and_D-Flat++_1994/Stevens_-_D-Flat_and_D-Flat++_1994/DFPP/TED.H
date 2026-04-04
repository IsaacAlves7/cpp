// --------- ted.h

#ifndef TED_H
#define TED_H

#include "applicat.h"
#include "fileopen.h"
#include "tedtools.h"
#include "editor.h"
#include "notice.h"

#define Df void (DFWindow::*)()
#define Ap void (Application::*)()

extern MenuBarItem TedMenu[];
extern MenuSelection InsertCmd, WordWrapCmd, ToolbarCmd, RulerCmd;
extern MenuSelection CutCmd, CopyCmd, PasteCmd;

// ------- Ted application definition
class Ted : public Application {
    MenuBar menubar;
    TedTools toolbar;
    String  fname;
    Editor editor;
	Bool ruler;
protected:
    virtual void Move(int x, int y);
    virtual void Size(int x, int y);
	virtual void Paint();
    void BuildTitle();
    void TestChanged();
    void CloseWindow();
	void Border();
public:
    Ted();
	void OpenFile(const String& fn);
    // ----- menu command functions
    void CmNew();
    void CmOpen();
    void CmSave();
    void CmSaveAs();
    void CmInsert();
    void CmCut();
    void CmCopy();
    void CmPaste();
    void CmPara();
	void CmToolbar();
	void CmRuler();
	void CmTabs();
    void CmExit()   { CloseWindow(); }
	void PrepEdit();
	void Execute();
};

#define MinimumWidth  40
#define MinimumHeight 10

#endif
