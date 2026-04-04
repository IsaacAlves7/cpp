// ---------- tedtools.cpp

#include "ted.h"
#include "tedtools.h"

// -------- the TED toolbar
TedTools::TedTools(DFWindow *par) : ToolBar(par),
            newtool("New",   (DFWindow *) this),
            opentool("Open", (DFWindow *) this),
            savetool("Save", (DFWindow *) this)
{
    newtool.SetButtonFunction(this->Parent(),
        (Df) (Ap) &Ted::CmNew);
    opentool.SetButtonFunction(this->Parent(),
        (Df) (Ap) &Ted::CmOpen);
    savetool.SetButtonFunction(this->Parent(),
        (Df) (Ap) &Ted::CmSave);
}


