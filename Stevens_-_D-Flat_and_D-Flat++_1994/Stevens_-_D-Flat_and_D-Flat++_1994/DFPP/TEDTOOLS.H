// ---------- tedtools.h

#ifndef TEDTOOLS_H
#define TEDTOOLS_H

#include "toolbar.h"

// -------- the TED toolbar
class TedTools : public ToolBar    {
    ToolButton newtool;
    ToolButton opentool;
    ToolButton savetool;
public:
    TedTools(DFWindow *par);
};

#endif

