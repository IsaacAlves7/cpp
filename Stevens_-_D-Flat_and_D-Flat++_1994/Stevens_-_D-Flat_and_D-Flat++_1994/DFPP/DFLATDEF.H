// -------- dflatdef.h

#ifndef DFLATDEF_H
#define DFLATDEF_H

// -------- window class types
enum WndType {
    DFlatWindow,
    ApplicationWindow,
    TextboxWindow,
	LabelWindow,
    FrameWindow,
    ScrollbarWindow,
    MenubarWindow,
    ListboxWindow,
    PopdownWindow,
    EditboxWindow,
    EditorWindow,
    DialogWindow,
    PushButtonWindow,
    RadioButtonWindow,
	CheckBoxWindow,
    ToolButtonWindow,
    ToolbarWindow,
    StatusbarWindow
};

enum Bool {False, True};

template<class T>
T min(T a, T b)
{
	return a < b ? a : b;
}

template<class T>
T max(T a, T b)
{
	return a > b ? a : b;
}

template<class T>
void swap(T& a, T& b)
{
	T tmp = a;
	a = b;
	b = tmp;
}

#endif

