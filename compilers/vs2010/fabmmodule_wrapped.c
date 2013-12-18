#ifdef _DEBUG
	#define _SAVE_DEBUG
	#undef _DEBUG
#endif

#include <stdarg.h>
#include "Python.h"
#include "fortranobject.h"
#include <math.h>

// Below we redefine _DEBUG (see above)
#ifdef _SAVE_DEBUG
	#define _DEBUG
	#undef _SAVE_DEBUG
#endif

#include "fabmmodule.c"