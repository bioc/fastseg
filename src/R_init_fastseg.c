#include "fastseg.h"
#include <R_ext/Rdynload.h>

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

static const R_CallMethodDef callMethods[] = {

/* segment.cpp */
	CALLMETHOD_DEF(segment, 7),

/* segmentCyberT.cpp */
	CALLMETHOD_DEF(segmentCyberT, 9),

	{NULL, NULL, 0}
};


void R_init_fastseg(DllInfo *info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	return;
}
