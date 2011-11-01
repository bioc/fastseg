#include <Rdefines.h>

/* segment.cpp */
SEXP segment(SEXP xS, SEXP epsS, SEXP deltaS, SEXP maxIntS,
		SEXP minSegS, SEXP squashingS, SEXP cyberWeightS);

/* segmentCyberT.cpp */
SEXP segmentCyberT(SEXP xS, SEXP epsS, SEXP maxSDS, SEXP deltaS, SEXP maxIntS,
		SEXP minSegS, SEXP squashingS, SEXP cyberWeightS, SEXP minVarPercS);
