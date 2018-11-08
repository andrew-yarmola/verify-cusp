#ifndef __Params_h
#define __Params_h
#include <math.h>
#include "SL2ACJ.h"

template<class N> struct Params {
	N lattice;
	N loxodromic_sqrt;
	N parabolic;
};

SL2ACJ construct_G(const Params<ACJ>& params);
ACJ construct_T(const Params<ACJ>& params, int M, int N);
SL2ACJ construct_word(const Params<ACJ>& params, const char* word);

#endif // __Params_h
