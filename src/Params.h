#ifndef __Params_h
#define __Params_h
#include <math.h>
#include "SL2ACJ.h"

typedef struct {
  XComplex lattice;
  XComplex loxodromic_sqrt;
  XComplex parabolic;
} XParams;

typedef struct {
  ACJ lattice;
  ACJ loxodromic_sqrt;
  ACJ parabolic;
} ACJParams;

SL2ACJ construct_G(const ACJParams& params);
ACJ construct_T(const ACJParams& params, int M, int N);
SL2ACJ construct_word(const ACJParams& params, const char* word);

#endif // __Params_h
