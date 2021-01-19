#ifndef _params_h_
#define _params_h_
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

#endif // _params_h_
