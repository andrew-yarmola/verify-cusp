#include "SL2ACJ.h"
#include "Params.h"
#include <string.h> 
#include <stdio.h> 

SL2ACJ construct_G(const ACJParams& params)
{
  ACJ I(XComplex(0, 1));
  const ACJ& sl = params.loxodromic_sqrt;
  return SL2ACJ(I*params.parabolic*sl, I/sl, I*sl, ACJ(0));
}

ACJ construct_T(const ACJParams& params, int M, int N)
{
  return params.lattice * double(N) + double(M);
}

SL2ACJ construct_word(const ACJParams& params, const char* word)
{
  ACJ one(1), zero(0);
  SL2ACJ w(one, zero, zero, one);
  SL2ACJ G(construct_G(params));
  SL2ACJ g(inverse(G));

  int M = 0;
  int N = 0;
  ACJ T;
  size_t pos;
  for (pos = strlen(word); pos > 0; --pos) {
    char h = word[pos-1];
    switch(h) {
      case 'm': --M; break;
      case 'M': ++M; break;
      case 'n': --N; break;
      case 'N': ++N; break;
      default: {
                 if (M != 0 || N != 0) {
                   T = construct_T(params, M, N);
                   w = SL2ACJ(w.a + T * w.c, w.b + T * w.d, w.c, w.d);
                   M = N = 0;
                 }
                 if (h == 'g')
                   w = g * w;
                 else if (h == 'G')
                   w = G * w;
               }
    }
  }
  if (M != 0 || N != 0) {
    T = construct_T(params, M, N);
    w = SL2ACJ(w.a + T * w.c, w.b + T * w.d, w.c, w.d);
  }
  return w;
}

