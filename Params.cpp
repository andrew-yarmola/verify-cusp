#include "SL2ACJ.h"
#include "Params.h"
#include <string.h> 
#include <stdio.h> 

SL2ACJ construct_G(const Params<ACJ>& params)
{
	ACJ I(XComplex(0., 1.));
	const ACJ& sl = params.loxodromic_sqrt;
	return SL2ACJ(I*params.parabolic*sl, I/sl, I*sl, ACJ(0.));
}

SL2ACJ construct_T(const Params<ACJ>& params, int x, int y)
{
	return SL2ACJ(ACJ(XComplex(1.)), params.lattice * double(y) + double(x), ACJ(XComplex(0.)), ACJ(XComplex(1.)));
}

SL2ACJ construct_word(const Params<ACJ>& params, char* word)
{
	ACJ one(1.), zero(0.);
	SL2ACJ w(one, zero, zero, one);
	SL2ACJ G(construct_G(params));
	SL2ACJ g(inverse(G));

	int x = 0;
	int y = 0;
	size_t pos;
	for (pos = strlen(word); pos > 0; --pos) {
		char h = word[pos-1];
		switch(h) {
			case 'm': --x; break;
			case 'M': ++x; break;
			case 'n': --y; break;
			case 'N': ++y; break;
			default: {
				if (x != 0 || y != 0) {
					w = construct_T(params, x, y) * w;
					x = y = 0;
				}
				if (h == 'g')
					w = g * w;
				else if (h == 'G')
					w = G * w;
			}
		}
	}
    if (x != 0 || y != 0) {
        w = construct_T(params, x, y) * w;
    }
    return w;
}

