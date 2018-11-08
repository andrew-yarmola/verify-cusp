#include "Box.h"
#include <stdlib.h>

double scale[6];
static bool scale_initialized = false; 
Box::Box(char * where) {
	if (!scale_initialized) {
		scale_initialized = true;
		for (int i = 0; i < 6; ++i) {
			scale[i] = pow(2, -i / 6.0);
		}
	}
	for (int i = 0; i < 6; ++i) {
		center_digits[i] = 0;
		size_digits[i] = 8;
	}
    int pos = 0;
    int idx = 0;
    int dir;
    while (where[idx] != '\0') {
        if (where[idx] == '0') {
            dir = 0;
        } else if (where[idx] == '1') {
            dir = 1;
        } else {
            fprintf(stderr, "verify: fatal error at %s\n", where);
            exit(1);
        }
	    size_digits[pos] *= 0.5;
	    center_digits[pos] += (2*dir-1)*size_digits[pos];
	    ++pos;
        if (pos == 6) {
            pos = 0;
        }
	    ++idx;
    }
    compute_center_and_size();
}

void Box::compute_center_and_size()
{
	for (int i = 0; i < 6; ++i) {
        // GMT paper page 419 of Annals
        // box_size guarantees that :
        // box_center - box_size <= true_center - true_size
        // box_center + box_size >= true_center + true_size
        // where box operations are floating point. 
        box_center[i] = scale[i]*center_digits[i];
        box_size[i]= (1+2*EPS)*(size_digits[i]*scale[i]+HALFEPS*fabs(center_digits[i]));
    }
}

Params<ACJ> Box::cover() const
{
	Params<ACJ> result;
	result.lattice = ACJ(
		XComplex(box_center[3], box_center[0]),
		XComplex(box_size[3], box_size[0]),
		0.,
		0.
	);
	result.loxodromic_sqrt = ACJ(
		XComplex(box_center[4], box_center[1]),
		0.,
		XComplex(box_size[4], box_size[1]),
		0.
	);
	result.parabolic = ACJ(
		XComplex(box_center[5], box_center[2]),
		0.,
		0.,
		XComplex(box_size[5], box_size[2])
	);
	return result;
}

Params<XComplex> Box::nearer() const
{
	double m[6];
	for (int i = 0; i < 6; ++i) {
        m[i] = 0; // inconclusive cases
        if (center_digits[i] > 0 && // center is positive 
            center_digits[i] > size_digits[i] &&  // true diff is positive
            box_center[i]    > box_size[i]) { // machine diff is >= 0
            // Want lower bound on true_center - true_size.  Assume no overflow or underflow 
            // Note, sign(center_digits) == sign(box_center), unless box_center == 0. Also, box_size is always >= 0. 
            // GMT paper page 419 of Annals gives with true arithmetic
            //      box_center - box_size <= true_center - true_size
            // Now, in machine arthimetric, by IEEE, if 
            //      box_center > box_size then box_center (-) box_size >= 0.
            // Lemma 7.0 gives,
            //      (1-EPS)(*)( box_center (-) box_size ) <= box_center - box_size <= true_center - box_size. 
            m[i] = (1-EPS)*(box_center[i] - box_size[i]);
        } else if (center_digits[i] < 0 && // center is negative
                   center_digits[i] < -size_digits[i] && // true sum is negative
                   box_center[i]    < -box_size[i]) {  // machine sum is negative
            // Want upper bound on true_center - true_size.  Assume no overflow or underflow
            // Note, sign(center_digits) == sign(box_center), unless box_center == 0. Also, box_size is always >= 0. 
            // GMT paper page 419 of Annals gives with true arithmetic
            //      true_center + true_size <= box_center + box_size
            // Now, in machine arthimetric, by IEEE, if 
            //      -box_center > box_size then (-box_center) (-) box_size >= 0.
            // Lemma 7.0 gives,
            //      (1-EPS)(*)( (-box_center) (-) box_size ) <= -box_center - box_size <= -true_center - true_size.
            // So,
            //      -((1-EPS)(*)( (-box_center) (-) box_size )) >= true_center + true_size.
            // Note, negation is exact for machine numbers
            m[i] = -((1-EPS)*((-box_center[i]) - box_size[i]));
        }
	}
	
	Params<XComplex> result;
	result.lattice = XComplex(m[3], m[0]);
	result.loxodromic_sqrt = XComplex(m[4], m[1]);
	result.parabolic = XComplex(m[5], m[2]);
	return result;
}

Params<XComplex> Box::further() const
{
	double m[6];
	for (int i = 0; i < 6; ++i) {
        m[i] = 0; // inconclusive cases
		if (center_digits[i] > -size_digits[i]) { // true sum is positive 
            // Want upper bound of true_center + true_size. Assume no overflow or underflow
            // Note, sign(center_digits) == sign(box_center), unless box_center == 0. Also, box_size is always >= 0. 
            // GMT paper page 419 of Annals gives with true arithmetic
            //      true_center + true_size <= box_center + box_size
            // By IEEE (+) and (-) resepct <= and >=, so box_center (+) box_size >=0 and
            // Lemma 7.0 for floating point arithmetic gives and upper bound
            //      (1+EPS)(*)(box_center (+) box_size) >= box_center + box_size >= true_center + true_size
		    m[i] = (1+EPS)*(box_center[i] + box_size[i]);
        } else { // true sum is <= 0
            // Want lower bound of true_center - true_size. Assume no overflow or underflow
            // Note, sign(center_digits) == sign(box_center), unless box_center == 0 
            // GMT paper page 419 of Annals gives with true arithmetic
            //      box_center - box_size <= true_center - true_size
            // By IEEE, (+) and (-) respects <= and >=, and negation is exact.
            // Thus, (-box_center) (+) box_size >=0 and Lemma 7.0 for floating point arithmetic gives
            //        (1+EPS)(*)( (-box_center) (+) box_size) ) >= (-box_center) + box_size
            // So,
            //      -((1+EPS)(*)( (-box_center) (+) box_size) ))<= box_center - box_size <= true_center - true_size
            m[i] = -((1+EPS)*((-box_center[i]) + box_size[i]));
        }
	}
	
	Params<XComplex> result;
	result.lattice = XComplex(m[3], m[0]);
	result.loxodromic_sqrt = XComplex(m[4], m[1]);
	result.parabolic = XComplex(m[5], m[2]);
	return result;
}

Params<XComplex> Box::greater() const
{
	double m[6];
	for (int i = 0; i < 6; ++i) {
        m[i] = 0; // inconclusive cases
		if (center_digits[i] > -size_digits[i]) { // true sum is positive
            // Want upper bound of true_center + true_size. Assume no overflow or underflow
            // Note, sign(center_digits) == sign(box_center), unless box_center == 0. Also, box_size is always >= 0. 
            // GMT paper page 419 of Annals gives with true arithmetic
            //      true_center + true_size <= box_center + box_size.
            // Notice that box_center + box_size >= true_center + true_size > 0.
            // By IEEE, box_center (+) box_size >=0, as it's guanrateed to evaluate to nearest representable.
            // Lemma 7.0 for floating point arithmetic gives and upper bound
            //      (1+EPS)(*)(box_center (+) box_size) >= box_center + box_size >= true_center + true_size
		    m[i] = (1+EPS)*(box_center[i] + box_size[i]);
        } else if (center_digits[i] < -size_digits[i] && // true sum is negative
                   box_center[i]    < -box_size[i]) { // machine sum is <= 0
            // Want upper bound of true_center + true_size. Assume no overflow or underflow
            // Note, sign(center_digits) == sign(box_center), unless box_center == 0. Also, box_size is always >= 0. 
            // GMT paper page 419 of Annals gives with true arithmetic
            //      true_center + true_size <= box_center + box_size.
            // Notice that box_center + box_size < 0.
            // By IEEE, box_center (+) box_size <= 0, as it's guanrateed to evaluate to nearest representable.
            // Lemma 7.0 for floating point arithmetic gives a bound
            //      (1-EPS)(*)| box_center (+) box_size | < | box_center + box_size |
            // So,
            //      -((1-EPS)(*)(-(box_center (+) box_size))) >= box_center + box_size >= true_center + true_size
            m[i] = -((1-EPS)*(-(box_center[i] + box_size[i])));
        }
	}
	
	Params<XComplex> result;
	result.lattice = XComplex(m[3], m[0]);
	result.loxodromic_sqrt = XComplex(m[4], m[1]);
	result.parabolic = XComplex(m[5], m[2]);
	return result;
}
