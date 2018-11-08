#ifndef __Box_h
#define __Box_h
#include "Params.h"

struct Box {
	Box(char* where);
	Params<ACJ> cover() const;
	Params<XComplex> nearer() const; // returns all values closer to 0 than in box or 0 if box overlaps
	Params<XComplex> further() const; // returns all values futher from 0 that in the box
	Params<XComplex> greater() const; // returns all values greater than in the box
private:
	double center_digits[6];
	double size_digits[6];
    double box_center[6];
    double box_size[6];
    void compute_center_and_size();
};

#endif // __Box_h
