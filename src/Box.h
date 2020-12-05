#ifndef __Box_h
#define __Box_h
#include "Params.h"

typedef struct {
  double center_digits[6];
  double size_digits[6];
  double center[6];
  double size[6];
  ACJParams cover;
  XParams nearer; // all values closer to 0 than in box or 0 if box overlaps
  XParams further; // all values futher from 0 that in the box
  XParams greater; // all values greater than in the box
} Box;

Box build_box(char* where);

#endif // __Box_h
