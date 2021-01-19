#include <stdio.h> 
#include <stdlib.h>
#include <string.h> 
#include "elimination.h"

extern double g_max_area;

// Helper functions

void check(bool inequalities, char* where)
{
    if (!inequalities) {       
        fprintf(stderr, "Fatal: verify error at %s\n", where);
        exit(3);
    }
}

void parse_word(char* code)
{
    char buf[MAX_CODE_LEN];
    strncpy(buf, code, MAX_CODE_LEN);
    char * start = strchr(buf,'(');
    char * end = strchr(buf,')');
    size_t len = end - start - 1;
    strncpy(code, start+1, len);
    code[len] = '\0'; 
}

word_pair parse_word_pair(char* code)
{
    word_pair result;
    char * start = strchr(code,'(');
    char * mid = strchr(code,',');
    char * end = strchr(code,')');
    size_t first_len = mid - start - 1;
    size_t second_len = end - mid - 1;
    strncpy(result.first, start+1, first_len);
    strncpy(result.second, mid+1, second_len);
    result.first[first_len] = '\0';
    result.second[second_len] = '\0';
    return result; 
}

// Elimination Tools

// Gives number of g and G in a word
int g_length(char* word)
{
    int g_len = 0;
    char* c = word;
    while (*c != '\0') {
      if (*c == 'g' || *c == 'G') {
        ++g_len;
      }
      ++c;
    }
    return g_len;
}

// Cusp area lower bound
inline const double areaLB(const XParams& nearer, char* where)
{
    // Area is |lox_sqrt|^2*|Im(lattice)|.
    XComplex lox_sqrt = nearer.loxodromic_sqrt;
    double lat_im     = nearer.lattice.im;
    if (lat_im  < 0) { // this should never happen
        fprintf(stderr, "Fatal: lattice imaginary part is negative at %s\n", where);
        exit(5);
    }
    // Apply Lemma 7.0 of GMT.
    double lox_re = (1-EPS)*(lox_sqrt.re*lox_sqrt.re);
    double lox_im = (1-EPS)*(lox_sqrt.im*lox_sqrt.im);
    double lox_norm = (1-EPS)*(lox_re + lox_im);
    return (1-EPS)*(lox_norm*lat_im);
}

// Our compact parameter space has the following bounds:
// 0. |lox_sqrt| >= 1 
// 1. 
//    a. Im(lox_sqrt) >= 0
//    b. Im(lattice) >= 0 
//    c. Im(parabolic) >= 0
//    d. Re(parabolic) >= 0
// 2. -1/2 <= Re(lattice) <= 1/2
// 3. |lattice| >= 1
// 4. Im(parabolic) <= Im(lattice)/2 
// 5. Re(parabolic) <= 1/2
// 6. |lox_sqrt^2| Im(L) <= g_max_area (area of fundamental paralleogram)
void verify_out_of_bounds(char* where, char bounds_code)
{
    Box box = build_box(where);
    switch(bounds_code) {
        case '0': {
            check(absUB(box.further.loxodromic_sqrt) < 1, where);
            break; } 
        case '1': {
            check(box.greater.loxodromic_sqrt.im < 0
             || box.greater.lattice.im < 0
             || box.greater.parabolic.im < 0
             || box.greater.parabolic.re < 0, where);
            break; } 
        case '2': {
            check(fabs(box.nearer.lattice.re) > 0.5, where);
            break; } 
        case '3': {
            check(absUB(box.further.lattice) < 1, where);
            break; } 
        case '4': {
            // Note: we can exclude the box if and only if the parabolic imag part is
            // bigger than half the lattice imag part over the WHOLE box
            // Multiplication by 0.5 is EXACT (if no underflow or overflow)
            check(box.nearer.parabolic.im > 0.5 * box.further.lattice.im, where);
            break; } 
        case '5': {
            check(box.nearer.parabolic.re > 0.5, where);
            break; } 
        case '6': {
            double area = areaLB(box.nearer, where);
            check(area > g_max_area, where);
            break;
        }
    }
}

// Check that the matrix is NOT of the forms
// 1 b  OR  -1  b
// 0 1       0 -1
// anywhere in the box
const int not_parabolic_at_inf(const SL2ACJ& x) {
    return absLB(x.c) > 0
        || ((absLB(x.a - 1) > 0 ||  absLB(x.d - 1) > 0) && (absLB(x.a + 1) > 0 || absLB(x.d + 1) > 0));
}

// Check that the matrix is NOT of the forms
// 1  OR  -1  0
// 0 1       0 -1
// anywhere in the box
const int not_identity(const SL2ACJ& x) {
    return absLB(x.b) > 0 || not_parabolic_at_inf(x);
}

// The infinity horoball has height t = 1/|loxodromic_sqrt|. An SL2C matrix
// a b
// c d
// Takes an infinity horoball of height t to a horoball of height 1/(t |c|^2)
// We want 1/(t |c|^2) > t. With t = 1/|loxodromic_sqrt|, this gives
// |c / loxodromic_sqrt| < 1.
const int large_horoball(const SL2ACJ& x, const ACJParams& p) {
    return absUB( x.c / p.loxodromic_sqrt ) < 1;
}

// Conditions checked:
//  1) word is not a parabolic fixing infinity anywhere in the box
//  2) word(infinity_horoball) intersects infinity_horoball
void verify_killed(char* where, char* word)
{
    Box box = build_box(where);
    SL2ACJ w = construct_word(box.cover, word);

    check(large_horoball(w, box.cover), where);
    check(not_parabolic_at_inf(w), where);
}

// Conditions checked:
//  1) word has g-length at most g_len 
//  2) word(infinity_horoball) intersects infinity_horoball
void verify_len(char* where, char* word, int g_len)
{
    Box box = build_box(where);
    SL2ACJ w = construct_word(box.cover, word); 

    check(g_length(word) <= g_len, where); 
    check(large_horoball(w, box.cover), where);
}

// Conditions checked:
//  1) the box is inside the variety neighborhood for given word 
void verify_variety(char* where, char* variety)
{
    Box box = build_box(where);
    SL2ACJ w = construct_word(box.cover, variety); 

    check((absUB(w.c) < 1) && (absUB(w.b) < 1 || absLB(w.c) > 0), where);
}
