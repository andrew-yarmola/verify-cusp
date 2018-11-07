#include <stdio.h> 
#include <stdlib.h>
#include <string.h> 
#include "roundoff.h"
#include "Box.h"

#define MAX_DEPTH 200
#define MAX_CODE_LEN 100
#define MAX_AREA 5.24

bool check(bool inequalities, char* where)
{
    //printf("%d\n", inequalities);
    if (!inequalities) {       
        fprintf(stderr, "verify: fatal error at %s\n", where);
        return false;
//        exit(1);
    }
    return true;
}

inline const double areaLB(const Params<XComplex>&nearer)
{
    // Area is |lox_sqrt|^2*|Im(lattice)|.
    XComplex lox_sqrt = nearer.loxodromic_sqrt;
    double lat_im     = nearer.lattice.im;
    if (lat_im  < 0) { // this should never happen
      printf("%s\n", "ERROR : Have negative Im(L)!!!");
      lat_im = -lat_im;
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
//    b. Im(L) >= 0 
//    c. Im(P) >= 0
//    d. Re(P) >= 0
// 2. -1/2 <= Re(L) <= 1/2
// 3. |L| >= 1
// 4. Im(P) <= Im(L)/2 
// 5. Re(P) <= 1/2
// 6. |lox_sqrt^2| Im(L) <= MAX_AREA (area of fundamental paralleogram)
void verify_out_of_bounds(char* where, char bounds_code)
{
    Box box(where);
    Params<XComplex> nearer = box.nearer();
    Params<XComplex> further = box.further();
    Params<XComplex> greater = box.greater();
    switch(bounds_code) {
        case '0': {
            check(absUB(further.loxodromic_sqrt) < 1, where);
            break; } 
        case '1': {
            check(greater.loxodromic_sqrt.im < 0
             || greater.lattice.im < 0
             || greater.parabolic.im < 0
             || greater.parabolic.re < 0, where);
            break; } 
        case '2': {
            check(fabs(nearer.lattice.re) > 0.5, where);
            break; } 
        case '3': {
            check(absUB(further.lattice) < 1, where);
            break; } 
        case '4': {
            // Note: we can exclude the box if and only if the parabolic imag part is
            // bigger than half the lattice imag part over the WHOLE box
            // Multiplication by 0.5 is EXACT (if no underflow or overflow)
            check(nearer.parabolic.im > 0.5*further.lattice.im, where);
            break; } 
        case '5': {
            check(nearer.parabolic.re > 0.5, where);
            break; } 
        case '6': {
            double area = areaLB(nearer);
            check(area > MAX_AREA, where);
            break;
        }
    }
}

// Check that the matrix is NOT of the forms
// 1 b  OR  -1  b
// 0 1       0 -1
// anywhere in the box
const int not_parabolic_at_inf(const SL2ACJ&x) {
    return absLB(x.c) > 0
        || ((absLB(x.a - 1.) > 0 ||  absLB(x.d - 1.) > 0) && (absLB(x.a + 1.) > 0 || absLB(x.d + 1.) > 0));
}

// Check that the matrix is NOT of the forms
// 1 0  OR  -1  0
// 0 1       0 -1
// anywhere in the box
const int not_identity(const SL2ACJ&x) {
    return absLB(x.b) > 0 || not_parabolic_at_inf(x);
}

// The infinity horoball has height t = 1/|loxodromic_sqrt|. An SL2C matrix
// a b
// c d
// Takes an infinity horoball of height t to a horoball of height 1/(t |c|^2)
// We want 1/(t |c|^2) > t. With t = 1/|loxodromic_sqrt|, this gives
// |c / loxodromic_sqrt| < 1.
const int large_horoball(const SL2ACJ&x, const Params<ACJ>&p) {
    return absUB( x.c / p.loxodromic_sqrt ) < 1;
}

void debug_info(char* where, char* word)
{
    Box box(where);
    Params<ACJ> params = box.cover();
    Params<XComplex> nearer = box.nearer();
    SL2ACJ w = construct_word(params, word);
    fprintf(stderr, "Params: %s\n", where);
    fprintf(stderr, "L: %f + I %f size %f\n", params.lattice.f.re, params.lattice.f.im, params.lattice.size);
    fprintf(stderr, "S: %f + I %f size %f\n", params.loxodromic_sqrt.f.re, params.loxodromic_sqrt.f.im, params.loxodromic_sqrt.size);
    fprintf(stderr, "P: %f + I %f size %f\n", params.parabolic.f.re, params.parabolic.f.im, params.parabolic.size);
    XComplex a = w.a.f;
    XComplex b = w.b.f;
    XComplex c = w.c.f;
    XComplex d = w.d.f;
    ACJ det = w.a * w.d - w.b * w.c;
    fprintf(stderr, "Word: %s\n", word);
    fprintf(stderr, "At the center is has coords\n");
    fprintf(stderr, "a: %f + I %f\n", a.re, a.im);
    fprintf(stderr, "b: %f + I %f\n", b.re, b.im);
    fprintf(stderr, "c: %f + I %f\n", c.re, c.im);
    fprintf(stderr, "d: %f + I %f\n", d.re, d.im);
    fprintf(stderr, "det : %f + I %f\n", det.f.re, det.f.im); 
    SL2ACJ G = construct_word(params, "G");
    a = G.a.f;
    b = G.b.f;
    c = G.c.f;
    d = G.d.f;
    det = G.a * G.d - G.b * G.c;
    fprintf(stderr, "Word: G\n");
    fprintf(stderr, "At the center is has coords\n");
    fprintf(stderr, "a: %f + I %f\n", a.re, a.im);
    fprintf(stderr, "b: %f + I %f\n", b.re, b.im);
    fprintf(stderr, "c: %f + I %f\n", c.re, c.im);
    fprintf(stderr, "d: %f + I %f\n", d.re, d.im);
    fprintf(stderr, "det : %f + I %f\n", det.f.re, det.f.im); 

    fprintf(stderr, "absLB(c) : %f \n absLB(a - 1) : %f\n absLB(d - 1) : %f\n absLB(a+1) : %f\n absLB(d+1) : %f\n",
                            absLB(w.c), absLB(w.a - 1.),  absLB(w.d - 1.), absLB(w.a + 1.), absLB(w.d + 1.));
    fprintf(stderr, "horo_ratio : %f\n", absUB( w.c / params.loxodromic_sqrt ) );
    
    fprintf(stderr, "det_norm LB : %f, UB : %f\n", absLB( w.a * w.d - w.b * w.c), absUB( w.a * w.d - w.b * w.c) );
    fprintf(stderr, "G det_norm LB : %f, UB : %f\n", absLB( G.a * G.d - G.b * G.c), absUB( G.a * G.d - G.b * G.c) );

    fprintf(stderr, "L : %f\n S : %f\n P : %f\n", absLB(params.lattice - ACJ(XComplex(0.5, 0.8660254037844386467637231707529361834714))),
                                                  absLB(params.loxodromic_sqrt - ACJ(XComplex(2., 0.))),
                                                  absLB(params.parabolic - ACJ(XComplex(0.25, 0.4330127018922193233818615853764680917357))));
    fprintf(stderr, "Area LB : %f\n", areaLB(nearer));
    
    double one = 1; // Exact
    ACJ T = (absUB(w.d - one) < 2 || absUB(w.a - one) < 2) ? w.b : -w.b;
    ACJ L = params.lattice;

    ACJ d1 = T / (L + one);
    ACJ d2 = d1 - one; // uses fewer operations
    ACJ d3 = (T - one) / (L - one);
    ACJ d4 = d3 - one; // better error estimate

    fprintf(stderr, "absUB(d1) = %f\n absUB(d2) = %f\n absUB(d3) = %f\n absUB(d4) = %f\n", absUB(d1), absUB(d2), absUB(d3), absUB(d4));
}

// Conditions checked:
//  1) word is not a parabolic fixing infinity anywhere in the box
//  2) word(infinity_horoball) intersects infinity_horoball
void verify_killed(char* where, char* word)
{
    Box box(where);
    Params<ACJ> params = box.cover();
    SL2ACJ w = construct_word(params, word);

    check(not_parabolic_at_inf(w), where);
    check(large_horoball(w, params), where);
}

// Conditions checked:
//  1) the box is inside the variety neighborhood for giver word 
void verify_variety(char* where, char* variety)
{
    Box box(where);
    Params<ACJ> params = box.cover();
    SL2ACJ w = construct_word(params, variety); 

    check((absUB(w.c) < 1) && (absUB(w.b) < 1 || absLB(w.c) > 0), where);
}

// Conditions checked:
//  1) word(infinity_horoball) intersects infinity_horoball
//  2) check word is one that cannot be a parabolic fixing infinity
void verify_parabolic_always_impossible(char* where, char* word)
{
    Box box(where);
    Params<ACJ> params = box.cover();
    SL2ACJ w = construct_word(params, word);

    check(large_horoball(w, params), where);

    // TODO Finish -- load file of impossible parabolics
    // fprintf(stderr, "verify: no implementation of checking impossible parabolic contradiction at %s\n", where);
}

// Conditions checked:
//  1) word(infinity_horoball) intersects infinity_horoball
//  2) checks that a parabolic translate of the word is a power of a subword that cannot be parabolic
void verify_parabolic_impossible(char* where, char* word, char* subword)
{
    Box box(where);
    Params<ACJ> params = box.cover();
    SL2ACJ w = construct_word(params, word);
    SL2ACJ v = construct_word(params, subword);

    check(large_horoball(w, params), where);
    check(not_parabolic_at_inf(v), where);
    // TODO Finish -- load file of impossible parabolics
    // fprintf(stderr, "verify: no implementation of checking impossible identity contradiction at %s\n", where);
}

// Conditions checked:
//  1) word(infinity_horoball) intersects infinity_horoball
//  2) if word is a parabolic fixinig infinity, it must be the idenity
//  3) word on the list of those that cannot be the idenity
void verify_identity_always_impossible(char* where, char* word)
{
    Box box(where);
    Params<ACJ> params = box.cover();
    SL2ACJ w = construct_word(params, word);

//    debug_info(where, word);
    check(large_horoball(w, params), where);
//    fprintf(stderr,"try horo - not ident\n"); 
//    fprintf(stderr,"pass horo - not ident\n");
//    fprintf(stderr, "absUB(w.b) = %f\n", absUB(w.b));
    check(absUB(w.b) < 1, where);
    // TODO Finish -- load file of impossible parabolics
    // fprintf(stderr, "verify: no implementation of checking impossible identity contradiction at %s\n", where);
}

// Conditions checked:
//  1) word(infinity_horoball) intersects infinity_horoball
//  2) if word is a parabolic fixinig infinity, it must be the idenity
//  3) word is not the idenity
void verify_indiscrete_lattice_simple(char* where, char* word)
{
    Box box(where);
    Params<ACJ> params = box.cover();
    SL2ACJ w = construct_word(params, word);

    check(large_horoball(w, params), where);
    check(absUB(w.b) < 1, where);
    check(not_identity(w), where);
//
//    bool h = check(large_horoball(w, params), where);
//    bool b = check(absUB(w.b) < 1, where);
//    bool id = check(not_identity(w), where);
//
//    if (!h || !b || ! id) {
//        fprintf(stderr,"Q FAILURE+++++++++++++++\n");
//        debug_info(where, word);
//    }
//    if (!h) {
//        fprintf(stderr,"fail horo - indiscrete_lattice_simple\n");
//    }
//    if (!b) {
//        fprintf(stderr,"fail b bound - indiscrete_lattice_simple\n");
//        fprintf(stderr, "absUB(w.b) = %f\n", absUB(w.b));
//    }
//    if (!id) {
//        fprintf(stderr,"fail not ident - indiscrete_lattice_simple\n");
//        fprintf(stderr," absLB(b) = %f\n absLB(c) = %f\n absLB(a-1) = %f\n absLB(d-1) = %f\n absLB(a+1) = %f\n absLB(d+1) = %f\n",
//                         absLB(w.b), absLB(w.c), absLB(w.a - 1.), absLB(w.d - 1.), absLB(w.a + 1.), absLB(w.d + 1.));
//    } 
}

// Conditions checked:
//  1) word(infinity_horoball) intersects infinity_horoball
//  2) at the points where the word is parabolic, it is not on the lattice
void verify_indiscrete_lattice(char* where, char* word)
{
    Box box(where);
    Params<ACJ> params = box.cover();
    SL2ACJ w = construct_word(params, word);
    double one = 1; // Exact

//    debug_info(where, word);

    check(large_horoball(w, params), where);
    
    // For all parabolic points in the box, we want verify
    // that none of them are lattice points. At such a point, the parabolic
    // translation will be given as +/- w.b.
    //
    // We check the box is small enough to determine the sign.
    check(absUB(w.d - one) < 2 || absUB(w.d + one) < 2 || 
          absUB(w.a - one) < 2 || absUB(w.a + one) < 2, where);
    
    ACJ T = (absUB(w.d - one) < 2 || absUB(w.a - one) < 2) ? w.b : -w.b;
    ACJ L = params.lattice;

    // There are now 4 equations to check corresponding to the intersection
    // of 4 circles :
    // |translation - 0          | < |1 + lattice|
    // |translation - (1+lattice)| < |1 + lattice|
    // |translation - 1          | < |1 - lattice|
    // |translation - lattice    | < |1 - lattice|
    // These inequailties show that transltion is not on the lattice (assuming
    // parameterd space constraitns). See proof in text.
    // 
    // To make the computation efficient, rearange and take absolute values at the end.

    ACJ d1 = T / (L + one);
    ACJ d2 = d1 - one; // uses fewer operations
    ACJ d3 = (T - one) / (L - one);
    ACJ d4 = d3 - one; // better error estimate

//    fprintf(stderr, "%f, %f, %f, %f\n", absLB(L + one) - absUB(T),
//                                        absLB(L + one) - absUB(T - L - one),
//                                        absLB(L - one) - absUB(T - one),
//                                        absLB(L - one) - absUB(T - L));

    check(absUB(d1) < 1 && absUB(d2) < 1 && absUB(d3) < 1 && absUB(d4) < 1, where);
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

void verify(char* where, size_t depth, size_t* count_ptr)
{
    check(depth < MAX_DEPTH, where);
    // TODO: Make a conditional list file and update the tree with conditions
    char code[MAX_CODE_LEN];
    fgets(code,MAX_CODE_LEN, stdin);
    //printf("%s CODE %s\n", where, code);
    switch(code[0]) {
        case 'X': { 
            where[depth] = '0';
            where[depth+1] = '\0';
            verify(where, depth+1, count_ptr);
            where[depth] = '1';
            where[depth+1] = '\0';
            verify(where, depth+1, count_ptr);
            break; }
        *count_ptr += 1;
        case '0': 
        case '1': 
        case '2': 
        case '3': 
        case '4': 
        case '5': 
        case '6': {
            verify_out_of_bounds(where, code[0]);
            break; }
        case 'K': { // Line has format  K(word) - killer word
            parse_word(code);
            verify_killed(where, code);
            break; }
        case 'V': { // Line has format V(word) - box in variety nhd
            parse_word(code);
            verify_variety(where, code);
            break; }
        case 'P': { // Line has format P(word) - word cannot be parabolic
            parse_word(code);
            verify_parabolic_always_impossible(where, code);
            break; } 
        case 'I': { // Line has format I(word) - impossible identity 
            parse_word(code);
            verify_identity_always_impossible(where, code);
            break; }
        case 'Q': { // Line has format Q(word) - failed quasi-relator 
            parse_word(code);
            verify_indiscrete_lattice_simple(where, code);
            break; }
        case 'L': { // Line has format L(word) - all parabolics indiscrete
            parse_word(code);
            verify_indiscrete_lattice(where, code);
            break; } 
        case 'E': { // Line has format E(word, subword) - a word that cannot be parabolic has parabolic power 
            parse_word(code);
            char * comma = strchr(code,',');
            check(comma != NULL, where);
            char * subword = comma + 1;
            comma[0] = '\0'; 
            char * word = code;
            verify_parabolic_impossible(where, word, subword);
            break; }
        case 'H' : {
            printf("HOLE - skipping\n");
            break;
        } 
        default: {
            check(false, where);
        }
    }
}

int main(int argc, char**argv)
{
    if(argc != 2) {
        fprintf(stderr,"Usage: %s position < data\n", argv[0]);
        exit(1);
    }
    char where[MAX_DEPTH];
    size_t depth = 0;
    while (argv[1][depth] != '\0') {
        if (argv[1][depth] != '0' && argv[1][depth] != '1'){
            fprintf(stderr,"bad position %s\n",argv[1]);
            exit(2);
        }
        where[depth] = argv[1][depth];
        depth++;
    }
    where[depth] = '\0';

    printf("verified %s - {\n", where);
    initialize_roundoff();
    size_t count = 0;
    verify(where, depth, &count);
    if(!roundoff_ok()){
        printf(". underflow may have occurred\n");
        exit(1);
    }
    printf("Verified %lu nodes\n", count);
    printf("}.\n");
    exit(0);
}
