/*
 * slopes.h
 *
 * Defines interface for slope computations in search.
 * Currently focused on slope distance
 *
 */

#include "slopes.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector>

using namespace std;

int gcd(int u, int v) {
    while (v != 0) {
        int r = u % v;
        u = v;
        v = r;
    }
    return u;
}

const int slope_dist(slope s1, slope s2) {
    int dist = s1.first * s2.second - s1.second * s2.first;
    if (dist < 0) { return -dist; }
    else { return dist; }
}

inline const double areaLB(const Params<XComplex>&nearer, char* where)
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

const slope short_slopes_dist_and_count(char* where) {
    Box box(where);
    slope error(-1,-1);
    vector< slope > short_slopes;
    Params<XComplex> nearer = box.nearer();
    Params<ACJ> cover = box.cover();
    ACJ S = cover.loxodromic_sqrt;
    ACJ L = cover.lattice;
    ACJ SL = S*L;
    fprintf(stderr, "Box %s has min area: %f nearer lat: %f + I %f lox: %f + I %f par: %f + I %f\n", where, areaLB(nearer, where), nearer.lattice.re, nearer.lattice.im, nearer.loxodromic_sqrt.re, nearer.loxodromic_sqrt.im, nearer.parabolic.re, nearer.parabolic.im);
    if (nearer.lattice.im <= 0) {
        fprintf(stderr, "Error: lattice with zero imaginary part in box %s\n", where);
        return error;
    }
    if (fabs(nearer.lattice.re) > 0.5) {
        fprintf(stderr, "Error: abs(lattice real part) > 1/2 in box %s\n", where);
        return error;
    }
    double b = 0;
    // fprintf(stderr, "absUB(S*6) is %f\n", absUB(S * 6)); 
    while (absLB((S*nearer.lattice.im)*b) <= 6) {
        double a = 0;  
        // fprintf(stderr, "b is %d\n", (int) b);
        while ((1-EPS)*(absLB((S*a)*(S*a)) + absLB((SL*b)*(SL*b))) <= (1+EPS)*(36 + absUB((S*a)*(S*b)))) {
            // fprintf(stderr, "a is %d and b is %d\n", (int) a, (int) b);
            for (int sgn = -1; sgn < 2; sgn += 2) {
                ACJ g = ((L * b) + (sgn * a))*S;
                // avoid long slopes, non-primatives, and overconting (-1,0) and (1,0) or (0,1)
                if ((a == 0 && b == 0) || (sgn == -1 && (a == 0 || b == 0)) ||
                    gcd((int) a, (int) b) > 1 || absLB(g) > 6) {
                    continue;
                } 
                slope p((int) sgn*a, (int) b);
                short_slopes.push_back(p);
            }
            a += 1;
        }
        b += 1;
        // fprintf(stderr, "absLB(b * nearer.lattice.im) is %f and  absUB(S*6) is %f\n", absLB(b * nearer.lattice.im), absUB(S * 6)); 
    }
    int num_slopes = short_slopes.size();
    int max_dist = 0;
    for ( auto it1 = short_slopes.begin(); it1 != short_slopes.end(); ++it1) {
        for ( auto it2 = it1; it2 != short_slopes.end(); ++it2) {
            if (it1 == it2) { continue; }
            max_dist = max(max_dist, slope_dist(*it1, *it2));
        }    
    }
    slope ans(max_dist, num_slopes);
    // box code, max_dist, num_short_slopes, the slopes
    fprintf(stdout, "%s,%d,%d", where, max_dist, num_slopes);
    for ( auto it = short_slopes.begin(); it != short_slopes.end(); ++it) {
        fprintf(stdout, ",{%d,%d}", it->first, it->second);
    }
    fprintf(stdout, "\n");
    return ans;
}

int main(int argc, char** argv)
{
    short_slopes_dist_and_count(argv[1]);	
}
/*
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

void walk(char* where, size_t depth, size_t* count_ptr)
{
    check(depth < MAX_DEPTH, where);
    *count_ptr += 1;
    char code[MAX_CODE_LEN];
    fgets(code, MAX_CODE_LEN, stdin);
    switch(code[0]) {
        case 'X': { 
            *count_ptr -= 1; // don't count branch nodes
            where[depth] = '0';
            where[depth + 1] = '\0';
            walk(where, depth + 1, count_ptr);
            where[depth] = '1';
            where[depth + 1] = '\0';
            walk(where, depth + 1, count_ptr);
            break; }
        case 'V': { // Line has format V(word) - box in variety nhd
            print_exceptional(where);
            break; }
        default: {
            return;
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
    initialize_roundoff();
    size_t count = 0;
    walk(where, depth, &count);
    if(!roundoff_ok()){
        printf(". underflow may have occurred\n");
        exit(1);
    }
    exit(0);
} */
