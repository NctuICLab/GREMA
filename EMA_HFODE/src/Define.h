#ifndef DEFINE_H
#define DEFINE_H
#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <fstream>
#include <cstring>
#include <cmath>
#include <climits>
double const MAXDOUBLE = pow(2.0,30.0);
double const MINDOUBL =  pow(2.0,-50);
double const PI = 3.14159265358979;
double const DELTA_S = 0.03f;
#define PAUSE printf("Press any key to continue...\n"); fgetc(stdin);
int const neg_zero = -9999;
double const XMIN = 0.0001;
int const XMAX	= 100000;
int const BUF_SIZE = 2048;
// Default values of IGA
double const DEFAULT_PC	= 0.8;
double const DEFAULT_PM1 = 0.02;
double const DEFAULT_PM2 = 0.05;
double const DEFAULT_PS	= 0.2;
double const DEFAULT_SA = 0.997;
// Maximum number of regulator slots in one chromosome.  The slot mask is stored as a
// bit set inside one double gene (exactly representable up to 53 bits) and is handled
// as unsigned long long, so at most 52 slots are supported.
int const MAX_SLOTS = 52;
// Fitness values are always kept finite: NaN/inf/overflow (e.g. from pow() of a
// degenerate parameter set) are mapped to this "worst" value so that every
// comparison and sort in the GA is well defined.
double const WORST_FITNESS = 1e30;
// NaN-safe view of a fitness value for comparisons: NaN counts as +infinity (worst).
inline double SafeFit(double v){ return (v == v) ? v : HUGE_VAL; }
#endif
