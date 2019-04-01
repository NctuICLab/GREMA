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
#endif
