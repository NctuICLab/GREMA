#ifndef stdafx_H
#define stdafx_H
//#pragma once
#include <iostream>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <fstream>
#include <string.h>
#include <string>
#include <vector>
#include <algorithm>
#include <unistd.h>
#define MAXDOUBLE		pow(2.0,30.0)
#define MINDOUBL        pow(2.0,-50)
#define PI				3.14159265358979
#define DELTA_S			0.03f
#define PAUSE printf("Press any key to continue...\n"); fgetc(stdin);
#define neg_zero        -9999
#define XMIN			0.0001
#define XMAX			100000
#define BUF_SIZE		2048
#define NAN_PENALTY		1e30	//fitness given to a model whose objective function is NaN (= worst)
#endif
