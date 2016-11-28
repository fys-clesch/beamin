#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <float.h>
#include <limits.h>
#define M_PI 3.1415926535897932
#define M_PIh 1.5707963267948966
#define M_EULER 2.7182818284590452

typedef struct
{
    double set[20];
    double eval;
} pset;

typedef struct
{
    double R, I;
} cplx;

typedef struct
{
    double set[10];
    double eval;
} p_set;

typedef struct
{
    char name[12];
    char unit[12];
    double init;
    double min, max;
    int log;
    int fit;
    double val;
} parameter;

typedef struct
{
    double val, wgt;
} target;

typedef struct
{
    double R, I;
    double lam;
} qbeam;
