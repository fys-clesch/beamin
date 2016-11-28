#include "include_me.h"
#include "msg.h"
#include "allocer.h"

double *alloc_matrix(int row, int col, int add, double init)
{
    int i;
    double *m;
    m = (double *)malloc((row * (col + add)) * sizeof(double));
    if(NULL == m)
    {
        error_msg("memory allocation failed in alloc_matrix. good bye.", __FILE__, __LINE__);
        free(m);
        exit(EXIT_FAILURE);
    }
    for(i = 0; i < (row * (col + add)); i++) m[i] = init;
    return m;
}

double *alloc_vector(int row, double init)
{
    int i;
    double *m;
    m = (double *)malloc(row * sizeof(double));
    if(NULL == m)
    {
        error_msg("memory allocation failed in alloc_vector. good bye.", __FILE__, __LINE__);
        free(m);
        exit(EXIT_FAILURE);
    }
    for(i = 0; i < row; i++) m[i] = init;
    return m;
}

pset *alloc_pset(int count, double init)
{
    int i, j;
    pset *m;
    m = (pset *)malloc((count) * sizeof(pset));
    if(NULL == m)
    {
        error_msg("memory allocation failed in alloc_pset. good bye.", __FILE__, __LINE__);
        free(m);
        exit(EXIT_FAILURE);
    }
    for(i = 0; i < count; i++)
    {
        for(j = 0; j < 20; j++) m[i].set[j] = init;
        m[i].eval = init;
    }
    return m;
}

target *alloc_target(int count)
{
    int i;
    target *m;
    m = (target *)malloc((count) * sizeof(target));
    if(NULL == m)
    {
        error_msg("memory allocation failed in alloc_target. good bye.", __FILE__, __LINE__);
        free(m);
        exit(EXIT_FAILURE);
    }
    for(i = 0; i < count; i++)
    {
        m[i].val = 0.;
        m[i].wgt = 1.;
    }
    return m;
}

cplx *alloc_cplx(int count)
{
    int i;
    cplx *m;
    m = (cplx *)malloc((count) * sizeof(cplx));
    if(NULL == m)
    {
        error_msg("memory allocation failed in alloc_cplx. good bye.", __FILE__, __LINE__);
        free(m);
        exit(EXIT_FAILURE);
    }
    for(i = 0; i < count; i++)
        m[i].R = m[i].I = 0.;
    return m;
}

parameter *alloc_parameter(int row)
{
    int i;
    parameter *m;
    m = (parameter *)malloc(row * sizeof(parameter));
    if(NULL == m)
    {
        error_msg("memory allocation failed in alloc_parameter. good bye.", __FILE__, __LINE__);
        free(m);
        exit(EXIT_FAILURE);
    }
    for(i = 0; i < row; i++)
    {
        strncpy(m[i].name, "", 12);
        strncpy(m[i].unit, "", 12);
        m[i].init = m[i].min = m[i].max = m[i].val = 0.;
        m[i].log = m[i].fit = 0;
    }
    return m;
}
