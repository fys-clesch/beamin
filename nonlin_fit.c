#include "include_me.h"
#include "allocer.h"
#include "msg.h"
#include "nonlin_fit.h"

//Adopted from Numerical Recipes, 2nd Ed.

#define TINY 1e-15 //For denumerator healing
#define NMAX 5000

//Calculate sum of each row of full_set: psum stores the sums of each of the "mdim-times" ndim parameters
#define GET_PSUM	for (j=0;j<ndim;j++)\
			{\
			for(sum=0.,i=0;i<mdim;i++) sum += p[i].set[j];\
			psum[j]=sum;\
			}
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

double small_amotry(pset *p, double *psum, int ndim, double (*func)(double *), int ihi, double fac)
{
    /*
     Extrapolates by "fac" through the face of the simplex across from the highest point,
     tries it, and replaces the high point if the new point is better.
     See explanation for parameters in function amoeba
    */

    int j;
    double fac1, fac2, res_try, ptry[ndim];

    fac1 = (1. - fac) / ndim;
    fac2 = fac1 - fac;

    for(j = 0; j < ndim; j++) ptry[j] = psum[j] * fac1 - p[ihi].set[j] * fac2; //Generate new pset set

    res_try = (*func)(ptry); //Evaluate the function at the trial point

    if(res_try < p[ihi].eval) //New one is better than the worst
    {
        p[ihi].eval = res_try;
        for(j = 0; j < ndim; j++)
        {
            psum[j] += ptry[j] - p[ihi].set[j]; //Correct psum
            p[ihi].set[j] = ptry[j]; //Set in the new parameters
        }
    }
    return res_try; //return func eval..
}

void small_amoeba(pset *p, int ndim, double ftol, double (*func)(double *), int *steps)
{
    /*
     Minimizes func
     x[ndim]: pset vector
     Matrix full_set: Starting matrix, formated like [ndim+1][ndim]
     res_vec[ndim+1]: Function evaluation at full_set
     steps: number of function evaluations
     ftol: Fractional convergence tolerance to be achieved
    */

    int i, j, ihi, ilo, sec_ihi;
    int mdim = ndim + 1;
    double rtol, res_try, res_save, *psum;
    double sum, swap;

    psum = alloc_vector(ndim, 0);

    *steps = 0;
    GET_PSUM

    for(;;)
    {
        //First, test were the best (lowest) and where the worst (highest eval.) is stored...
        ilo = 0;
        ihi = (p[0].eval) > (p[1].eval) ? (sec_ihi = 1, 0) : (sec_ihi = 0, 1); //If [0]>[1], set ihi=0 and sec_ihi=1, otherwise opposite

        for(i = 0; i < mdim; i++)
        {
            if(p[i].eval <= p[ilo].eval) ilo = i; //New minimum in function eval
            if(p[i].eval > p[ihi].eval)
            {
                sec_ihi = ihi; //Set second worst
                ihi = i; //Set worst
            }
            else if(p[i].eval > p[sec_ihi].eval && i != ihi) sec_ihi = i;
        }

        //"Contrast" between high and lo. Sth. like the relative size of the simplex.
        rtol = 2.*fabs(p[ihi].eval - p[ilo].eval) / (fabs(p[ihi].eval) + fabs(p[ilo].eval) + TINY);

        if(rtol < ftol) //If this is true, break
        {
            SWAP(p[0].eval, p[ilo].eval) //Store best one in slot 0
            for(i = 0; i < ndim; i++)
                SWAP(p[0].set[i], p[ilo].set[i])
                fprintf(stdout, "steps: %i\n", *steps);
            break; //Here's the end
        }
        if(*steps >= NMAX)
        {
            msg("NMAX exceeded in small_amoeba");
            break;
        }

        res_try = small_amotry(p, psum, ndim, func, ihi, -1.); //Reflection
        (*steps)++;

        if(res_try <= p[ilo].eval) //better than the best eval?
        {
            res_try = small_amotry(p, psum, ndim, func, ihi, 2.); //Expand further
            (*steps)++;
        }
        else if(res_try >= p[sec_ihi].eval) //worse than second worst?
        {
            res_save = p[ihi].eval;
            res_try = small_amotry(p, psum, ndim, func, ihi, .5); //Shrink
            (*steps)++;
            if(res_try >= res_save) //worse than worst? //if(res_try >= res_save)
            {
                for(i = 0; i < mdim; i++)
                    if(i != ilo)
                    {
                        //If so, then shrink the whole simplex around ilo. psum just for temp save.
                        for(j = 0; j < ndim; j++)
                            p[i].set[j] = psum[j] = .5 * (p[i].set[j] + p[ilo].set[j]);
                        p[i].eval = (*func)(psum); //func eval
                    }
                *steps += ndim; //Record all func evals
                GET_PSUM //Get psum
            }
        }
    }
    free(psum);
}

#undef GET_PSUM
#undef NMAX
#undef TINY
#undef SWAP
