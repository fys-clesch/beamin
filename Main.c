/*
 clemens schaefermeier
 clemens@fh-muenster.de
*/

#include "include_me.h"
#include "msg.h"
#include "optic.h"
#include "ff.h"

int main()
{
//Find best parameters to yield a beam waist of 30 um


    double lambda = 632.8e-6;
    double w0 = 1e-1;
    double target = 3e-2;

    manff_simplex(lambda, w0, target);

    /*
    The output should be read as follows

    MESSAGE:
    NMAX exceeded in small_amoeba
    528.3232857859  222.2128733950  961.1164778513  1.412123305e-026
    528.3232857859  222.2128733950  961.1164778513  1.36018537e-026
    528.3232857859  222.2128733950  961.1164778513  1.36018537e-026
    528.3232857859  222.2128733950  961.1164778513  1.412123305e-026
    */

//Ray propagation through optical elements

    /*
    double lambda=632.8e-6;
    double w0=1e-1;
    double M1[4],M2[4];
    qbeam q,qp;
    setq(&q,0.,rayleighr(w0,lambda),lambda);

    q_out(&q,"in");

    propagation(M1,400.);
    thin_lens(M2,258.);
    M_mult_ip(M2,M1);

    propagation(M1,1400.);
    M_mult_ip(M1,M2);

    thin_lens(M2,140.4);
    M_mult_ip(M2,M1);

    propagation(M1,172.1);
    M_mult_ip(M1,M2);

    qp=qtrans(M1,&q);

    q_out(&qp,"out");
    */
    /*
    The output should be read as follows

    q.in    :      z[e-3m]     w0[e-6m]
                         0          100
                   R[e-3m]     wz[e-6m]
                       inf          100
    q.out   :      z[e-3m]     w0[e-6m]
                -0.0721956      40.0079
                   R[e-3m]     wz[e-6m]
                  -874.728      40.0095
    */

//A linear resonator, left mirror is plane, distance to right mirror is 400 mm, curvature of right mirror is 2 m

    /*
    double M1[4],M2[4];
    qbeam q;
    cplx c;
    q.lam=1064e-6;

    propagation(M1,400.);
    reflection(M2,2000.);
    M_mult_ip(M2,M1);

    propagation(M1,400.);
    M_mult_ip(M1,M2);

    c=eigenmod_res(M1);
    q.R=c.R;
    q.I=c.I;

    q_out(&q,"res");
    */
    /*
    The output should be read as follows

    q.res   :      z[e-3m]     w0[e-6m]
                         0      520.524
                   R[e-3m]     wz[e-6m]
                       inf      520.524
    */
    return 0;
}
