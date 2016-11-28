#include "include_me.h"
#include "msg.h"
#include "optic.h"
#include "allocer.h"
#include "ff.h"
#include "nonlin_fit.h"

static qbeam q, qp;
static double wtar;

double ff_trip_simplex(double *x)
{
    double ls1, ls2, M1[4], M2[4];

    propagation(M1, x[0]);
    thin_lens(M2, 400.);
    M_mult_ip(M2, M1);

    propagation(M1, x[1]);
    M_mult_ip(M1, M2);

    thin_lens(M2, 400.);
    M_mult_ip(M2, M1);

    propagation(M1, x[2]);
    M_mult_ip(M1, M2);

    thin_lens(M2, 200.);
    M_mult_ip(M2, M1);

    propagation(M1, 2000. - (x[0] + x[1] + x[2])); //Overall length: 2000. mm
    M_mult_ip(M1, M2);

    qp = qtrans(M1, &q);
    ls1 = get_w0(&qp) - wtar;
    ls2 = get_z(&qp) - 0.;

    if(x[0] <= 1. || x[1] <= 1. || x[2] <= 1.) return 1000.*(ls1 * ls1 + ls2 * ls2);
    else return ls1 * ls1 + ls2 * ls2;
}

void manff_simplex(const double lam, const double w0, const double wf)
{
    wtar = wf;
    int i, Ndim = 3, Mdim = 4, steps = 0;
    pset *testO = alloc_pset(Mdim, 0.);
    double *init = alloc_vector(Ndim, 0.);
    setq(&q, 0., rayleighr(w0, lam), lam);

    init[0] = testO[0].set[0] = 500.;
    init[1] = testO[0].set[1] = 200.;
    init[2] = testO[0].set[2] = 920.;
    testO[0].eval = ff_trip_simplex(init);
    init[0] = testO[1].set[0] = 550.;
    init[1] = testO[1].set[1] = 200.;
    init[2] = testO[1].set[2] = 921.;
    testO[1].eval = ff_trip_simplex(init);
    init[0] = testO[2].set[0] = 550.;
    init[1] = testO[2].set[1] = 201.;
    init[2] = testO[2].set[2] = 920.;
    testO[2].eval = ff_trip_simplex(init);
    init[0] = testO[3].set[0] = 551.;
    init[1] = testO[3].set[1] = 200.;
    init[2] = testO[3].set[2] = 920.;
    testO[3].eval = ff_trip_simplex(init);

    small_amoeba(testO, Ndim, 1e-14, ff_trip_simplex, &steps);

    for(i = 0; i < Mdim; i++)
        printf("%12.10f  %12.10f  %12.10f  %12.10g\n", testO[i].set[0], testO[i].set[1], testO[i].set[2], testO[i].eval);

    free(testO);
    free(init);
}
