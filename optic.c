#include "include_me.h"
#include "allocer.h"
#include "msg.h"
#include "optic.h"

void thin_lens(double *M, const double f)
{
    M[0] = 1.;
    M[1] = 0.;
    M[2] = -1. / f;
    M[3] = 1.;
}

void thick_lens(double *M, const double r1, const double r2, const double d, const double n)
{
    double ni = 1. / n, dr1, r1i, r2i;
    if(r1 != 0.) r1i = 1. / r1;
    else r1i = 0.;
    dr1 = d * r1i;
    if(r2 != 0.) r2i = 1. / r2;
    else r2i = 0.;
    M[0] = 1. + dr1 * (ni - 1.);
    M[1] = d * ni;
    M[2] = (1. - n) * (r1i - r2i) + dr1 * r2i * (2. - ni - n);
    M[3] = 1. - d * r2i * (ni - 1.);
}

void thick_p_ck_lens(double *M, const double r2, const double d, const double n)
{
    double ni = 1. / n;
    M[0] = 1.;
    M[1] = d * ni;
    M[2] = -(1. - n) / r2;
    M[3] = 1. - d / r2 * (ni - 1.);
}

void propagation(double *M, const double d)
{
    M[0] = 1.;
    M[1] = d;
    M[2] = 0.;
    M[3] = 1.;
}

void refraction(double *M, const double n1, const double n2, const double r)
{
    M[0] = 1.;
    M[1] = 0.;
    if(r != 0.) M[2] = (n1 - n2) / (n2 * r);
    else M[2] = 0.;
    M[3] = n1 / n2;
}

void reflection(double *M, const double r)
{
    M[0] = 1.;
    M[1] = 0.;
    M[2] = -2. / r;
    M[3] = 1.;
}

void grin(double *M, const double d, const double n0, const double nmin)
{
    double t1, t2, t3;
    t1 = sqrt(nmin / n0);
    t2 = t1 * d;
    t3 = sin(t2);
    M[0] = cos(t2);
    M[1] = t3 / t1;
    M[2] = -t1 * t3;
    M[3] = M[0];
}

void M_mult_ip(double *Mr, const double *M)
{
    double t1, t2;
    t1 = Mr[0] * M[0] + Mr[1] * M[2];
    t2 = Mr[0] * M[1] + Mr[1] * M[3];
    Mr[0] = t1;
    Mr[1] = t2;
    t1 = Mr[2] * M[0] + Mr[3] * M[2];
    t2 = Mr[2] * M[1] + Mr[3] * M[3];
    Mr[2] = t1;
    Mr[3] = t2;
}

qbeam qtrans(const double *M, const qbeam *q)
{
    qbeam qp;
    qp.lam = (*q).lam;
    double t1, t2, t3, t4, t5, denom;

    t1 = M[2] * (*q).R + M[3];
    t2 = t1 * t1;
    t3 = M[0] * (*q).R + M[1];
    t4 = M[2] * M[2];
    t5 = (*q).I * (*q).I;

    if((denom = t2 + t4 * t5) == 0.)
    {
        error_msg("singularity in qtrans. healing.", __FILE__, __LINE__);
        denom = 1e-15;
    }

    qp.R = t1 * t3 + t5 * M[0] * M[2];
    qp.R /= denom;
    qp.I = (*q).I * (M[0] * t1 - M[2] * t3);
    qp.I /= denom;
    return qp;
}

double wz(const double w0, const double z, const double zr)
{
    double t = z / zr;
    t *= t;
    return w0 * sqrt(1. + t);
}

double rayleighr(const double w0, const double lam)
{
    return M_PI * w0 * w0 / lam;
}

double stab_parm(const double *M)
{
    return ((M[0] + M[3]) / 2.);
}

cplx *eigenval_res(const double *M)
{
    cplx *c = alloc_cplx(2);
    double g = stab_parm(M), t;
    t = g * g - 1.;
    if(t < 0.)
    {
        t = sqrt(-t);
        c[0].R = c[1].R = g;
        c[0].I = t;
        c[1].I = -t;
    }
    else if(t > 0.)
    {
        t = sqrt(t);
        c[0].I = c[1].I = 0.;
        c[0].R = g + t;
        c[1].R = g - t;
    }
    else
    {
        c[0].I = c[1].I = 0.;
        c[0].R = c[1].R = g;
    }
    return c;
}

cplx eigenmod_res(const double *M)
{
    cplx q;
    double t1, t2, t3;
    t1 = M[3] - M[0];
    if(t1 != 0.) t3 = t1 * t1;
    else
    {
        q.R = 0.;
        if((t3 = M[1] / M[2]) >= 0.)
        {
            fprintf(stdout, "\nfound no eigenmode solution");
            q.I = 0.;
        }
        else q.I = sqrt(-t3);
        return q;
    }
    t2 = 2.*M[2];
    q.R = t1 / t2;
    t2 *= t2;
    if((t3 = (t3 + 4.*M[1] * M[2]) / t2) >= 0.)
    {
        fprintf(stdout, "\nfound no eigenmode solution");
        q.R = q.I = 0.;
    }
    else q.I = sqrt(-t3);
    return q;
}

void setq(qbeam *q, const double z, const double zr, const double lam)
{
    (*q).R = z;
    (*q).I = zr;
    (*q).lam = lam;
}

double get_w0(const qbeam *q)
{
    return sqrt((*q).I * (*q).lam / M_PI);
}

double get_z(const qbeam *q)
{
    return (*q).R;
}

double get_zr(const qbeam *q)
{
    return (*q).I;
}

double get_beamrad(const qbeam *q)
{
    double r = (*q).R;
    if(r != 0.) r += ((*q).I * (*q).I / r);
    else return DBL_MAX;
    return r;
}

double get_wz(const qbeam *q)
{
    double wz = (*q).I;
    wz += ((*q).R * (*q).R / (*q).I);
    wz *= ((*q).lam / M_PI);
    wz = sqrt(wz);
    return wz;
}

double get_div(const qbeam *q)
{
    return atan(sqrt((*q).lam / ((*q).I * M_PI)));
}

void q_out(const qbeam *q, const char *s)
{
    double r = get_beamrad(q);
    (r == DBL_MAX) ?
    fprintf(stdout, "q.%-6.6s: %12s %12s\n%9s %12g %12g\n"
            "%9s %12s %12s\n%9s %12s %12g\n"
            , s, "z[e-3m]", "w0[e-6m]", "", get_z(q), get_w0(q) * 1e3
            , "", "R[e-3m]", "wz[e-6m]", "", "inf", get_wz(q) * 1e3) :
    fprintf(stdout, "q.%-6.6s: %12s %12s\n%9s %12g %12g\n"
            "%9s %12s %12s\n%9s %12g %12g\n"
            , s, "z[e-3m]", "w0[e-6m]", "", get_z(q), get_w0(q) * 1e3
            , "", "R[e-3m]", "wz[e-6m]", "", r, get_wz(q) * 1e3);
}

void q_fout(const qbeam *q, const char *s)
{
    time_t tt;
    time(&tt);
    FILE *fout;
    fout = fopen("q_fout.txt", "aw");
    if(NULL == fout)
    {
        error_msg("can not open 'q_fout.txt'", __FILE__, __LINE__);
        fclose(fout);
    }
    double r = get_beamrad(q);
    (r == DBL_MAX) ?
    fprintf(fout, "%s\nq.%-6.6s: %12s %12s\n%9s %12g %12g\n"
            "%9s %12s %12s\n%9s %12s %12g\n", ctime(&tt)
            , s, "z[e-3m]", "w0[e-6m]", "", get_z(q), get_w0(q) * 1e3
            , "", "R[e-3m]", "wz[e-6m]", "", "inf", get_wz(q) * 1e3) :
    fprintf(fout, "%s\nq.%-6.6s: %12s %12s\n%9s %12g %12g\n"
            "%9s %12s %12s\n%9s %12g %12g\n", ctime(&tt)
            , s, "z[e-3m]", "w0[e-6m]", "", get_z(q), get_w0(q) * 1e3
            , "", "R[e-3m]", "wz[e-6m]", "", r, get_wz(q) * 1e3);
    fclose(fout);
}
