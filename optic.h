void thin_lens(double *M, const double f);
void thick_lens(double *M, const double r1, const double r2, const double d, const double n);
void thick_p_ck_lens(double *M, const double r2, const double d, const double n);
void propagation(double *M, const double d);
void refraction(double *M, const double n1, const double n2, const double r);
void reflection(double *M, const double r);
void grin(double *M, const double d, const double n0, const double nmin);

void M_mult_ip(double *Mr, const double *M);

qbeam qtrans(const double *M, const qbeam *q);

double wz(const double w0, const double z, const double zr);
double rayleighr(const double w0, const double lam);
double stab_parm(const double *M);
cplx *eigenval_res(const double *M);
cplx eigenmod_res(const double *M);
void setq(qbeam *q, const double z, const double zr, const double lam);

double get_w0(const qbeam *q);
double get_z(const qbeam *q);
double get_zr(const qbeam *q);
double get_beamrad(const qbeam *q);
double get_wz(const qbeam *q);
double get_div(const qbeam *q);

void q_out(const qbeam *q, const char *s);
void q_fout(const qbeam *q, const char *s);
