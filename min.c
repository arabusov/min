#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define ERR_NOMEM -1

double diff (double (*fcn) (double *), int i, double *p, double eps)
{
    double tmp = p[i];
    double f1, f2;
    p[i] += eps;
    f2=(*fcn) (p);
    p[i] -= 2*eps;
    f1=(*fcn) (p);
    p[i] = tmp;
    return (f2-f1)/2.;
}

void report (int iter, double val, double *p, int n_par)
{
    int i;
    printf ("%d %lf", iter, val);
    for (i = 0; i < n_par; i++) {
        printf (" %lf", p[i]);
    }
    printf ("\n");
}

int minimize (double (*fcn) (double *), double *p, int n_par)
{
    int i, iter=0;
    double reps=1e-3, stop_eps=1e-6;
    if (n_par <= 0) return 0;
    double *grad = malloc (sizeof(double)*n_par);
    double init_val, mn[2];
    if (grad == NULL) return ERR_NOMEM;
    if (mn == NULL) return ERR_NOMEM;
    mn[0] = (*fcn) (p);
    report (0, mn[0], p, n_par);
    do {
        for (i=0; i < n_par; i++) {
            double eps = reps;
            if (fabs (p[i]) > eps) eps *= fabs (p[i]);
            grad [i] = diff (fcn, i, p, eps);
        }
        for (i=0; i < n_par; i++) {
            p[i] -= grad [i];
        }
        mn[1] = mn [0];
        mn[0] = (*fcn) (p);
        iter++;
        report (iter, mn[0], p, n_par);
    } while (mn[0] < (mn[1]-stop_eps));
    free (grad);
    return 0;
}

void bw (double s, double m, double g, int l, double *re, double *im)
{
    double mpi=0.13957;
    double mp;
    double msq;
    double w;
    double wgs, denom;
    mp = 4*mpi*mpi;
    msq = m*m;
    w = sqrt (s);
    wgs = 0.0;
    if (w>2*mpi) {
        double qs, qm;
        int ipower;
        qs = sqrt (fabs ((s-mp)*s))/w;
        qm = sqrt (fabs ((msq-mp)*msq))/m;
        ipower = 2*l+1;
        wgs = g*(msq/w)*pow(qs/qm, ipower);
    }
/*    return complex (msq, 0.)/ complex (msq-s, -wgs);*/
    denom = (msq-s)*(msq-s) + wgs*wgs;
    *re = msq*(msq-s)/denom;
    *im = msq*wgs/denom;
}

double chi2 (double m, double G, int l,
        double *x, double *Re, double *Im, int n_x)
{
    double *re, *im;
    double chi2=0.;
    int i;
    re = malloc (sizeof(double)*n_x);
    im = malloc (sizeof(double)*n_x);
    if (re == NULL || im == NULL) return -1.;
    for (i = 0; i < n_x; i++) {
        bw (x[i]*x[i], m, G, l, re+i, im+i);
        chi2 += (re[i]-Re[i])*(re[i]-Re[i]) +
            (im[i]-Im[i])*(im[i]-Im[i]);
    }
    free (im);
    free (re);
    return chi2;
}

double *x, *Re, *Im;
int n_x;
void init_x ()
{
    int i;
    n_x = 19;
    x = malloc (sizeof (double)*n_x);
    Re = malloc (sizeof (double)*n_x);
    Im = malloc (sizeof (double)*n_x);
    for (i = 0; i < n_x; i++) {
        x [i] = 0.140 + i * 0.1;
        bw (x[i]*x[i], 0.770, 0.1501, 1, &Re[i], &Im[i]);
    }
}

double fcn (double *p)
{
/*    return (p[0]-2.)*p[0]+(p[1]-2.)*p[1];  */
    return chi2 (p[0], p[1], 1, x, Re, Im, n_x);
}

int main ()
{
    init_x ();
    double p[2] = {1.0, 0.24};
    minimize (fcn, p, 2);
    return 0;
}
