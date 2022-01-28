/*
Copyright 2022 Andrei Rabusov

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "m.h"

void bw (real s, real m, real g, int l, real *re, real *im)
{
    real mpi=0.13957;
    real mp;
    real msq;
    real w;
    real wgs, denom;
    mp = 4*mpi*mpi;
    msq = m*m;
    w = sqrt (s);
    wgs = 0.0;
    if (w>2*mpi) {
        real qs, qm;
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

real chi2 (real m, real G, int l,
        real *x, real *Re, real *Im, int n_x)
{
    real *re, *im;
    real chi2=0.;
    int i;
    re = malloc (sizeof(real)*n_x);
    im = malloc (sizeof(real)*n_x);
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

real *x, *Re, *Im;
int n_x;
void init_x (void)
{
    int i;
    n_x = 19;
    x = malloc (sizeof (real)*n_x);
    Re = malloc (sizeof (real)*n_x);
    Im = malloc (sizeof (real)*n_x);
    for (i = 0; i < n_x; i++) {
        x [i] = 0.140 + i * 0.1;
        bw (x[i]*x[i], 0.770, 0.1501, 1, &Re[i], &Im[i]);
    }
}

real fcn (real *p)
{
/*    return (p[0]-2.)*p[0]+(p[1]-2.)*p[1];  */
    return chi2 (p[0], p[1], 1, x, Re, Im, n_x);
}

int main (void)
{
    real p[2] = {1.0, 0.24};
    init_x ();
    minimize (fcn, p, 2);
    return 0;
}
