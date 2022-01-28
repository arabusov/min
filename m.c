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
#define ERR_NOMEM -1

real diff (real (*fcn) (real *), int i, real *p, real eps)
{
    real tmp = p[i];
    real f1, f2;
    p[i] += eps;
    f2=(*fcn) (p);
    p[i] -= 2*eps;
    f1=(*fcn) (p);
    p[i] = tmp;
    return (f2-f1)/2.;
}

void report (int iter, real val, real *p, int n_par)
{
    int i;
    printf ("%d %lf", iter, val);
    for (i = 0; i < n_par; i++) {
        printf (" %lf", p[i]);
    }
    printf ("\n");
}

extern int minimize (real (*fcn) (real *), real *p, int n_par)
{
    int i, iter=0;
    real reps=1e-3, stop_eps=1e-6;
    real *grad;
    real init_val, mn[2];
    if (n_par <= 0) return 0;
    grad = malloc (sizeof(real)*n_par);
    if (grad == NULL) return ERR_NOMEM;
    if (mn == NULL) return ERR_NOMEM;
    mn[0] = (*fcn) (p);
    report (0, mn[0], p, n_par);
    do {
        for (i=0; i < n_par; i++) {
            real eps = reps;
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
