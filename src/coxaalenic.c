#include "coxinterval.h"
#include "rcplex.h"

static int n, p, q, d, dq;
static double *grad1c, *grad2c, *grad3c, *grad1L, *grad2L;

static double
loglik(double *c, double *L, double *z, double *w, int *J)
{
  int i, j, k, l, m;
  double zc[n], wL[n*d], A[n*d], S[n*d], U[n*d], prob[n], grad1[n], grad2[n],
    grad3[n], ll = 0, temp;
  /* initialize derivative vectors, matrices */
  for (j = 0; j < p; j++) {
    grad1c[j] = 0;
    grad3c[j] = 0;
    for (k = 0; k < p; k++)
      grad2c[j + k*p] = 0;
  }
  for (j = 0; j < dq; j++) {
    grad1L[j] = 0;
    for (k = 0; k < dq; k++)
      grad2L[j + k*dq] = 0;
  }
  /* for each contribution */
  for (i = 0; i < n; i++) {
    zc[i] = 0;
    for (j = 0; j < p; j++)
      zc[i] += z[i + j*n] * c[j];
    /* for each maximal support point */
    for (k = 0; k < d; k++) {
      wL[i + k*n] = 0;
      for (j = 0; j < q; j++)
        wL[i + k*n] += w[i + j*n] * L[j + k*q];
      /* cumulative intensity */
      A[i + k*n] = wL[i + k*n] * exp(zc[i]);
      /* survivor function */
      S[i + k*n] = exp(-A[i + k*n]);
      /* 2nd order derivative of contribution wrt regression function,
         without scaling by w^{\otimes 2} exp(z'coef) */
      U[i + k*n] = (J[i + (k+1)*n] - J[i + k*n]) * S[i + k*n];
    }
    /* left-censored */
    prob[i] = J[i] * (1 - S[i]);
    /* rth order derivatives of contribution wrt regression coefficient c,
       without scaling by z^{\otimes r} */
    grad1[i] = J[i] * A[i] * S[i];
    grad2[i] = J[i] * (1 - A[i]) * A[i] * S[i];
    grad3[i] = J[i] * (1 - 3*A[i] + pow(A[i], 2)) * A[i] * S[i];
    for (k = 1; k < d; k++) {
      /* interval-censored */
      prob[i] += J[i + k*n] * (S[i + (k-1)*n] - S[i + k*n]);
      grad1[i] += J[i + k*n]
        * (A[i + k*n] * S[i + k*n] - A[i + (k-1)*n] * S[i + (k-1)*n]);
      grad2[i] += J[i + k*n] * (A[i + k*n] * (1 - A[i + k*n]) * S[i + k*n]
                               - A[i + (k-1)*n] * (1 - A[i + (k-1)*n])
                               * S[i + (k-1)*n]);
      grad3[i] += J[i + k*n] * ((1 - 3*A[i + k*n] + pow(A[i + k*n], 2))
                              * A[i + k*n] * S[i + k*n]
                              - (1 - 3*A[i + (k-1)*n] + pow(A[i + (k-1)*n], 2))
                              * A[i + (k-1)*n] * S[i + (k-1)*n]);
    }
    /* right-censored */
    prob[i] += J[i+d*n] * S[i + (d-1)*n];
    grad1[i] -= J[i+d*n] * A[i + (d-1)*n] * S[i + (d-1)*n];
    grad2[i] -= J[i+d*n] * A[i + (d-1)*n] * (1 - A[i + (d-1)*n])
      * S[i + (d-1)*n];
    grad3[i] -= J[i+d*n] * (1 - 3*A[i + (d-1)*n] + pow(A[i + (d-1)*n], 2))
      * A[i + (d-1)*n] * S[i + (d-1)*n];
    ll += log(prob[i]) / n;
    /* negative derivatives wrt regression coefficient c */
    for (j = 0; j < p; j++) {
      grad1c[j] -= z[i + j*n] * grad1[i] / (n * prob[i]);
      for (k = 0; k < p; k++) {
        grad2c[j + k*p] -= z[i + j*n] * z[i + k*n]
          * (grad2[i] / prob[i] - pow(grad1[i] / prob[i], 2)) / n;
        grad3c[j] -= z[i + j*n] * pow(z[i + k*n], 2)
          * (grad3[i] / prob[i] - 3 * grad2[i] * grad1[i] / pow(prob[i], 2)
             + 2 * pow(grad1[i] / prob[i], 3)) / n;
      }
    }
    /* negative derivatives wrt cumulative regression function L */
    for (j = 0; j < q; j++)
      for (l = 0; l < d; l++) {
        grad1L[j + l*q] += w[i + j*n] * exp(zc[i]) * U[i + l*n]
          / (n * prob[i]);
        for (k = j; k < q; k++)
          for (m = l; m < d; m++) {
            temp = w[i + j*n] * w[i + k*n] * exp(2*zc[i]) * U[i + l*n]
              * (U[i + m*n] - ((l == m) ? prob[i] : 0))
              / (n * pow(prob[i], 2));
            grad2L[j + l*q + (k + m*q)*dq] += temp;
            grad2L[k + m*q + (j + l*q)*dq] += temp;
          }
      }
  }
  return ll;
}

void
coxaalenic(double *c, double *L, int *nrowL, double *z, int *nrowz, int *ncolz,
           double *w, int *ncolw, int *J, double *A, int *nrowA, double *eps,
           int *epsnorm, int *maxiter, double *armijo, double *typc,
           double *supc, int *trace, int *maxthread, double *var, double *ll,
           int *numiter, double *maxnorm, double *gradnorm, double *cputime,
           int *flag)
{
  clock_t begtime, endtime;
  char uplo = 'U', errmsg[1024];
  int i, j, k, l, m, status = 0, r = *nrowA, iter = 0, steppow, screen = *trace,
    threads = *maxthread;
  double stepinc, *stepc, *candc, *fixc, *stepL, *candL, *pL, *lp, candll, pll,
    *pllvec, *pllmat, *curv, *delta;
  n = *nrowz;
  p = *ncolz;
  q = *ncolw;
  d = *nrowL;
  dq = d*q;
  r = *nrowA;
  grad1c = Calloc(p, double);
  grad2c = Calloc(p*p, double);
  grad3c = Calloc(p, double);
  grad1L = Calloc(dq, double);
  grad2L = Calloc(dq*dq, double);
  stepc = Calloc(p, double);
  candc = Calloc(p, double);
  lp = Calloc(dq, double);
  stepL = Calloc(dq, double);
  candL = Calloc(dq, double);
  fixc = Calloc(p, double);
  pL = Calloc(dq, double);
  curv = Calloc(p*p, double);
  delta = Calloc(p*p, double);
  pllvec = Calloc(p, double);
  pllmat = Calloc(p*p, double);
  begtime = clock();
  ll[iter] = loglik(c, L, z, w, J);
  do { /* constrained Netwon-Raphson method */
    stepinc = 0;
    /* Netwon-Raphson step for regression coefficient c */
    F77_CALL(dpotrf)(&uplo, &p, grad2c, &p, &status);
    if (status) {
      *flag = 1;
      goto deallocate;
    }
    F77_CALL(dpotri)(&uplo, &p, grad2c, &p, &status);
    if (status) {
      *flag = 1;
      goto deallocate;
    }
    for (i = 0; i < p; i++) {
      for (j = 0; j < i; j++)
        grad2c[i + j*p] = grad2c[j + i*p];
      stepc[i] = 0;
      for (j = 0; j < p; j++)
        stepc[i] -= grad2c[i + j*p] * grad1c[j];
      stepinc -= grad1c[i] * stepc[i];
    }
    /* (negative) linear objective coefficient in QP */
    for (i = 0; i < dq; i++) {
      lp[i] = grad1L[i];
      for (j = 0; j < dq; j++)
        lp[i] -= grad2L[i + j*dq] * L[j];
    }
    /* candidate step for cumulative regression function L */
    status = qpcplex(dq, r, lp, grad2L, A, candL, screen, threads, errmsg);
    if (status) {
      *flag = 2;
      goto deallocate;
    }
    for (i = 0; i < dq; i++) {
      stepL[i] = candL[i] - L[i];
      stepinc -= grad1L[i] * stepL[i];
    }
    stepinc *= *armijo;
    steppow = -1;
    do { /* line search */
      ++steppow;
      for (i = 0; i < p; i++)
        candc[i] = c[i] + stepc[i] * pow(0.5, steppow);
      for (i = 0; i < dq; i++)
        candL[i] = L[i] + stepL[i] * pow(0.5, steppow);
      candll = loglik(candc, candL, z, w, J);
    } while (stepinc * pow(0.5, steppow) > candll - ll[iter]);
    /* update estimates */
    *maxnorm = 0;
    *gradnorm = 0;
    for (i = 0; i < p; i++) {
      *maxnorm = max(*maxnorm, fabs(c[i] - candc[i]));
      c[i] = candc[i];
      *gradnorm -= grad1c[i] * c[i];
    }
    for (i = 0; i < dq; i++) {
      *maxnorm = max(*maxnorm, fabs(L[i] - candL[i]));
      L[i] = candL[i];
      *gradnorm -= grad1L[i] * L[i];
    }
    ++iter;
    ll[iter] = candll;
  } while ((*epsnorm ? *maxnorm : fabs(*gradnorm)) > *eps && iter < *maxiter);
  *numiter = iter;
  for (i = 0; i < p; i++) { /* curvature scale */
    fixc[i] = c[i];
    for (j = 0; j <= i; j++)
      curv[i + j*p] = cbrt(-2 * ll[*numiter]
                           / (fabs(grad3c[i]) + (i < j) * fabs(grad3c[j])));
  }
  for (i = 0; i < p; i++)
    for (j = 0; j <= i; j++)
      delta[i + j*p] =
        copysign(max(max(min(curv[i + j*p], *supc), *typc),
                     max(fabs(c[i]), fabs(c[j]))), curv[i + j*p]) / sqrt(n);
  for (i = 0; i < p; i++) { /* approximate profile information */
    for (j = 0; j <= i; j++) {
      fixc[i] += delta[i + j*p];
      for (k = 0; k <= (i == j); k++) {
        fixc[j] += delta[i + j*p];
        for (l = 0; l < dq; l++)
          pL[l] = L[l];
        iter = 0;
        candll = loglik(fixc, pL, z, w, J);
        do { /* until convergence in profile loglikelihood */
          pll = candll;
          /* (negative) linear objective coefficient in QP */
          for (l = 0; l < dq; l++) {
            lp[l] = grad1L[l];
            for (m = 0; m < dq; m++)
              lp[l] -= grad2L[l + m*dq] * pL[m];
          }
          /* candidate step for profile likelihood maximizer pL */
          status =
            qpcplex(dq, r, lp, grad2L, A, candL, screen, threads, errmsg);
          if (status) goto outer;
          stepinc = 0;
          for (l = 0; l < dq; l++) {
            stepL[l] = candL[l] - pL[l];
            stepinc -= grad1L[l] * stepL[l];
          }
          stepinc *= *armijo;
          steppow = -1;
          do { /* line search */
            ++steppow;
            for (l = 0; l < dq; l++)
              candL[l] = pL[l] + stepL[l] * pow(0.5, steppow);
            candll = loglik(fixc, candL, z, w, J);
          } while (stepinc * pow(0.5, steppow) > candll - pll);
          for (l = 0; l < dq; l++)
            pL[l] = candL[l];
          ++iter;
        } while (fabs(1 - pll / candll) > *eps && iter < *maxiter);
        if (k == 0) pllmat[i + j*p] = candll;
        else pllvec[i] = candll;
        fixc[j] = c[j];
      }
      fixc[i] = c[i];
    }
  }
  for (i = 0; i < p; i++)
    for (j = 0; j <= i; j++) {
      var[i + j*p] = (pllvec[i] - ll[*numiter]) / pow(delta[i + i*p], 2)
        + (pllvec[j] - ll[*numiter]) / pow(delta[j + j*p], 2)
        - (pllmat[i + j*p] - ll[*numiter]) / pow(delta[i + j*p], 2);
      if (j < i) var[j + i*p] = var[i + j*p];
    }
  F77_CALL(dpotrf)(&uplo, &p, var, &p, &status);
  if (status) {
    *flag = 2;
    goto deallocate;
  }
  F77_CALL(dpotri)(&uplo, &p, var, &p, &status);
  if (status) {
    *flag = 2;
    goto deallocate;
  }
outer:
  if (status) {
    *flag = 3;
    goto deallocate;
  }
  for (i = 0; i < p; i++) {
    var[i + i*p] /= n;
    for (j = 0; j < i; j++) {
      var[j + i*p] /= n;
      var[i + j*p] = var[j + i*p];
    }
  }
  endtime = clock();
  *cputime = ((double) (endtime - begtime)) / CLOCKS_PER_SEC;
deallocate:
  Free(grad1c);
  Free(grad2c);
  Free(grad3c);
  Free(grad1L);
  Free(grad2L);
  Free(stepc);
  Free(candc);
  Free(lp);
  Free(stepL);
  Free(candL);
  Free(fixc);
  Free(pL);
  Free(curv);
  Free(delta);
  Free(pllvec);
  Free(pllmat);
  return;
}
