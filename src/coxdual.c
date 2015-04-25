#include "coxinterval.h"

static int sv, n, p, K, d[3], D[3 + 1], *sidx, *lidx, *ridx, *vidx;
static double eps, *newh, *grad1c, *grad2c, *grad3c, *grad1h;

/* under sieve: length of overlap between kth subinterval
   on type j partition and (beg, end]
   otherwise: indicator that right-endpoint of kth subinterval
   lies in (beg, end] */
static double
L(double *t, int j, int k, double beg, double end)
{
  if (sv == 0) return (beg < t[k + D[j]] && t[k + D[j]] <= end);
  else return max(0, min(t[k + D[j]], end) - max(t[k - 1 + D[j]], beg));
}

/* piecewise function for type j baseline intensity */
static double
H(double *h, double *t, int j, double beg, double end)
{
  int i;
  double val = 0;
  for (i = 1; i < d[j]; i++)
    val += h[i + D[j]] * L(t, j, i, beg, end);
  return val;
}

static double
loglik(double *c, double *h, double *z, double *t, double *s, double *left,
       double *right, double *u, double *v, int *contrib, int *absorb)
{
  int i, j, k, l;
  double ll = 0, prob, prob1, prob2, A01, A02, A12, a01, a02, a12, h02, h12,
    rsk[3], g1c1[p], g1c2[p], g2c1[p*p], g2c2[p*p], g3c1[p], g3c2[p], g1cr[p],
    g2cr[p*p], g3cr[p], g1cn[p], g2cn[p*p], g3cn[p], g1cd[p], g2cd[p*p],
    g3cd[p], g1c[p], g2c[p*p], g3c[p], beg, end, len, pseg, dseg, rseg, s0, s1,
    g1h1[D[3]], g1h2[D[3]], EdN[D[3]], EY[D[3]], newhn[D[3]], newhd[D[3]], y;
  /* negative derivatives wrt coefficient c */
  for (i = 0; i < p; i++) {
    grad1c[i] = 0;
    grad3c[i] = 0;
    for (j = 0; j < p; j++)
      grad2c[i + j*p] = 0;
  }
  for (i = 0; i < D[3]; i++) {
    /* derivative wrt increments in cumulative baseline intensities h */
    grad1h[i] = 0;
    /* numerator and denominator in EM estimator update */
    newhn[i] = 0;
    newhd[i] = 0;
  }
  /* for each individual */
  for (i = 0; i < n; i++) {
    /* likelihood contribution through 0 -> 1 -> 2 */
    prob1 = 0;
    /* likelihood contribution through 0 -> 2 */
    prob2 = 0;
    for (j = 0; j < 3; j++) {
      rsk[j] = 1;
      for (k = 0; k < p; k++)
        rsk[j] *= exp(z[i + k*n + j*n*p] * c[k]);
    }
    /* kth order derivatives of prob1, prob2 wrt c; k = 1, 2, 3 */
    for (j = 0; j < p; j++) {
      g1c1[j] = 0;
      g1c2[j] = 0;
      for (k = 0; k < p; k++) {
        g2c1[j + k*p] = 0;
        g2c2[j + k*p] = 0;
      }
      g3c1[j] = 0;
      g3c2[j] = 0;
    }
    for (j = 0; j < D[3]; j++) {
      /* derivatives of prob1, prob2 wrt h */
      g1h1[j] = 0;
      g1h2[j] = 0;
      /* conditional expectation of increment in N */
      EdN[j] = 0;
      /* conditional expectation of risk set size */
      EY[j] = 0;
    }
    /* contribution via 0 -> 1 -> 2 */
    a12 = (absorb[i] == 1) ? (h[vidx[i + 2*n] + D[2]] * rsk[2]) : 1;
    /* for each partition point in (L, R] */
    for (j = lidx[i]; contrib[i] != 2 && a12 > 0 && j <= ridx[i]; j++) {
      /* contribution in which 0 -> 1 transition takes place over jth
         subinterval, beginning at 'beg' and ending at 'end' */
      a01 = h[sidx[j]] * rsk[0];
      if (a01 <= 0) continue;
      end = min(right[i], s[j]);
      beg = (sv == 0) ? end : max(left[i], s[j - 1]);
      len = end - beg;
      A01 = H(h, t, 0, u[i], beg) * rsk[0];
      A02 = H(h, t, 1, u[i], beg) * rsk[1];
      A12 = H(h, t, 2, end, v[i]) * rsk[2];
      pseg = exp(-(A01 + A02)) * a01 * exp(-A12) * a12;
      /* derivatives of pseg wrt c, short of scaling by pseg */
      for (k = 0; k < p; k++)
        g1c[k] = z[i + k*n] * (1 - A01)
          - z[i + k*n + n*p] * A02
          + z[i + k*n + 2*n*p] * (absorb[i] - A12);
      for (k = 0; k < p; k++) {
        g3c[k] = 0;
        for (l = 0; l < p; l++) {
          g2c[k + l*p] = -z[i + k*n] * z[i + l*n] * A01
            - z[i + k*n + n*p] * z[i + l*n + n*p] * A02
            - z[i + k*n + 2*n*p] * z[i + l*n + 2*n*p] * A12
            + g1c[k] * g1c[l];
          g3c[k] += -z[i + k*n] * pow(z[i + l*n], 2) * A01
            - z[i + k*n + n*p] * pow(z[i + l*n + n*p], 2) * A02
            - z[i + k*n + 2*n*p] * pow(z[i + l*n + 2*n*p], 2) * A12
            + 3 * g2c[k + l*p] * g1c[l] + g1c[k] * pow(g1c[l], 2);
        }
      }
      /* T01 concentrates at 'end' */
      if (len == 0) {
        prob1 += pseg;
        for (k = 0; k < p; k++) {
          g1c1[k] += g1c[k] * pseg;
          g3c1[k] += g3c[k] * pseg;
          for (l = 0; l < p; l++)
            g2c1[k + l*p] += g2c[k + l*p] * pseg;
        }
        g1h1[sidx[j]] += rsk[0] * pseg / a01;
        EdN[sidx[j]] += pseg;
        for (k = 0; k < 2; k++)
          /* at risk for 0 -> 1, 2 */
          for (l = 1; l <= d[k]; l++) {
            if (h[l + D[k]] <= 0) continue;
            g1h1[l + D[k]] -= rsk[k] * pseg * L(t, k, l, u[i], beg);
            EY[l + D[k]] += pseg * L(t, k, l, u[i], beg);
          }
        /* at risk for 1 -> 2 */
        for (k = 1; k <= d[2]; k++) {
          if (h[k + D[2]] <= 0) continue;
          g1h1[k + D[2]] -= rsk[2] * pseg * L(t, 2, k, end, v[i]);
          EY[k + D[2]] += pseg * L(t, 2, k, end, v[i]);
        }
        g1h1[vidx[i + 2*n] + D[2]] += absorb[i] * rsk[2] * pseg / a12;
        EdN[vidx[i + 2*n] + D[2]] += absorb[i] * pseg;
      }
      /* T01 uniform on (beg, end) */
      else {
        h02 = h[sidx[j + K] + D[1]] * rsk[1];
        h12 = h[sidx[j + 2*K] + D[2]] * rsk[2];
        s0 = exp(-len * (a01 + h02));
        s1 = exp(-len * h12);
        dseg = a01 + h02 - h12;
        rseg = (s1 - s0) / dseg;
        prob1 += pseg * rseg;
        /* derivatives of rseg wrt c */
        for (k = 0; k < p; k++) {
          g1cd[k] = (z[i + k*n] * a01 + z[i + k*n + n*p] * h02
                     - z[i + k*n + 2*n*p] * h12) / dseg;
          g1cn[k] = (z[i + k*n] * a01 + z[i + k*n + n*p] * h02) * len * s0
            - z[i + k*n + 2*n*p] * h12 * len * s1;
          g1cr[k] = g1cn[k] / dseg - rseg * g1cd[k];
        }
        for (k = 0; k < p; k++) {
          g3cd[k] = 0;
          g3cn[k] = 0;
          for (l = 0; l < p; l++) {
            g2cd[k + l*p] =
              (z[i + k*n] * z[i + l*n] * a01
               + z[i + k*n + n*p] * z[i + l*n + n*p] * h02
               - z[i + k*n + 2*n*p] * z[i + l*n + 2*n*p] * h12) / dseg;
            g2cn[k + l*p] = len * s0
              * (z[i + k*n] * z[i + l*n] * a01
                 + z[i + k*n + n*p] * z[i + l*n + n*p] * h02)
              - (z[i + k*n] * a01 + z[i + k*n + n*p] * h02)
              * (z[i + l*n] * a01 + z[i + l*n + n*p] * h02) * pow(len, 2) * s0
              - (z[i + k*n + 2*n*p] * z[i + l*n + 2*n*p] * h12)
              * len * s1 * (1 - h12 * len);
            g2cr[k + l*p] = g2cn[k + l*p] / dseg - g1cn[k] / dseg * g1cd[l]
              - g1cd[k] * g1cn[l] / dseg + 2 * rseg * g1cd[k] * g1cd[l]
              - rseg * g2cd[k + l*p];
            g3cd[k] +=
              (z[i + k*n] * pow(z[i + l*n], 2) * a01
               + z[i + k*n + n*p] * pow(z[i + l*n + n*p], 2) * h02
               - z[i + k*n + 2*n*p] * pow(z[i + l*n + 2*n*p], 2) * h12) / dseg;
            g3cn[k] += len * s0
              * (z[i + k*n] * pow(z[i + l*n], 2) * a01
                 + z[i + k*n + n*p] * pow(z[i + l*n + n*p], 2) * h02
                 - 3 * (z[i + k*n] * z[i + l*n] * a01
                        + z[i + k*n + n*p] * z[i + l*n + n*p] * h02)
                 * (z[i + l*n] * a01 + z[i + l*n + n*p] * h02) * len
                 + (z[i + k*n] * a01 + z[i + k*n + n*p] * h02)
                 * pow(len * (z[i + l*n] * a01 + z[i + l*n + n*p] * h02), 2))
              - len * s1 * z[i + k*n + 2*n*p]
              * (1 - 3 * h12 * len + pow(h12 * len, 2))
              * pow(z[i + l*n + 2*n*p], 2) * h12;
          }
        }
        for (k = 0; k < p; k++) {
          g3cr[k] = g3cn[k] / dseg - rseg * g3cd[k];
          for (l = 0; l < p; l++)
            g3cr[k] += -3 * g2cn[k + l*p] / dseg * g1cd[l]
              - 3 * g2cd[k + l*p] * g1cn[l] / dseg
              + 6 * g1cd[k] * g1cd[l] * g1cn[l] / dseg
              + 6 * rseg * g2cd[k + l*p] * g1cd[l]
              - 6 * rseg * g1cd[k] * pow(g1cd[l], 2);
        }
        for (k = 0; k < p; k++) {
          g1c1[k] += (g1c[k] * rseg + g1cr[k]) * pseg;
          g3c1[k] += (g3c[k] * rseg + g3cr[k]) * pseg;
          for (l = 0; l < p; l++) {
            g2c1[k + l*p] += pseg * (g2c[k + l*p] * rseg + g2cr[k + l*p]
                                     + g1c[k] * g1cr[l] + g1c[l] * g1cr[k]);
            g3c1[k] += 3 * (g2c[k + l*p] * rseg * g1cr[l]
                            + g2cr[k + l*p] * g1c[l] * rseg) * pseg;
          }
        }
        g1h1[sidx[j]] += rsk[0] * pseg * rseg / a01;
        EdN[sidx[j]] += pseg * rseg;
        for (k = 0; k < 2; k++) {
          if (h[sidx[j + k*K] + D[k]] > 0) {
            g1h1[sidx[j + k*K] + D[k]] += rsk[k] * pseg * len * s0 / dseg
              - rsk[k] * pseg * rseg / dseg;
            EY[sidx[j + k*K] + D[k]] += pseg * rseg / dseg
              - pseg * len * s0 / dseg;
          }
          /* at risk for 0 -> 1, 2 */
          for (l = 1; l <= d[k]; l++) {
            if (h[l + D[k]] <= 0) continue;
            g1h1[l + D[k]] -= rsk[k] * pseg * rseg * L(t, k, l, u[i], beg);
            EY[l + D[k]] += pseg * rseg * L(t, k, l, u[i], beg);
          }
        }
        if (h[sidx[j + 2*K] + D[2]] > 0) {
          g1h1[sidx[j + 2*K] + D[2]] += rsk[2] * pseg * rseg / dseg
            - rsk[2] * pseg * len * s1 / dseg;
          EY[sidx[j + 2*K] + D[2]] += pseg * len * s1 / dseg
            - pseg * rseg / dseg;
        }
        /* at risk for 1 -> 2 */
        for (k = 1; k <= d[2]; k++) {
          if (h[k + D[2]] <= 0) continue;
          g1h1[k + D[2]] -= rsk[2] * pseg * rseg * L(t, 2, k, end, v[i]);
          EY[k + D[2]] += pseg * rseg * L(t, 2, k, end, v[i]);
        }
        g1h1[vidx[i + 2*n] + D[2]] += absorb[i] * rsk[2] * pseg * rseg / a12;
        EdN[vidx[i + 2*n] + D[2]] += absorb[i] * pseg * rseg;
      }
    }
    /* contribution via 0 -> 2 */
    a02 = (absorb[i]) ? (h[vidx[i + n] + D[1]] * rsk[1]) : 1;
    if (contrib[i] != 1 && a02 > 0) {
      A01 = H(h, t, 0, u[i], v[i]) * rsk[0];
      A02 = H(h, t, 1, u[i], v[i]) * rsk[1];
      prob2 = exp(-(A01 + A02)) * a02;
      for (j = 0; j < p; j++)
        g1c[j] = -z[i + j*n] * A01 + z[i + j*n + n*p] * (absorb[i] - A02);
      for (j = 0; j < p; j++)
        for (k = 0; k < p; k++)
          g2c[j + k*p] = -z[i + j*n] * z[i + k*n] * A01
            - z[i + j*n + n*p] * z[i + k*n + n*p] * A02 + g1c[j] * g1c[k];
      for (j = 0; j < p; j++) {
        g3c[j] = 0;
        for (k = 0; k < p; k++) {
          g3c[j] += -z[i + j*n] * pow(z[i + k*n], 2) * A01
            - z[i + j*n + n*p] * pow(z[i + k*n + n*p], 2) * A02
            + 3 * g2c[j + k*p] * g1c[k] + g1c[j] * pow(g1c[k], 2);
          g2c2[j + k*p] = g2c[j + k*p] * prob2;
        }
        g1c2[j] = g1c[j] * prob2;
        g3c2[j] = g3c[j] * prob2;
      }
      g1h2[vidx[i + n] + D[1]] += absorb[i] * rsk[1] * prob2 / a02;
      EdN[vidx[i + n] + D[1]] += absorb[i] * prob2;
      for (j = 0; j < 2; j++)
        for (k = 1; k <= d[j]; k++) {
          if (h[k + D[j]] <= 0) continue;
          g1h2[k + D[j]] -= rsk[j] * prob2 * L(t, j, k, u[i], v[i]);
          EY[k + D[j]] += prob2 * L(t, j, k, u[i], v[i]);
        }
    }
    /* combine contributions from each trajectory */
    prob = prob1 + prob2;
    ll += log(prob) / n;
    for (j = 0; j < p; j++) {
      grad1c[j] -= (g1c1[j] + g1c2[j]) / (n * prob);
      grad3c[j] -= (g3c1[j] + g3c2[j]) / (n * prob);
      for (k = 0; k < p; k++) {
        grad2c[j + k*p] -= (g2c1[j + k*p] + g2c2[j + k*p]) / (n * prob)
          - (g1c1[j] + g1c2[j]) * (g1c1[k] + g1c2[k]) / (n * pow(prob, 2));
        grad3c[j] -= -3 * (g2c1[j + k*p] + g2c2[j + k*p])
          * (g1c1[k] + g1c2[k]) / (n * pow(prob, 2))
          + (g1c1[j] + g1c2[j]) * pow(g1c1[k] + g1c2[k], 2) / (n * pow(prob, 3));
      }
    }
    for (j = 0; j < 3; j++)
      for (k = 1; k < d[j]; k++) {
        grad1h[k + D[j]] += (h[k + D[j]] > 0 && h[k + D[j]] <= 1/eps)
          * (g1h1[k + D[j]] + g1h2[k + D[j]]) / (n * prob);
        newhn[k + D[j]] += EdN[k + D[j]] / prob;
        newhd[k + D[j]] += rsk[j] * EY[k + D[j]] / prob;
      }
  }
  /* update EM estimator */
  for (i = 0; i < 3; i++)
    for (j = 1; j < d[i]; j++) {
      if (newhn[j + D[i]] > 0)
        newh[j + D[i]] = newhn[j + D[i]] / newhd[j + D[i]];
      else newh[j + D[i]] = 0;
    }
  return ll;
}

void
coxdual(double *c, double *h, int *dimc, int *dimh, double *t, double *s,
        int *dims, double *z, int *nrow, double *left, double *right, double *u,
        double *v, int *contrib, int *absorb, double *varc, double *ll,
        double *epsilon, int *maxiter, double *typc, double *supc, int *zeroc,
        int *sieve, int *numiter, double *maxnorm, double *gradnorm,
        double *cputime, int *flag)
{
  clock_t begtime, endtime;
  char uplo = 'U';
  int i, j, k, l, m, status = 0, iter = 0, *ipiv, lwork, halving;
  double oldll, newll, *candc, *stepc, *fixc, *candh, *steph, *ph, *curv,
    *delta, *pllvec, *pllmat, *work;
  sv = *sieve;
  eps = *epsilon;
  n = *nrow;
  p = *dimc * (1 - *zeroc);
  K = *dims;
  for (i = 0; i < 3; i++) d[i] = dimh[i];
  D[0] = 0;
  for (i = 1; i <= 3; i++) D[i] = D[i - 1] + d[i - 1];
  grad1c = Calloc(p, double);
  grad2c = Calloc(p * p, double);
  grad3c = Calloc(p, double);
  grad1h = Calloc(D[3], double);
  candc = Calloc(p, double);
  stepc = Calloc(p, double);
  fixc = Calloc(p, double);
  newh = Calloc(D[3] + d[1], double);
  candh = Calloc(D[3] + d[1], double);
  steph = Calloc(D[3] + d[1], double);
  ph = Calloc(D[3] + d[1], double);
  curv = Calloc(p * p, double);
  delta = Calloc(p * p, double);
  pllvec = Calloc(p, double);
  pllmat = Calloc(p * p, double);
  ipiv = Calloc(p, int);
  work = Calloc(p, double);
  /* index of s, from common partition, in type j partition */
  sidx = Calloc(3 * K, int);
  for (i = 0; i < 3; i++)
    for (j = 0; j < d[i]; j++)
      for (k = 0; k < K; k++)
        sidx[k + i*K] += s[k] > t[j + D[i]];
  lidx = Calloc(n, int);
  ridx = Calloc(n, int);
  vidx = Calloc(3 * n, int);
  for (i = 0; i < n; i++) {
    /* indices of smallest L <= s and R <= s, s in common partition */
    for (k = 0; k < K; k++) {
      lidx[i] += contrib[i] != 2 && left[i] > s[k];
      ridx[i] += contrib[i] != 2 && right[i] > s[k];      
    }
    /* index V <= t, t in type-specific partition */
    for (j = 0; j < 3; j++)
      for (k = 0; k < d[j]; k++) vidx[i + j*n] += v[i] > t[k + D[j]];
  }
  begtime = clock();
  ll[iter] = loglik(c, h, z, t, s, left, right, u, v, contrib, absorb);
  do { /* estimate parameters */
    if (*zeroc) goto emstep;
    /* Netwon-Raphson step for regression coefficient c */
    lwork = -1;
    F77_CALL(dsytrf)(&uplo, &p, grad2c, &p, ipiv, work, &lwork, &status);
    if (status) { /* can't query DSYTRF workspace */
      *flag = 1;
      goto deallocate;
    }
    lwork = (int) work[0];
    work = Realloc(work, lwork, double);
    F77_CALL(dsytrf)(&uplo, &p, grad2c, &p, ipiv, work, &lwork, &status);
    if (status) { /* can't factorize coefficient Hessian */
      *flag = 1;
      goto deallocate;
    }
    F77_CALL(dsytri)(&uplo, &p, grad2c, &p, ipiv, work, &status);
    if (status) { /* can't invert coefficient Hessian */
      *flag = 1;
      goto deallocate;
    }
    for (i = 0; i < p; i++) {
      for (j = 0; j < i; j++)
        grad2c[i + j*p] = grad2c[j + i*p];
      stepc[i] = 0;
      for (j = 0; j < p; j++)
        stepc[i] -= grad2c[i + j*p] * grad1c[j];
      candc[i] = c[i] + stepc[i];
    }
  emstep:
    /* EM step for baseline transition intensities h */
    loglik(candc, h, z, t, s, left, right, u, v, contrib, absorb);
    for (i = 0; i < 3; i++)
      for (j = 1; j < d[i]; j++) {
        steph[j + D[i]] = newh[j + D[i]] - h[j + D[i]];
        candh[j + D[i]] = newh[j + D[i]];
      }
    newll = loglik(candc, candh, z, t, s, left, right, u, v, contrib, absorb);
    /* overshoot also characterized by exp(z*c) = Inf => log-likelihood NaN */
    halving = 0;
    while (newll < ll[iter] || ISNAN(newll)) { /* step halving */
      halving = 1;
      for (i = 0; i < p; i++) {
        stepc[i] *= 0.5;
        candc[i] = c[i] + stepc[i];
      }
      for (i = 0; i < 3; i++)
        for (j = 1; j < d[i]; j++) {
          steph[j + D[i]] *= 0.5;
          candh[j + D[i]] = h[j + D[i]] + steph[j + D[i]];
        }
      newll = loglik(candc, candh, z, t, s, left, right, u, v, contrib, absorb);
    }
    ++iter;
    ll[iter] = newll;
    *maxnorm = 0;
    *gradnorm = 0;
    for (i = 0; i < p; i++) {
      *maxnorm = max(*maxnorm, fabs(stepc[i]));
      c[i] = candc[i];
      *gradnorm = max(*gradnorm, fabs(c[i] * grad1c[i]));
    }
    for (i = 0; i < 3; i++)
      for (j = 1; j < d[i]; j++) {
        *maxnorm = max(*maxnorm, fabs(steph[j + D[i]]));
        h[j + D[i]] = candh[j + D[i]];
        *gradnorm = max(*gradnorm, fabs(h[j + D[i]] * grad1h[j + D[i]]));
      }
  } while (*maxnorm > eps && iter < *maxiter);
  if (iter == 1) {
    REprintf("'Converged' after one step. Try another starting value.\n");
    *flag = 1;
    goto deallocate;
  }
  *numiter = iter;
  if (*zeroc) goto deallocate;
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
        for (l = 0; l < 3; l++)
          for (m = 0; m < d[l]; m++)
            ph[m + D[l]] = h[m + D[l]];
        iter = 0;
        newll = loglik(fixc, ph, z, t, s, left, right, u, v, contrib, absorb);
        do { /* evaluate profile loglikelihood */
          oldll = newll;
          for (l = 0; l < 3; l++)
            for (m = 0; m < d[l]; m++) {
              steph[m + D[l]] = newh[m + D[l]] - ph[m + D[l]];
              candh[m + D[l]] = newh[m + D[l]];
            }
          newll =
            loglik(fixc, candh, z, t, s, left, right, u, v, contrib, absorb);
          while (newll < oldll || ISNAN(newll)) { /* step-halving */
            for (l = 0; l < 3; l++)
              for (m = 0; m < d[l]; m++) {
                steph[m + D[l]] *= 0.5;
                candh[m + D[l]] = ph[m + D[l]] + steph[m + D[l]];
              }
            newll =
              loglik(fixc, candh, z, t, s, left, right, u, v, contrib, absorb);
          }
          ++iter;
          for (l = 0; l < 3; l++)
            for (m = 0; m < d[l]; m++)
              ph[m + D[l]] = candh[m + D[l]];
        } while (fabs(1 - oldll / newll) > eps && iter < *maxiter);
        if (k == 0) pllmat[i + j*p] = newll;
        else pllvec[i] = newll;
        fixc[j] = c[j];
      }
      fixc[i] = c[i];
    }
  }
  for (i = 0; i < p; i++)
    for (j = 0; j <= i; j++) {
      varc[i + j*p] = (pllvec[i] - ll[*numiter]) / pow(delta[i + i*p], 2)
        + (pllvec[j] - ll[*numiter]) / pow(delta[j + j*p], 2)
        - (pllmat[i + j*p] - ll[*numiter]) / pow(delta[i + j*p], 2);
      if (j < i) varc[j + i*p] = varc[i + j*p];
    }
  F77_CALL(dpotrf)(&uplo, &p, varc, &p, &status);
  if (status) { /* can't factorize profile information */
    *flag = 2;
    goto deallocate;
  }
  F77_CALL(dpotri)(&uplo, &p, varc, &p, &status);
  if (status) { /* can't invert profile information */
    *flag = 2;
    goto deallocate;
  }
  for (i = 0; i < p; i++) {
    varc[i + i*p] /= n;
    for (j = 0; j < i; j++) {
      varc[j + i*p] /= n;
      varc[i + j*p] = varc[j + i*p];
    }
  }
deallocate:
  endtime = clock();
  *cputime = ((double) (endtime - begtime)) / CLOCKS_PER_SEC;
  Free(sidx);
  Free(lidx);
  Free(ridx);
  Free(vidx);
  Free(grad1c);
  Free(grad2c);
  Free(grad3c);
  Free(grad1h);
  Free(candc);
  Free(stepc);
  Free(fixc);
  Free(newh);
  Free(candh);
  Free(steph);
  Free(ph);
  Free(curv);
  Free(delta);
  Free(pllvec);
  Free(pllmat);
  Free(ipiv);
  Free(work);
  return;
}
