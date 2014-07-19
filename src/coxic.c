#include "coxinterval.h"

#define M 3

static int n, p, K, d[M], D[M + 1], *sidx, *lidx, *ridx, *vidx;
static double *newh, *grad1c, *grad2c, *grad3c, *grad1h;

/* piecewise function for baseline intensity */
static double
H(double *h, double *t, int j, double beg, double end)
{
  int i;
  double val = 0;
  for (i = 1; i < d[j]; i++)
    val +=
      h[i + D[j]] * max(0, min(t[i + D[j]], end) - max(t[i - 1 + D[j]], beg));
  return val;
}

/* ... length */
static double
L(double *t, int j, int k, double beg, double end)
{
  double val;
  val = max(0, min(t[k + D[j]], end) - max(t[k - 1 + D[j]], beg));
  return val;
}

static double
loglik(double *c, double *h, double *z, double *t, double *s, double *left,
       double *right, double *u, double *v, int *contrib, int *absorb)
{
  int i, j, k, l;
  double ll = 0, prob, prob1, prob2, A01, A02, A12, a01, a02, a12, h02, h12,
    rsk[M], g1c1[p], g1c2[p], g2c1[p*p], g2c2[p*p], g3c1[p], g3c2[p], g1cr[p],
    g2cr[p*p], g3cr[p], g1cn[p], g2cn[p*p], g3cn[p], g1cd[p], g2cd[p*p],
    g3cd[p], g1c[p], g2c[p*p], g3c[p], beg, end, len, pseg, dseg, rseg, s0, s1,
    g1h1[D[M]], g1h2[D[M]], EdN[D[M]], EY[D[M]], newhn[D[M]], newhd[D[M]], y;
  /* initialize derivatives */
  for (i = 0; i < p; i++) {
    grad1c[i] = 0;
    grad3c[i] = 0;
    for (j = 0; j < p; j++)
      grad2c[i + j*p] = 0;
  }
  for (i = 0; i < D[M]; i++) {
    grad1h[i] = 0;
    newhn[i] = 0;
    newhd[i] = 0;
  }
  /* likelihood contributions */
  for (i = 0; i < n; i++) {
    prob1 = 0;
    prob2 = 0;
    for (j = 0; j < M; j++) {
      rsk[j] = 1;
      for (k = 0; k < p; k++)
        rsk[j] *= exp(z[i + k*n + j*n*p] * c[k]);
    }
    for (j = 0; j < p; j++) {
      g1c1[j] = 0; g1c2[j] = 0;
      for (k = 0; k < p; k++) {
        g2c1[j + k*p] = 0;
        g2c2[j + k*p] = 0;
      }
      g3c1[j] = 0;
      g3c2[j] = 0;
    }
    for (j = 0; j < D[M]; j++) {
      g1h1[j] = 0;
      g1h2[j] = 0;
      EdN[j] = 0;
      EY[j] = 0;
    }
    /* contribution via 0 -> 1 -> 2 */
    for (j = lidx[i]; contrib[i] != 2 && j <= ridx[i]; j++) {
      beg = max(left[i], s[j - 1]);
      end = min(right[i], s[j]);
      len = end - beg;
      A01 = H(h, t, 0, u[i], beg) * rsk[0];
      A02 = H(h, t, 1, u[i], beg) * rsk[1];
      A12 = H(h, t, 2, end, v[i]) * rsk[2];
      a01 = h[sidx[j]] * rsk[0];
      a12 = (absorb[i]) ? (h[vidx[i + 2*n] + D[2]] * rsk[2]) : 1;
      for (k = 0; k < p; k++)
        g1c[k] = z[i + k*n] * (1 - A01)
          - z[i + k*n + n*p] * A02 + z[i + k*n + 2*n*p] * (absorb[i] - A12);
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
      /* L = R = T01 */
      if (len == 0) {
        prob1 += exp(-(A01 + A02)) * a01 * exp(-A12) * a12;
        for (k = 0; k < p; k++) {
          g1c1[k] += g1c[k] * prob1;
          g3c1[k] += g3c[k] * prob1;
          for (l = 0; l < p; l++)
            g2c1[k + l*p] += g2c[k + l*p] * prob1;
        }
        g1h1[sidx[j]] += rsk[0] * exp(-(A01 + A02)) * exp(-A12) * a12;
        EdN[sidx[j]] += prob1;
        for (k = 0; k < M - 1; k++)
          for (l = 1; l <= sidx[j + k*K]; l++) {
            g1h1[l + D[k]] -= rsk[k] * prob1 * L(t, k, l, u[i], end);
            EY[l + D[k]] += prob1 * L(t, k, l, u[i], end);
          }
        for (k = 1; k <= vidx[i + 2*n]; k++) {
          g1h1[k + D[2]] -= rsk[2] * prob1 * L(t, 2, k, end, v[i]);
          EY[k + D[2]] += prob1 * L(t, 2, k, end, v[i]);
        }
        g1h1[vidx[i + 2*n] + D[2]] +=
          absorb[i] * rsk[2] * exp(-(A01 + A02)) * a01 * exp(-A12);
        EdN[vidx[i + 2*n] + D[2]] += absorb[i] * prob1;
      }
      /* T01 in (L, R] */
      else {
        pseg = exp(-(A01 + A02)) * a01 * exp(-A12) * a12;
        h02 = h[sidx[j + K] + D[1]] * rsk[1];
        h12 = h[sidx[j + 2*K] + D[2]] * rsk[2];
        s0 = exp(-len * (a01 + h02));
        s1 = exp(-len * h12);
        dseg = a01 + h02 - h12;
        rseg = (s1 - s0) / dseg;
        prob1 += pseg * rseg;
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
        g1h1[sidx[j]] += rsk[0] * exp(-(A01 + A02)) * rseg * exp(-A12) * a12;
        EdN[sidx[j]] += pseg * rseg;
        for (k = 0; k < M; k++)
          for (l = 1; l < d[k]; l++) {
            if (L(t, k, l, beg, end) > 0) {
              y = (k < M - 1) ? t[l - 1 + D[k]] : min(v[i], t[l + D[k]]);
              g1h1[l + D[k]] -= ((k < M - 1) ? 1 : -1)
                * rsk[k] * pseg / pow(dseg, 2)
                * (s1 * (1 + (beg - y) * dseg) - s0 * (1 + (end - y) * dseg));
              EY[l + D[k]] += ((k < M - 1) ? 1 : -1) * pseg / pow(dseg, 2)
                * (s1 * (1 + (beg - y) * dseg) - s0 * (1 + (end - y) * dseg));
            }
            else {
              y = (t[l + D[k]] <= beg && k < M - 1) * L(t, k, l, u[i], beg)
                + ((end <= t[l - 1 + D[k]] && k == M - 1)
                   * L(t, k, l, u[i], v[i]));
              g1h1[l + D[k]] -= rsk[k] * y * pseg * rseg;
              EY[l + D[k]] += y * pseg * rseg;
            }
          }
        g1h1[vidx[i + 2*n] + D[2]] +=
          absorb[i] * rsk[2] * exp(-(A01 + A02)) * a01 * rseg * exp(-A12);
        EdN[vidx[i + 2*n] + D[2]] += absorb[i] * pseg * rseg;
      }
    }
    /* contribution via 0 -> 2 */
    if (contrib[i] != 1) {
      A01 = H(h, t, 0, u[i], v[i]) * rsk[0];
      A02 = H(h, t, 1, u[i], v[i]) * rsk[1];
      a02 = (absorb[i]) ? (h[vidx[i + n] + D[1]] * rsk[1]) : 1;
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
      g1h2[vidx[i + n] + D[1]] += absorb[i] * rsk[1] * exp(-(A01 + A02));
      EdN[vidx[i + n] + D[1]] += absorb[i] * prob2;
      for (j = 0; j < M - 1; j++)
        for (k = 1; k <= vidx[i + j*n]; k++) {
          g1h2[k + D[j]] -= rsk[j] * prob2 * L(t, j, k, u[i], v[i]);
          EY[k + D[j]] += prob2 * L(t, j, k, u[i], v[i]);
        }
    }
    /* overall contribution */
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
          + (g1c1[j] + g1c2[j]) * pow(g1c1[k] + g1c2[k], 2)
          / (n * pow(prob, 3));
      }
    }
    for (j = 0; j < M; j++)
      for (k = 1; k < d[j]; k++) {
        grad1h[k + D[j]] += (g1h1[k + D[j]] + g1h2[k + D[j]]) / (n * prob);
        newhn[k + D[j]] += EdN[k + D[j]] / prob;
        newhd[k + D[j]] += rsk[j] * EY[k + D[j]] / prob;
      }
  }
  /* EM estimator */
  for (i = 0; i < M; i++)
    for (j = 1; j < d[i]; j++) {
      if (newhn[j + D[i]] > 0)
        newh[j + D[i]] = newhn[j + D[i]] / newhd[j + D[i]];
      else newh[j + D[i]] = 0;
    }
  return ll;
}

void
coxic(double *c, double *h, int *dimc, int *dimh, double *t, double *s,
      int *dims, double *z, int *nrow, double *left, double *right, double *u,
      double *v, int *contrib, int *absorb, double *varc, double *ll,
      double *eps, int *maxiter, double *typc, double *supc, int *zeroc,
      int *numiter, double *maxnorm, double *gradnorm, double *cputime,
      int *flag)
{
  clock_t begtime, endtime;
  char uplo = 'U';
  int i, j, k, l, m, status = 0, iter = 0, *ipiv, lwork;
  double oldll, newll, *candc, *stepc, *fixc, *candh, *steph, *ph, *curv,
    *delta, *pllvec, *pllmat, *work;
  n = *nrow;
  p = *dimc;
  K = *dims;
  for (i = 0; i < M; i++)
    d[i] = dimh[i];
  D[0] = 0;
  for (i = 1; i <= M; i++)
    D[i] = D[i - 1] + d[i - 1];
  grad1c = Calloc(p, double);
  grad2c = Calloc(p * p, double);
  grad3c = Calloc(p, double);
  grad1h = Calloc(D[M], double);
  candc = Calloc(p, double);
  stepc = Calloc(p, double);
  fixc = Calloc(p, double);
  newh = Calloc(D[M] + d[1], double);
  candh = Calloc(D[M] + d[1], double);
  steph = Calloc(D[M] + d[1], double);
  ph = Calloc(D[M] + d[1], double);
  curv = Calloc(p * p, double);
  delta = Calloc(p * p, double);
  pllvec = Calloc(p, double);
  pllmat = Calloc(p * p, double);
  ipiv = Calloc(p, int);
  work = Calloc(p, double);
  /* index of h(s), s in common partition */
  sidx = Calloc(M * K, int);
  for (i = 0; i < M; i++)
    for (j = 0; j < d[i]; j++)
      for (k = 0; k < K; k++)
        sidx[k + i*K] += t[j + D[i]] < s[k];
  /* indices of smallest L <= s and R <= s, s in common partition */
  lidx = Calloc(n, int);
  ridx = Calloc(n, int);
  /* index V <= t, t in type-specific partition */
  vidx = Calloc(M * n, int);
  for (i = 0; i < n; i++) {
    if (contrib[i] != 2) {
      for (k = K - 1; k > 0 && right[i] <= s[k]; k--) ridx[i] = k;
      for (k = ridx[i]; k > 0 && left[i] <= s[k]; k--) lidx[i] = k;
    }
    for (j = 0; j < M; j++)
      for (k = d[j] - 1; k > 0 && v[i] <= t[k + D[j]]; k--) vidx[i + j*n] = k;
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
    for (i = 0; i < M; i++)
      for (j = 1; j < d[i]; j++) {
        steph[j + D[i]] = newh[j + D[i]] - h[j + D[i]];
        candh[j + D[i]] = newh[j + D[i]];
      }
    newll = loglik(candc, candh, z, t, s, left, right, u, v, contrib, absorb);
    /* overshoot also characterized by exp(z*c) = Inf => log-likelihood NaN */
    while (newll < ll[iter] || ISNAN(newll)) { /* step halving */
      for (i = 0; i < p; i++) {
        stepc[i] *= 0.5;
        candc[i] = c[i] + stepc[i];
      }
      for (i = 0; i < M; i++)
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
      *gradnorm = max(*gradnorm, fabs(grad1c[i]));
    }
    for (i = 0; i < M; i++)
      for (j = 1; j < d[i]; j++) {
        *maxnorm = max(*maxnorm, fabs(steph[j + D[i]]));
        h[j + D[i]] = candh[j + D[i]];
        *gradnorm = max(*gradnorm, fabs(grad1h[j + D[i]]));
      }
  } while (*maxnorm > *eps && iter < *maxiter);
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
        for (l = 0; l < M; l++)
          for (m = 0; m < d[l]; m++)
            ph[m + D[l]] = h[m + D[l]];
        iter = 0;
        newll = loglik(fixc, ph, z, t, s, left, right, u, v, contrib, absorb);
        do { /* evaluate profile loglikelihood */
          oldll = newll;
          for (l = 0; l < M; l++)
            for (m = 0; m < d[l]; m++) {
              steph[m + D[l]] = newh[m + D[l]] - ph[m + D[l]];
              candh[m + D[l]] = newh[m + D[l]];
            }
          newll =
            loglik(fixc, candh, z, t, s, left, right, u, v, contrib, absorb);
          while (newll < oldll || ISNAN(newll)) { /* step-halving */
            for (l = 0; l < M; l++)
              for (m = 0; m < d[l]; m++) {
                steph[m + D[l]] *= 0.5;
                candh[m + D[l]] = ph[m + D[l]] + steph[m + D[l]];
              }
            newll =
              loglik(fixc, candh, z, t, s, left, right, u, v, contrib, absorb);
          }
          ++iter;
          for (l = 0; l < M; l++)
            for (m = 0; m < d[l]; m++)
              ph[m + D[l]] = candh[m + D[l]];
        } while (fabs(1 - oldll / newll) > *eps && iter < *maxiter);
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
