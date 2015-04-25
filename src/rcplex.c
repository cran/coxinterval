#include <R.h>
#include "rcplex.h"

void
freecplex(void)
{
  int status = 0;
  char errmsg[1024];
  if (lp != NULL) {
    status = CPXfreeprob(env, &lp);
    lp = NULL;
    if (status) {
      CPXgeterrorstring(env, status, errmsg);
      error("Could not free CPLEX problem.\n%s\n", errmsg);
    }
  }
  if (env != NULL && closecplex) {
    status = CPXcloseCPLEX(&env);
    env = NULL;
    if (status) {
      CPXgeterrorstring(env, status, errmsg);
      error("Could not close CPLEX.\n%s\n", errmsg);
    }
    else
      Rprintf("Closed CPLEX.\n");
  }
  return;
}

int
qpcplex(const int ncol, const int nrow, const double *c, const double *Q,
        const double *A, double *soln, const int screen, const int threads,
        char *errmsg)
{
  int status = 0, i, j, Abeg[ncol], Aind[nrow * ncol], Acnt[ncol], Qbeg[ncol],
    Qind[ncol * ncol], Qcnt[ncol];
  char *probname = "icsurv", sense[ncol];
  double lb[ncol], ub[ncol], b[nrow];
  if (env == NULL) {
    env = CPXopenCPLEX(&status);
    if (env == NULL) {
      CPXgeterrorstring(env, status, errmsg);
      warning("Could not open CPLEX.\n%s\n", errmsg);
    }
  }
  if (lp == NULL) {
    status = CPXsetintparam(env, CPX_PARAM_SCRIND, screen ? CPX_ON : CPX_OFF);
    if (status)
      warning("Failed to set trace parameter. Screen ouput suppressed.\n");
    status = CPXsetintparam(env, CPX_PARAM_THREADS, threads);
    if (status)
      warning("Failed to set thread parameter. Default value used instead.\n");
    lp = CPXcreateprob(env, &status, probname);
    if (lp == NULL) {
      closecplex = 1;
      warning("Could not create CPLEX problem.\n");
    }
  }
  for (i = 0; i < nrow; i++) {
    sense[i] = 'G';
    b[i] = 0;
  }
  for (j = 0; j < ncol; j++) {
    lb[j] = -CPX_INFBOUND;
    ub[j] = CPX_INFBOUND;
    Abeg[j] = j * nrow;
    Qbeg[j] = j * ncol;
    Acnt[j] = nrow;
    Qcnt[j] = ncol;
    for (i = 0; i < nrow; i++)
      Aind[i + j*nrow] = i;
    for (i = 0; i < ncol; i++)
      Qind[i + j*ncol] = i;
  }
  status = CPXcopylp(env, lp, ncol, nrow, 1, c, b, sense, Abeg, Acnt, Aind, A,
                     lb, ub, NULL);
  if (status) goto geterrmsg;
  status = CPXcopyquad(env, lp, Qbeg, Qcnt, Qind, Q);
  if (status) goto geterrmsg;
  status = CPXqpopt(env, lp);
  if (status) goto geterrmsg;
  status = CPXgetx(env, lp, soln, 0, CPXgetnumcols(env, lp) - 1);
  if (status) goto geterrmsg;
geterrmsg:
  CPXgeterrorstring(env, status, errmsg);
  return status;
}
