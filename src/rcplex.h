#include <R.h>
#include <ilcplex/cplex.h>

CPXENVptr env;
CPXLPptr lp;
int closecplex;

void freecplex(void);

int qpcplex(const int ncol, const int nrow, const double *c, const double *Q,
            const double *A, double *soln, const int screen, const int threads,
            char *errmsg);
