#ifndef TINKER_H
#define TINKER_H 1

#include <gsl/gsl_integration.h>

struct MF {
  double z;
  double alpha;
};

MF* mf_alloc();
MF* mf_alloc(const double z);
void mf_free(MF* const mf);
void mf_set_redshift(MF* const mf, const double a);
double mf_f(MF const * const mf, const double nu);
double mf_b(const double nu);

#endif

