#ifndef UTIL_H
#define UTIL_H 1

#include <cmath>

#ifdef NDEBUG
#define assert_double(x, y, eps)
#else
#define assert_double(x, y, eps) \
  ((void) (util::assert_double(__FILE__, __LINE__, x, y, eps)))
#endif

namespace util {
  
static inline float norm(const float x[]) {
  return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}

static inline float dot(float const * const x, float const * const y) {
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

void assert_double(const char file[], const unsigned int line, const double x, const double x_expected, const double eps);

}

class ErrorFile {

};

#endif
