#ifndef __fadamath_h
#define __fadamath_h

#include "math.h"

namespace GascoigneMath
{
inline double pi()
{
  return 3.14159265358979323846;
}

inline double max(double a, double b) 
{
  if (a>b) return a;
  return b;
}

inline int max_int(int a, int b) 
{
  if (a>b) return a;
  return b;
}

inline double min(double a, double b) 
{
  if (a<b) return a;
  return b;
}

inline int min_int(int a, int b) 
{
  if (a<b) return a;
  return b;
}

inline int abs_int(int a)
{
  if (a>0) return a;
  return -a;
}

/* #define PI 3.14159265358979323846 */

/* #ifndef MIN */
/* #define MIN(a,b) ( ((a)>(b)) ? (b) : (a) ) */
/* #endif */

/* #ifndef ABS */
/* #define ABS(a)   ( ((a)>(0)) ? (a) : (-a) ) */
/* #endif */

/* #ifndef FRAC */
/* #define FRAC(x)  (fabs(x-int(x))) */
/* #endif */

/* int ggt(int n1, int n2); */
}

#endif

