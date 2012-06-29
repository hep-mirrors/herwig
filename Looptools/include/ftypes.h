#ifndef FTYPES_H
#define FTYPES_H

#define FORTRAN(s) s##_

typedef int INTEGER;
typedef const INTEGER CINTEGER;
typedef double DOUBLE_PRECISION;
typedef const DOUBLE_PRECISION CDOUBLE_PRECISION;
typedef struct { DOUBLE_PRECISION re, im; } DOUBLE_COMPLEX;
typedef const DOUBLE_COMPLEX CDOUBLE_COMPLEX;
typedef char CHARACTER;
typedef const CHARACTER CCHARACTER;

#ifdef __cplusplus

#include <complex>
typedef std::complex<double> double_complex;
#define ToComplex(c) double_complex(c.re, c.im)
#define ToComplex2(r,i) double_complex(r, i)
#define Re(x) std::real(x)
#define Im(x) std::imag(x)

#else

typedef DOUBLE_COMPLEX double_complex;
#define ToComplex(c) c
#define ToComplex2(r,i) (double_complex){r, i}
#define Re(x) (x).re
#define Im(x) (x).im

#endif

#endif

