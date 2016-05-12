
#ifndef OCML_H
#define OCML_H


// This C header declares the functions defined by the OCML library
// The behavior of some of these functions can be controlled or optimized
// by the following set, which must be provided by the implementation
// possibly on a per-program basis
//
//    int __ocml_have_fast_fma32(void)
//        - fma on float arguments is faster than separate multiply then add
//
//    int __ocml_have_fast_fma64(void)
//        - fma on double arguments is faster than separate multiply then add
//
//    int __ocml_finite_only_opt(void)
//        - the application will only pass finite arguments and expects only finite results
//
//    int __ocml_fast_relaxed_opt(void)
//        - the aopplication accepts optimizations that may lower the accuracy of the results
//
//    int __ocml_daz_opt(void)
//        - the application allows subnormal inputs or outputs to be flushed to zero
//
//    int __ocml_amd_opt(void)
//        - teh applicaiton is running on an AMD device
//
//    int __ocml_correctly_rounded_sqrt32(void)
//        - the application is expecting sqrt(float) to produce a correctly rounded result
//
// it is expected that the implementation provides these as if compiled from the following
// C code:
//
//     __attribute__((always_inline, const)) int __ocml_...(void) { return 0; /* or 1 */ }
// 
// allowing them and any control flow associated with them to be optimized away


// Define here the return values from fpclassify
// These match most host definitions
#define FP_NAN 0
#define FP_INFINITE 1
#define FP_ZERO 2
#define FP_SUBNORMAL 3
#define FP_NORMAL 4

#define _MANGLE3(P,N,S) P##_##N##_##S
#define MANGLE3(P,N,S) _MANGLE3(P,N,S)
#define OCML_MANGLE_F32(N) MANGLE3(__ocml, N, f32)
#define OCML_MANGLE_F64(N) MANGLE3(__ocml, N, f64)
#define OCML_MANGLE_F16(N) MANGLE3(__ocml, N, f16)

#define DECL_OCML_UNARY_F32(N) extern float OCML_MANGLE_F32(N)(float);
#define _DECL_X_OCML_UNARY_F32(A,N) extern __attribute__((A)) float OCML_MANGLE_F32(N)(float);
#define DECL_PURE_OCML_UNARY_F32(N) _DECL_X_OCML_UNARY_F32(pure, N)
#define DECL_CONST_OCML_UNARY_F32(N) _DECL_X_OCML_UNARY_F32(const, N)

#define DECL_CONST_OCML_UNARYPRED_F32(N) extern __attribute__((const)) int OCML_MANGLE_F32(N)(float);

#define DECL_OCML_BINARY_F32(N) extern float OCML_MANGLE_F32(N)(float, float);
#define _DECL_X_OCML_BINARY_F32(A,N) extern __attribute__((A)) float OCML_MANGLE_F32(N)(float, float);
#define DECL_PURE_OCML_BINARY_F32(N) _DECL_X_OCML_BINARY_F32(pure, N)
#define DECL_CONST_OCML_BINARY_F32(N) _DECL_X_OCML_BINARY_F32(const, N)

#define DECL_CONST_OCML_BINARYPRED_F32(N) extern __attribute__((const)) int OCML_MANGLE_F32(N)(float, float);

#define _DECL_X_OCML_TERNARY_F32(A,N) extern __attribute__((A)) float OCML_MANGLE_F32(N)(float, float, float);
#define DECL_PURE_OCML_TERNARY_F32(N) _DECL_X_OCML_TERNARY_F32(pure, N)
#define DECL_CONST_OCML_TERNARY_F32(N) _DECL_X_OCML_TERNARY_F32(const, N)

#define DECL_OCML_UNARY_F64(N) extern double OCML_MANGLE_F64(N)(double);
#define _DECL_X_OCML_UNARY_F64(A,N) extern __attribute__((A)) double OCML_MANGLE_F64(N)(double);
#define DECL_PURE_OCML_UNARY_F64(N) _DECL_X_OCML_UNARY_F64(pure, N)
#define DECL_CONST_OCML_UNARY_F64(N) _DECL_X_OCML_UNARY_F64(const, N)

#define DECL_CONST_OCML_UNARYPRED_F64(N) extern __attribute__((const)) int OCML_MANGLE_F64(N)(double);

#define DECL_OCML_BINARY_F64(N) extern double OCML_MANGLE_F64(N)(double, double);
#define _DECL_X_OCML_BINARY_F64(A,N) extern __attribute__((A)) double OCML_MANGLE_F64(N)(double, double);
#define DECL_PURE_OCML_BINARY_F64(N) _DECL_X_OCML_BINARY_F64(pure, N)
#define DECL_CONST_OCML_BINARY_F64(N) _DECL_X_OCML_BINARY_F64(const, N)

#define DECL_CONST_OCML_BINARYPRED_F64(N) extern __attribute__((const)) int OCML_MANGLE_F64(N)(double, double);

#define _DECL_X_OCML_TERNARY_F64(A,N) extern __attribute__((A)) double OCML_MANGLE_F64(N)(double, double, double);
#define DECL_PURE_OCML_TERNARY_F64(N) _DECL_X_OCML_TERNARY_F64(pure, N)
#define DECL_CONST_OCML_TERNARY_F64(N) _DECL_X_OCML_TERNARY_F64(const, N)

#define DECL_OCML_UNARY_F16(N) extern half OCML_MANGLE_F16(N)(half);
#define _DECL_X_OCML_UNARY_F16(A,N) extern __attribute__((A)) half OCML_MANGLE_F16(N)(half);
#define DECL_PURE_OCML_UNARY_F16(N) _DECL_X_OCML_UNARY_F16(pure, N)
#define DECL_CONST_OCML_UNARY_F16(N) _DECL_X_OCML_UNARY_F16(const, N)

#define DECL_CONST_OCML_UNARYPRED_F16(N) extern __attribute__((const)) int OCML_MANGLE_F16(N)(half);

#define DECL_OCML_BINARY_F16(N) extern half OCML_MANGLE_F16(N)(half, half);
#define _DECL_X_OCML_BINARY_F16(A,N) extern __attribute__((A)) half OCML_MANGLE_F16(N)(half, half);
#define DECL_PURE_OCML_BINARY_F16(N) _DECL_X_OCML_BINARY_F16(pure, N)
#define DECL_CONST_OCML_BINARY_F16(N) _DECL_X_OCML_BINARY_F16(const, N)

#define DECL_CONST_OCML_BINARYPRED_F16(N) extern __attribute__((const)) int OCML_MANGLE_F16(N)(half, half);

#define _DECL_X_OCML_TERNARY_F16(A,N) extern __attribute__((A)) half OCML_MANGLE_F16(N)(half, half, half);
#define DECL_PURE_OCML_TERNARY_F16(N) _DECL_X_OCML_TERNARY_F16(pure, N)
#define DECL_CONST_OCML_TERNARY_F16(N) _DECL_X_OCML_TERNARY_F16(const, N)

DECL_CONST_OCML_UNARY_F32(acos)
DECL_CONST_OCML_UNARY_F32(acospi)
DECL_PURE_OCML_UNARY_F32(acosh)
DECL_CONST_OCML_UNARY_F32(asin)
DECL_CONST_OCML_UNARY_F32(asinpi)
DECL_PURE_OCML_UNARY_F32(asinh)
DECL_CONST_OCML_BINARY_F32(atan2)
DECL_CONST_OCML_BINARY_F32(atan2pi)
DECL_CONST_OCML_UNARY_F32(atan)
DECL_PURE_OCML_UNARY_F32(atanh)
DECL_CONST_OCML_UNARY_F32(atanpi)
DECL_PURE_OCML_UNARY_F32(cbrt)
DECL_CONST_OCML_UNARY_F32(ceil)
DECL_OCML_UNARY_F32(cos)
DECL_PURE_OCML_UNARY_F32(cosh)
DECL_OCML_UNARY_F32(cospi)
DECL_CONST_OCML_BINARY_F32(copysign)
DECL_PURE_OCML_UNARY_F32(erf)
DECL_PURE_OCML_UNARY_F32(erfc)
DECL_PURE_OCML_UNARY_F32(erfinv)
DECL_PURE_OCML_UNARY_F32(erfcinv)
DECL_PURE_OCML_UNARY_F32(exp)
DECL_PURE_OCML_UNARY_F32(exp2)
DECL_PURE_OCML_UNARY_F32(exp10)
DECL_PURE_OCML_UNARY_F32(expm1)
DECL_CONST_OCML_UNARY_F32(fabs)
DECL_CONST_OCML_BINARY_F32(fdim)
DECL_CONST_OCML_UNARY_F32(floor)
DECL_CONST_OCML_TERNARY_F32(fma)
DECL_CONST_OCML_BINARY_F32(fmax)
DECL_CONST_OCML_BINARY_F32(fmin)
DECL_CONST_OCML_BINARY_F32(fmod)
DECL_CONST_OCML_UNARYPRED_F32(fpclassify)
extern float OCML_MANGLE_F32(fract)(float, __private float *);
extern float OCML_MANGLE_F32(frexp)(float, __private int *);
DECL_CONST_OCML_BINARY_F32(hypot)
DECL_CONST_OCML_UNARYPRED_F32(ilogb)
DECL_CONST_OCML_UNARYPRED_F32(isfinite)
DECL_CONST_OCML_UNARYPRED_F32(isinf)
DECL_CONST_OCML_UNARYPRED_F32(isnan)
DECL_CONST_OCML_UNARYPRED_F32(isnormal)
DECL_OCML_UNARY_F32(j0)
DECL_OCML_UNARY_F32(j1)
extern __attribute__((const)) float OCML_MANGLE_F32(ldexp)(float, int);
DECL_OCML_UNARY_F32(lgamma)
extern float OCML_MANGLE_F32(lgamma_r)(float, __private int *);
DECL_PURE_OCML_UNARY_F32(log)
DECL_PURE_OCML_UNARY_F32(log2)
DECL_PURE_OCML_UNARY_F32(log10)
DECL_PURE_OCML_UNARY_F32(log1p)
DECL_CONST_OCML_UNARY_F32(logb)
DECL_CONST_OCML_TERNARY_F32(mad)
DECL_CONST_OCML_BINARY_F32(max)
DECL_CONST_OCML_BINARY_F32(min)
DECL_CONST_OCML_BINARY_F32(maxmag)
DECL_CONST_OCML_BINARY_F32(minmag)
extern float OCML_MANGLE_F32(modf)(float, __private float *);
extern __attribute__((const)) float OCML_MANGLE_F32(nan)(uint);
DECL_CONST_OCML_UNARY_F32(nearbyint)
DECL_CONST_OCML_BINARY_F32(nextafter)
DECL_PURE_OCML_BINARY_F32(pow)
DECL_PURE_OCML_BINARY_F32(powr)
extern __attribute__((pure)) float OCML_MANGLE_F32(pown)(float, int);
extern __attribute__((pure)) float OCML_MANGLE_F32(rootn)(float, int);
DECL_CONST_OCML_BINARY_F32(remainder)
extern float OCML_MANGLE_F32(remquo)(float, float, __private int *);
DECL_CONST_OCML_UNARY_F32(rint)
DECL_CONST_OCML_UNARY_F32(round)
DECL_PURE_OCML_UNARY_F32(rcbrt)
DECL_PURE_OCML_UNARY_F32(rsqrt)
DECL_CONST_OCML_BINARY_F32(scalb)
extern __attribute__((const)) float OCML_MANGLE_F32(scalbn)(float, int);
DECL_CONST_OCML_UNARYPRED_F32(signbit)
DECL_OCML_UNARY_F32(sin)
DECL_PURE_OCML_UNARY_F32(sinh)
DECL_OCML_UNARY_F32(sinpi)
extern float OCML_MANGLE_F32(sincos)(float, __private float *);
DECL_CONST_OCML_UNARY_F32(sqrt)
DECL_OCML_UNARY_F32(tan)
DECL_CONST_OCML_UNARY_F32(tanpi)
DECL_PURE_OCML_UNARY_F32(tanh)
DECL_OCML_UNARY_F32(tgamma)
DECL_CONST_OCML_UNARY_F32(trunc)
DECL_OCML_UNARY_F32(y0)
DECL_OCML_UNARY_F32(y1)

DECL_CONST_OCML_BINARY_F32(add_rte)
DECL_CONST_OCML_BINARY_F32(add_rtp)
DECL_CONST_OCML_BINARY_F32(add_rtn)
DECL_CONST_OCML_BINARY_F32(add_rtz)

DECL_CONST_OCML_BINARY_F32(div_rte)
DECL_CONST_OCML_BINARY_F32(div_rtp)
DECL_CONST_OCML_BINARY_F32(div_rtn)
DECL_CONST_OCML_BINARY_F32(div_rtz)

DECL_CONST_OCML_TERNARY_F32(fma_rte)
DECL_CONST_OCML_TERNARY_F32(fma_rtp)
DECL_CONST_OCML_TERNARY_F32(fma_rtn)
DECL_CONST_OCML_TERNARY_F32(fma_rtz)

DECL_CONST_OCML_BINARY_F32(mul_rte)
DECL_CONST_OCML_BINARY_F32(mul_rtp)
DECL_CONST_OCML_BINARY_F32(mul_rtn)
DECL_CONST_OCML_BINARY_F32(mul_rtz)

DECL_CONST_OCML_UNARY_F32(sqrt_rte)
DECL_CONST_OCML_UNARY_F32(sqrt_rtp)
DECL_CONST_OCML_UNARY_F32(sqrt_rtn)
DECL_CONST_OCML_UNARY_F32(sqrt_rtz)

DECL_CONST_OCML_BINARY_F32(sub_rte)
DECL_CONST_OCML_BINARY_F32(sub_rtp)
DECL_CONST_OCML_BINARY_F32(sub_rtn)
DECL_CONST_OCML_BINARY_F32(sub_rtz)


DECL_CONST_OCML_UNARY_F64(acos)
DECL_PURE_OCML_UNARY_F64(acosh)
DECL_CONST_OCML_UNARY_F64(acospi)
DECL_CONST_OCML_UNARY_F64(asin)
DECL_PURE_OCML_UNARY_F64(asinh)
DECL_CONST_OCML_UNARY_F64(asinpi)
DECL_CONST_OCML_UNARY_F64(atan)
DECL_PURE_OCML_UNARY_F64(atanh)
DECL_CONST_OCML_UNARY_F64(atanpi)
DECL_CONST_OCML_BINARY_F64(atan2)
DECL_CONST_OCML_BINARY_F64(atan2pi)
DECL_PURE_OCML_UNARY_F64(cbrt)
DECL_CONST_OCML_UNARY_F64(ceil)
DECL_CONST_OCML_BINARY_F64(copysign)
DECL_OCML_UNARY_F64(cos)
DECL_PURE_OCML_UNARY_F64(cosh)
DECL_OCML_UNARY_F64(cospi)
DECL_PURE_OCML_UNARY_F64(erf)
DECL_PURE_OCML_UNARY_F64(erfc)
DECL_PURE_OCML_UNARY_F64(erfinv)
DECL_PURE_OCML_UNARY_F64(erfcinv)
DECL_PURE_OCML_UNARY_F64(exp)
DECL_PURE_OCML_UNARY_F64(exp2)
DECL_PURE_OCML_UNARY_F64(exp10)
DECL_PURE_OCML_UNARY_F64(expm1)
DECL_CONST_OCML_UNARY_F64(fabs)
DECL_CONST_OCML_BINARY_F64(fdim)
DECL_CONST_OCML_UNARY_F64(floor)
DECL_CONST_OCML_TERNARY_F64(fma)
DECL_CONST_OCML_BINARY_F64(fmax)
DECL_CONST_OCML_BINARY_F64(fmin)
DECL_CONST_OCML_BINARY_F64(fmod)
DECL_CONST_OCML_UNARYPRED_F64(fpclassify)
extern double OCML_MANGLE_F64(fract)(double, __private double *);
extern double OCML_MANGLE_F64(frexp)(double, __private int *);
DECL_CONST_OCML_BINARY_F64(hypot)
DECL_CONST_OCML_UNARYPRED_F64(ilogb)
DECL_CONST_OCML_UNARYPRED_F64(isfinite)
DECL_CONST_OCML_UNARYPRED_F64(isinf)
DECL_CONST_OCML_UNARYPRED_F64(isnan)
DECL_CONST_OCML_UNARYPRED_F64(isnormal)
DECL_OCML_UNARY_F64(j0)
DECL_OCML_UNARY_F64(j1)
extern __attribute__((const)) double OCML_MANGLE_F64(ldexp)(double, int);
DECL_OCML_UNARY_F64(lgamma)
extern double OCML_MANGLE_F64(lgamma_r)(double, __private int *);
DECL_PURE_OCML_UNARY_F64(log)
DECL_PURE_OCML_UNARY_F64(log2)
DECL_PURE_OCML_UNARY_F64(log10)
DECL_PURE_OCML_UNARY_F64(log1p)
DECL_CONST_OCML_UNARY_F64(logb)
DECL_CONST_OCML_TERNARY_F64(mad)
DECL_CONST_OCML_BINARY_F64(max)
DECL_CONST_OCML_BINARY_F64(min)
DECL_CONST_OCML_BINARY_F64(maxmag)
DECL_CONST_OCML_BINARY_F64(minmag)
extern double OCML_MANGLE_F64(modf)(double, __private double *);
extern __attribute__((const)) double OCML_MANGLE_F64(nan)(ulong);
DECL_CONST_OCML_UNARY_F64(nearbyint)
DECL_CONST_OCML_BINARY_F64(nextafter)
DECL_PURE_OCML_BINARY_F64(pow)
DECL_PURE_OCML_BINARY_F64(powr)
extern __attribute__((pure)) double OCML_MANGLE_F64(pown)(double, int);
extern __attribute__((pure)) double OCML_MANGLE_F64(rootn)(double, int);
DECL_CONST_OCML_BINARY_F64(remainder)
extern double OCML_MANGLE_F64(remquo)(double, double, __private int *);
DECL_CONST_OCML_UNARY_F64(rint)
DECL_CONST_OCML_UNARY_F64(round)
DECL_PURE_OCML_UNARY_F64(rcbrt)
DECL_PURE_OCML_UNARY_F64(rsqrt)
DECL_CONST_OCML_BINARY_F64(scalb)
extern __attribute__((const)) double OCML_MANGLE_F64(scalbn)(double, int);
DECL_CONST_OCML_UNARYPRED_F64(signbit)
DECL_OCML_UNARY_F64(sin)
extern double OCML_MANGLE_F64(sincos)(double, __private double *);
DECL_PURE_OCML_UNARY_F64(sinh)
DECL_OCML_UNARY_F64(sinpi)
DECL_CONST_OCML_UNARY_F64(sqrt)
DECL_OCML_UNARY_F64(tan)
DECL_PURE_OCML_UNARY_F64(tanh)
DECL_CONST_OCML_UNARY_F64(tanpi)
DECL_OCML_UNARY_F64(tgamma)
DECL_CONST_OCML_UNARY_F64(trunc)
DECL_OCML_UNARY_F64(y0)
DECL_OCML_UNARY_F64(y1)

DECL_CONST_OCML_BINARY_F64(add_rte)
DECL_CONST_OCML_BINARY_F64(add_rtp)
DECL_CONST_OCML_BINARY_F64(add_rtn)
DECL_CONST_OCML_BINARY_F64(add_rtz)

DECL_CONST_OCML_BINARY_F64(div_rte)
DECL_CONST_OCML_BINARY_F64(div_rtp)
DECL_CONST_OCML_BINARY_F64(div_rtn)
DECL_CONST_OCML_BINARY_F64(div_rtz)

DECL_CONST_OCML_TERNARY_F64(fma_rte)
DECL_CONST_OCML_TERNARY_F64(fma_rtp)
DECL_CONST_OCML_TERNARY_F64(fma_rtn)
DECL_CONST_OCML_TERNARY_F64(fma_rtz)

DECL_CONST_OCML_BINARY_F64(mul_rte)
DECL_CONST_OCML_BINARY_F64(mul_rtp)
DECL_CONST_OCML_BINARY_F64(mul_rtn)
DECL_CONST_OCML_BINARY_F64(mul_rtz)

DECL_CONST_OCML_UNARY_F64(sqrt_rte)
DECL_CONST_OCML_UNARY_F64(sqrt_rtp)
DECL_CONST_OCML_UNARY_F64(sqrt_rtn)
DECL_CONST_OCML_UNARY_F64(sqrt_rtz)

DECL_CONST_OCML_BINARY_F64(sub_rte)
DECL_CONST_OCML_BINARY_F64(sub_rtp)
DECL_CONST_OCML_BINARY_F64(sub_rtn)
DECL_CONST_OCML_BINARY_F64(sub_rtz)

#pragma OPENCL EXTENSION cl_khr_fp16 : enable
DECL_CONST_OCML_UNARY_F16(acos)
DECL_PURE_OCML_UNARY_F16(acosh)
DECL_CONST_OCML_UNARY_F16(acospi)
DECL_CONST_OCML_UNARY_F16(asin)
DECL_PURE_OCML_UNARY_F16(asinh)
DECL_CONST_OCML_UNARY_F16(asinpi)
DECL_CONST_OCML_UNARY_F16(atan)
DECL_PURE_OCML_UNARY_F16(atanh)
DECL_CONST_OCML_UNARY_F16(atanpi)
DECL_CONST_OCML_BINARY_F16(atan2)
DECL_CONST_OCML_BINARY_F16(atan2pi)
DECL_PURE_OCML_UNARY_F16(cbrt)
DECL_CONST_OCML_UNARY_F16(ceil)
DECL_CONST_OCML_BINARY_F16(copysign)
DECL_OCML_UNARY_F16(cos)
DECL_PURE_OCML_UNARY_F16(cosh)
DECL_OCML_UNARY_F16(cospi)
DECL_PURE_OCML_UNARY_F16(erf)
DECL_PURE_OCML_UNARY_F16(erfc)
DECL_PURE_OCML_UNARY_F16(erfinv)
DECL_PURE_OCML_UNARY_F16(erfcinv)
DECL_PURE_OCML_UNARY_F16(exp)
DECL_PURE_OCML_UNARY_F16(exp2)
DECL_PURE_OCML_UNARY_F16(exp10)
DECL_PURE_OCML_UNARY_F16(expm1)
DECL_CONST_OCML_UNARY_F16(fabs)
DECL_CONST_OCML_BINARY_F16(fdim)
DECL_CONST_OCML_UNARY_F16(floor)
DECL_CONST_OCML_TERNARY_F16(fma)
DECL_CONST_OCML_TERNARY_F16(fma_rte)
DECL_CONST_OCML_TERNARY_F16(fma_rtp)
DECL_CONST_OCML_TERNARY_F16(fma_rtn)
DECL_CONST_OCML_TERNARY_F16(fma_rtz)
DECL_CONST_OCML_BINARY_F16(fmax)
DECL_CONST_OCML_BINARY_F16(fmin)
DECL_CONST_OCML_BINARY_F16(fmod)
DECL_CONST_OCML_UNARYPRED_F16(fpclassify)
extern half OCML_MANGLE_F16(fract)(half, __private half *);
extern half OCML_MANGLE_F16(frexp)(half, __private int *);
DECL_CONST_OCML_BINARY_F16(hypot)
DECL_CONST_OCML_UNARYPRED_F16(ilogb)
DECL_CONST_OCML_UNARYPRED_F16(isfinite)
DECL_CONST_OCML_UNARYPRED_F16(isinf)
DECL_CONST_OCML_UNARYPRED_F16(isnan)
DECL_CONST_OCML_UNARYPRED_F16(isnormal)
DECL_OCML_UNARY_F16(j0)
DECL_OCML_UNARY_F16(j1)
extern __attribute__((const)) half OCML_MANGLE_F16(ldexp)(half, int);
DECL_OCML_UNARY_F16(lgamma)
extern half OCML_MANGLE_F16(lgamma_r)(half, __private int *);
DECL_PURE_OCML_UNARY_F16(log)
DECL_CONST_OCML_UNARY_F16(logb)
DECL_PURE_OCML_UNARY_F16(log2)
DECL_PURE_OCML_UNARY_F16(log10)
DECL_PURE_OCML_UNARY_F16(log1p)
DECL_CONST_OCML_TERNARY_F16(mad)
DECL_CONST_OCML_BINARY_F16(max)
DECL_CONST_OCML_BINARY_F16(min)
DECL_CONST_OCML_BINARY_F16(maxmag)
DECL_CONST_OCML_BINARY_F16(minmag)
extern half OCML_MANGLE_F16(modf)(half, __private half *);
extern __attribute__((const)) half OCML_MANGLE_F16(nan)(ushort);
DECL_CONST_OCML_UNARY_F16(nearbyint)
DECL_CONST_OCML_BINARY_F16(nextafter)
DECL_PURE_OCML_BINARY_F16(pow)
DECL_PURE_OCML_BINARY_F16(powr)
extern __attribute__((pure)) half OCML_MANGLE_F16(pown)(half, int);
extern __attribute__((pure)) half OCML_MANGLE_F16(rootn)(half, int);
DECL_PURE_OCML_UNARY_F16(rcbrt)
DECL_CONST_OCML_BINARY_F16(remainder)
extern half OCML_MANGLE_F16(remquo)(half, half, __private int *);
DECL_CONST_OCML_UNARY_F16(rint)
DECL_CONST_OCML_UNARY_F16(round)
DECL_CONST_OCML_UNARY_F16(rsqrt)
DECL_CONST_OCML_BINARY_F16(scalb)
extern __attribute__((const)) half OCML_MANGLE_F16(scalbn)(half, int);
DECL_CONST_OCML_UNARYPRED_F16(signbit)
DECL_OCML_UNARY_F16(sin)
DECL_PURE_OCML_UNARY_F16(sinh)
DECL_OCML_UNARY_F16(sinpi)
extern half OCML_MANGLE_F16(sincos)(half, __private half *);
DECL_CONST_OCML_UNARY_F16(sqrt)
DECL_CONST_OCML_UNARY_F16(sqrt_rte)
DECL_CONST_OCML_UNARY_F16(sqrt_rtp)
DECL_CONST_OCML_UNARY_F16(sqrt_rtn)
DECL_CONST_OCML_UNARY_F16(sqrt_rtz)
DECL_OCML_UNARY_F16(tan)
DECL_CONST_OCML_UNARY_F16(tanpi)
DECL_PURE_OCML_UNARY_F16(tanh)
DECL_OCML_UNARY_F16(tgamma)
DECL_CONST_OCML_UNARY_F16(trunc)
DECL_OCML_UNARY_F16(y0)
DECL_OCML_UNARY_F16(y1)

DECL_CONST_OCML_BINARY_F16(add_rte)
DECL_CONST_OCML_BINARY_F16(add_rtp)
DECL_CONST_OCML_BINARY_F16(add_rtn)
DECL_CONST_OCML_BINARY_F16(add_rtz)

DECL_CONST_OCML_BINARY_F16(div_rte)
DECL_CONST_OCML_BINARY_F16(div_rtp)
DECL_CONST_OCML_BINARY_F16(div_rtn)
DECL_CONST_OCML_BINARY_F16(div_rtz)

DECL_CONST_OCML_TERNARY_F16(fma_rte)
DECL_CONST_OCML_TERNARY_F16(fma_rtp)
DECL_CONST_OCML_TERNARY_F16(fma_rtn)
DECL_CONST_OCML_TERNARY_F16(fma_rtz)

DECL_CONST_OCML_BINARY_F16(mul_rte)
DECL_CONST_OCML_BINARY_F16(mul_rtp)
DECL_CONST_OCML_BINARY_F16(mul_rtn)
DECL_CONST_OCML_BINARY_F16(mul_rtz)

DECL_CONST_OCML_UNARY_F16(sqrt_rte)
DECL_CONST_OCML_UNARY_F16(sqrt_rtp)
DECL_CONST_OCML_UNARY_F16(sqrt_rtn)
DECL_CONST_OCML_UNARY_F16(sqrt_rtz)

DECL_CONST_OCML_BINARY_F16(sub_rte)
DECL_CONST_OCML_BINARY_F16(sub_rtp)
DECL_CONST_OCML_BINARY_F16(sub_rtn)
DECL_CONST_OCML_BINARY_F16(sub_rtz)
#pragma OPENCL EXTENSION cl_khr_fp16 : disable

#endif // OCML_H