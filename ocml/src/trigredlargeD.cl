/*===--------------------------------------------------------------------------
 *                   ROCm Device Libraries
 *
 * This file is distributed under the University of Illinois Open Source
 * License. See LICENSE.TXT for details.
 *===------------------------------------------------------------------------*/

#include "mathD.h"
#include "trigredD.h"

struct f64x2 {
    double hi;
    double lo;
};

static CONSTATTR INLINEATTR struct f64x2
prod2(double a, double b)
{
    double p = a * b;
    struct f64x2 r = { p, BUILTIN_FMA_F64(a, b, -p) };
    return r;
}

static CONSTATTR INLINEATTR struct f64x2
sum2(double a, double b)
{
    double s = a + b;
    double d = s - a;
    struct f64x2 r = { s, (a - (s - d)) + (b - d) };
    return r;
}

static CONSTATTR INLINEATTR struct f64x2
fsum2(double a, double b)
{
    double s = a + b;
    struct f64x2 r = { s, b - (s - a) };
    return r;
}

CONSTATTR struct redret
MATH_PRIVATE(trigredlarge)(double x)
{
    // Scale x by relevant part of 2/pi
    double p2 = BUILTIN_TRIG_PREOP_F64(x, 0);
    double p1 = BUILTIN_TRIG_PREOP_F64(x, 1);
    double p0 = BUILTIN_TRIG_PREOP_F64(x, 2);

    x = x >= 0x1.0p+945 ? BUILTIN_FLDEXP_F64(x, -128) : x;

    struct f64x2 rp0 = prod2(p0, x);
    struct f64x2 rp1 = prod2(p1, x);
    struct f64x2 rp2 = prod2(p2, x);

    struct f64x2 v2 = sum2(rp2.lo, rp1.hi);
    struct f64x2 v1 = sum2(rp1.lo, rp0.hi);
    struct f64x2 w2 = sum2(v2.lo, v1.hi);

    double e3 = rp2.hi;
    double e2 = v2.hi;
    double e1 = w2.hi;
    double e0 = w2.lo + v1.lo + rp0.lo;

    struct f64x2 e32 = fsum2(e3, e2);
    struct f64x2 e21 = fsum2(e32.lo, e1);
    struct f64x2 e10 = fsum2(e21.lo, e0);

    double f2 = e32.hi;
    double f1 = e21.hi;
    double f0 = e10.hi;

    f2 = BUILTIN_FLDEXP_F64(BUILTIN_FRACTION_F64(BUILTIN_FLDEXP_F64(f2, -2)), 2);
    f2 += f2+f1 < 0.0 ? 4.0 : 0.0;

    int i = (int)(f2 + f1);
    f2 -= (double)i;

    struct f64x2 f21 = fsum2(f2, f1);
    struct f64x2 f10 = fsum2(f21.lo, f0);
    f2 = f21.hi;
    f1 = f10.hi;
    f0 = f10.lo;

    int g = f2 >= 0.5;
    i += g;
    f2 -= g ? 1.0 : 0.0;

    struct f64x2 rf = fsum2(f2, f1);
    f2 = rf.hi;
    f1 = rf.lo;

    const double pio2h  = 0x1.921fb54442d18p+0;
    const double pio2t  = 0x1.1a62633145c07p-54;

    double rh = f2 * pio2h;
    double rt = BUILTIN_FMA_F64(f1, pio2h, BUILTIN_FMA_F64(f2, pio2t, BUILTIN_FMA_F64(f2, pio2h, -rh)));

    struct f64x2 r = fsum2(rh, rt);

    struct redret ret;
    ret.hi = r.hi;
    ret.lo = r.lo;
    ret.i = i & 0x3;
    return ret;
}

