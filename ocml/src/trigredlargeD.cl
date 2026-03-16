/*===--------------------------------------------------------------------------
 *                   ROCm Device Libraries
 *
 * This file is distributed under the University of Illinois Open Source
 * License. See LICENSE.TXT for details.
 *===------------------------------------------------------------------------*/

#include "mathD.h"
#include "trigredD.h"

CONSTATTR struct redret
MATH_PRIVATE(trigredlarge)(double x)
{
    // Scale x by relevant part of 2/pi
    double p2 = BUILTIN_TRIG_PREOP_F64(x, 0);
    double p1 = BUILTIN_TRIG_PREOP_F64(x, 1);
    double p0 = BUILTIN_TRIG_PREOP_F64(x, 2);

    x = x >= 0x1.0p+945 ? BUILTIN_FLDEXP_F64(x, -128) : x;

    // PROD2(p0, x, p0h, p0l)
    double p0h = p0 * x;
    double p0l = BUILTIN_FMA_F64(p0, x, -p0h);

    // PROD2(p1, x, p1h, p1l)
    double p1h = p1 * x;
    double p1l = BUILTIN_FMA_F64(p1, x, -p1h);

    // PROD2(p2, x, p2h, p2l)
    double p2h = p2 * x;
    double p2l = BUILTIN_FMA_F64(p2, x, -p2h);

    // SUM2(p2l, p1h, v2h, v2l)
    double v2h = p2l + p1h;
    double v2h_aa = v2h - p1h;
    double v2h_bb = v2h - v2h_aa;
    double v2l = (p2l - v2h_aa) + (p1h - v2h_bb);

    // SUM2(p1l, p0h, v1h, v1l)
    double v1h = p1l + p0h;
    double v1h_aa = v1h - p0h;
    double v1h_bb = v1h - v1h_aa;
    double v1l = (p1l - v1h_aa) + (p0h - v1h_bb);

    // SUM2(v2l, v1h, w2h, w2l)
    double w2h = v2l + v1h;
    double w2h_aa = w2h - v1h;
    double w2h_bb = w2h - w2h_aa;
    double w2l = (v2l - w2h_aa) + (v1h - w2h_bb);

    double e3 = p2h;
    double e2 = v2h;
    double e1 = w2h;
    double e0 = w2l + v1l + p0l;

    // FSUM2(e3, e2, e3, e2)
    double e3_e2_sum = e3 + e2;
    e2 = e2 - (e3_e2_sum - e3);
    e3 = e3_e2_sum;

    // FSUM2(e2, e1, e2, e1)
    double e2_e1_sum = e2 + e1;
    e1 = e1 - (e2_e1_sum - e2);
    e2 = e2_e1_sum;

    // FSUM2(e1, e0, e1, e0)
    double e1_e0_sum = e1 + e0;
    e0 = e0 - (e1_e0_sum - e1);
    e1 = e1_e0_sum;

    double f2 = e3;
    double f1 = e2;
    double f0 = e1;

    f2 = BUILTIN_FLDEXP_F64(BUILTIN_FRACTION_F64(BUILTIN_FLDEXP_F64(f2, -2)), 2);
    f2 += f2+f1 < 0.0 ? 4.0 : 0.0;

    int i = (int)(f2 + f1);
    f2 -= (double)i;

    // FSUM2(f2, f1, f2, f1)
    double f2_f1_sum = f2 + f1;
    f1 = f1 - (f2_f1_sum - f2);
    f2 = f2_f1_sum;

    // FSUM2(f1, f0, f1, f0)
    double f1_f0_sum = f1 + f0;
    f0 = f0 - (f1_f0_sum - f1);
    f1 = f1_f0_sum;

    int g = f2 >= 0.5;
    i += g;
    f2 -= g ? 1.0 : 0.0;

    // FSUM2(f2, f1, f2, f1)
    double f2_f1_sum2 = f2 + f1;
    f1 = f1 - (f2_f1_sum2 - f2);
    f2 = f2_f1_sum2;

    const double pio2h  = 0x1.921fb54442d18p+0;
    const double pio2t  = 0x1.1a62633145c07p-54;

    double rh = f2 * pio2h;
    double rt = BUILTIN_FMA_F64(f1, pio2h, BUILTIN_FMA_F64(f2, pio2t, BUILTIN_FMA_F64(f2, pio2h, -rh)));

    // FSUM2(rh, rt, rh, rt)
    double rh_rt_sum = rh + rt;
    rt = rt - (rh_rt_sum - rh);
    rh = rh_rt_sum;

    struct redret ret;
    ret.hi = rh;
    ret.lo = rt;
    ret.i = i & 0x3;
    return ret;
}

