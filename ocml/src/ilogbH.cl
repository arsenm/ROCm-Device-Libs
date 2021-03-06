/*===--------------------------------------------------------------------------
 *                   ROCm Device Libraries
 *
 * This file is distributed under the University of Illinois Open Source
 * License. See LICENSE.TXT for details.
 *===------------------------------------------------------------------------*/

#include "mathH.h"

CONSTATTR INLINEATTR int2
MATH_MANGLE2(ilogb)(half2 x)
{
    return (int2)(MATH_MANGLE(ilogb)(x.lo), MATH_MANGLE(ilogb)(x.hi));
}

CONSTATTR INLINEATTR int
MATH_MANGLE(ilogb)(half x)
{
    int r = (int)BUILTIN_FREXP_EXP_F16(x) - 1;

    if (!FINITE_ONLY_OPT()) {
        r = BUILTIN_CLASS_F16(x, CLASS_QNAN|CLASS_SNAN) ? FP_ILOGBNAN : r;
        r = BUILTIN_CLASS_F16(x, CLASS_PINF|CLASS_NINF) ? INT_MAX : r;
    }

    r = x == 0.0h ? FP_ILOGB0 : r;
    return r;
}

