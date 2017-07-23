/*===--------------------------------------------------------------------------
 *                   ROCm Device Libraries
 *
 * This file is distributed under the University of Illinois Open Source
 * License. See LICENSE.TXT for details.
 *===------------------------------------------------------------------------*/

#include "mathF.h"

ALWAYS_INLINE CONSTATTR float
MATH_MANGLE(nearbyint)(float x)
{
    return BUILTIN_RINT_F32(x);
}

