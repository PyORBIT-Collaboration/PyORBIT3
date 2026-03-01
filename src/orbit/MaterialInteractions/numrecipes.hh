//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    numrecipes.hh
//
// CREATED
//    10/11/2011
//
// DESCRIPTION
//   Some numerical root finding recipes
///////////////////////////////////////////////////////////////////////////
#ifndef NUMRECIPES_H
#define NUMRECIPES_H

namespace OrbitUtils{

	float fstep(float s, float r_o, float pr_o, float theta);
	float rfunc(float x, float p, float fac1);

	/**
	 * @brief Bracket subintervals where a function changes sign.
	 *
	 * Scans the interval [x1, x2] by dividing it into @p n equal subintervals
	 * and records subinterval endpoints where the function @p fx changes sign
	 * (i.e., a root is bracketed).
	 *
	 * The function @p fx must have the signature:
	 *     float fx(float x, float p1, float p2, float p3);
	 * It is evaluated at the grid points x1, x1 + dx, x1 + 2*dx, ..., x2,
	 * where dx = (x2 - x1) / n.
	 *
	 * When an interval [xb1[k], xb2[k]] is found such that fx(xb1[k]) and
	 * fx(xb2[k]) have opposite signs (fc*fp <= 0), the endpoints are stored
	 * in the output arrays xb1 and xb2 at the same index k (1-based indexing
	 * is used internally; the caller should account for this).
	 *
	 * @param fx       Pointer to the function to examine. Signature:
	 *                 float fx(float x, float param1, float param2, float param3).
	 * @param x1       Left end of the search interval.
	 * @param x2       Right end of the search interval.
	 * @param n        Number of equal subintervals to split [x1, x2] into.
	 * @param xb1      Output array to receive left endpoints of bracketing subintervals.
	 *                 Must be large enough to hold up to the capacity indicated by the
	 *                 integer referenced by @p nb (see @p nb description).
	 * @param xb2      Output array to receive right endpoints of bracketing subintervals.
	 *                 Same required capacity as @p xb1.
	 * @param nb       Reference to an integer indicating the maximum number of bracket
	 *                 intervals the caller can accept on input; on return this integer
	 *                 is set to the actual number of brackets found (nb <= original nb).
	 * @param param1   First user parameter passed through to @p fx.
	 * @param param2   Second user parameter passed through to @p fx.
	 * @param param3   Third user parameter passed through to @p fx.
	 *
	 * @return Returns 0 on completion. If the number of found brackets reaches the
	 *         input capacity (nb), the function returns immediately with xb1/xb2
	 *         filled up to that capacity and nb set equal to that capacity.
	 *
	 * @note
	 * - The function uses sequential evaluations of @p fx and treats a zero or a
	 *   sign change (fc*fp <= 0.0) as a bracket. Adjacent subintervals that both
	 *   satisfy this condition will each be reported.
	 * - The implementation uses 1-based indexing internally when filling xb1/xb2.
	 *   The caller should allocate arrays accordingly and interpret the filled
	 *   entries consistent with the calling convention used in the surrounding code.
	 */
	int zbrak(float (*fx)(float, float, float, float), float x1, float x2, int n, float *xb1, float *xb2, int &nb, float param1,
			  float param2, float param3);
	float rtbis(float (*func)(float, float, float, float), float x1, float x2, float xacc, float param1, float param2, float param3);
	float bessj0(float x);
	float bessj1(float x);
	float qsimp(float (*func)(float, float, float), float a, float b, float p, float fac1);
	float trapzd(float (*rfunc)(float, float, float), float a, float b, int n, float p, float fac1);
}

#endif
