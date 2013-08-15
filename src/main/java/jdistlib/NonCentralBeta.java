/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998   Ross Ihaka
 *  Copyright (C) 2000-9 The R Development Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */
package jdistlib;

import static java.lang.Math.*;
import static jdistlib.Constants.*;
import static jdistlib.MathFunctions.*;
import jdistlib.rng.QRandomEngine;

public class NonCentralBeta {
	public static final double density(double x, double a, double b, double ncp, boolean give_log)
	{
		final double eps = 1.e-15;

		int kMax;
		double k, ncp2, dx2, d, D;
		double sum, term, p_k, q;

		if (Double.isNaN(x) || Double.isNaN(a) || Double.isNaN(b) || Double.isNaN(ncp)) return x + a + b + ncp;
		if (ncp < 0 || a <= 0 || b <= 0) return Double.NaN;

		if (Double.isInfinite(a) || Double.isInfinite(b) || Double.isInfinite(ncp))
			return Double.NaN;

		if (x < 0 || x > 1) return(give_log ? Double.NEGATIVE_INFINITY : 0.);
		if(ncp == 0)
			return Beta.density(x, a, b, give_log);

		/* New algorithm, starting with *largest* term : */
		ncp2 = 0.5 * ncp;
		dx2 = ncp2*x;
		d = (dx2 - a - 1)/2;
		D = d*d + dx2 * (a + b) - a;
		if(D <= 0) {
			kMax = 0;
		} else {
			D = ceil(d + sqrt(D));
			kMax = (D > 0) ? (int)D : 0;
		}

		/* The starting "middle term" --- first look at it's log scale: */
		term = Beta.density(x, a + kMax, b, /* log = */ true);
		p_k = Poisson.density_raw(kMax, ncp2,true);
		if(x == 0. || Double.isInfinite(term) || Double.isInfinite(p_k)) /* if term = +Inf */
			//return R_D_exp(p_k + term);
			return (give_log ? (p_k + term) : exp(p_k + term));

		/* Now if s_k := p_k * t_k  {here = exp(p_k + term)} would underflow,
		 * we should rather scale everything and re-scale at the end:*/

		p_k += term; /* = log(p_k) + log(t_k) == log(s_k) -- used at end to rescale */
		/* mid = 1 = the rescaled value, instead of  mid = exp(p_k); */

		/* Now sum from the inside out */
		sum = term = 1. /* = mid term */;
		/* middle to the left */
		k = kMax;
		while(k > 0 && term > sum * eps) {
			k--;
			q = /* 1 / r_k = */ (k+1)*(k+a) / (k+a+b) / dx2;
			term *= q;
			sum += term;
		}
		/* middle to the right */
		term = 1.;
		k = kMax;
		do {
			q = /* r_{old k} = */ dx2 * (k+a+b) / (k+a) / (k+1);
			k++;
			term *= q;
			sum += term;
		} while (term > sum * eps);

		//return R_D_exp(p_k + log(sum));
		return (give_log ? (p_k + log(sum)) : exp(p_k + log(sum)));
	}

	public static final double cumulative_raw (double x, double o_x, double a, double b, double ncp)
	{
		/* o_x  == 1 - x  but maybe more accurate */

		/* change errmax and itrmax if desired;
		 * original (AS 226, R84) had  (errmax; itrmax) = (1e-6; 100) */
		final double errmax = 1.0e-9;
		final int    itrmax = 10000;  /* 100 is not enough for pf(ncp=200)
					     see PR#11277 */

		double a0, lbeta, c, errbd, temp, x0; //, tmp_c;
		int j; //, ierr;

		/*long*/ double ans, ax, gx, q, sumq; // TODO long double

		if (ncp < 0. || a <= 0. || b <= 0.) return Double.NaN;

		if(x < 0. || o_x > 1. || (x == 0. && o_x == 1.)) return 0.;
		if(x > 1. || o_x < 0. || (x == 1. && o_x == 0.)) return 1.;

		c = ncp / 2.;

		/* initialize the series */

		x0 = floor(max(c - 7. * sqrt(c), 0.));
		a0 = a + x0;
		lbeta = lgammafn(a0) + lgammafn(b) - lgammafn(a0 + b);
		/* temp = pbeta_raw(x, a0, b, TRUE, FALSE), but using (x, o_x): */
		double[] tt = bratio(a0, b, x, o_x, false);
		temp = tt[0]; // tmp_c = tt[1]; ierr = (int) tt[2];

		gx = exp(a0 * log(x) + b * (x < .5 ? log1p(-x) : log(o_x))
				- lbeta - log(a0));
		if (a0 > a)
			q = exp(-c + x0 * log(c) - lgammafn(x0 + 1.));
		else
			q = exp(-c);

		sumq = 1. - q;
		ans = ax = q * temp;

		/* recurse over subsequent terms until convergence is achieved */
		j = (int) x0;
		do {
			j++;
			temp -= gx;
			gx *= x * (a + b + j - 1.) / (a + j);
			q *= c / j;
			sumq -= q;
			ax = temp * q;
			ans += ax;
			errbd = (temp - gx) * sumq;
		}
		while (errbd > errmax && j < itrmax + x0);

		if (errbd > errmax) {
			// ML_ERROR(ME_PRECISION, "pnbeta");
			System.err.println("Precision error NonCentralBeta.cumulative");
		}
		if (j >= itrmax + x0) {
			//ML_ERROR(ME_NOCONV, "pnbeta");
			System.err.println("Non-convergence error NonCentralBeta.cumulative");
		}

		return ans;
	}

	static final double pnbeta2(double x, double o_x, double a, double b, double ncp, boolean lower_tail, boolean log_p)
	/* o_x  == 1 - x  but maybe more accurate */
	{
		/* long */ double ans= cumulative_raw(x, o_x, a,b, ncp); // TODO long double

		/* return R_DT_val(ans), but we want to warn about cancellation here */
		if(lower_tail) return log_p	? log(ans) : ans;
		if(ans > 1 - 1e-10) {
			// ML_ERROR(ME_PRECISION, "pnbeta");
			System.err.println("Precision error NonCentralBeta.cumulative");
		}
		ans = min(ans, 1.0);  /* Precaution */
		return log_p ? log1p(-ans) : (1 - ans);
	}

	public static final double cumulative(double x, double a, double b, double ncp, boolean lower_tail, boolean log_p)
	{
	    if (Double.isNaN(x) || Double.isNaN(a) || Double.isNaN(b) || Double.isNaN(ncp)) return x + a + b + ncp;
	    // R_P_bounds_01(x, 0., 1.);
	    if(x <= 0) return (lower_tail ? (log_p ? Double.NEGATIVE_INFINITY : 0.) : (log_p ? 0. : 1.));
	    if(x >= 1) return (lower_tail ? (log_p ? 0. : 1.) : (log_p ? Double.NEGATIVE_INFINITY : 0.));

	    return pnbeta2(x, 1-x, a, b, ncp, lower_tail, log_p);
	}

	public static final double quantile(double p, double a, double b, double ncp, boolean lower_tail, boolean log_p)
	{
		final double accu = 1e-15;
		final double Eps = 1e-14; /* must be > accu */

		double ux, lx, nx, pp;

		if (Double.isNaN(p) || Double.isNaN(a) || Double.isNaN(b) || Double.isNaN(ncp))
			return p + a + b + ncp;
		if (Double.isInfinite(a)) return Double.NaN;

		if (ncp < 0. || a <= 0. || b <= 0.) return Double.NaN;

		// R_Q_P01_boundaries(p, 0, 1);
		if (log_p) {
			if(p > 0)
				return Double.NaN;
			if(p == 0) /* upper bound*/
				return lower_tail ? 1 : 0;
			if(p == Double.NEGATIVE_INFINITY)
				return lower_tail ? 0 : 1;
		}
		else { /* !log_p */
			if(p < 0 || p > 1)
				return Double.NaN;
			if(p == 0)
				return lower_tail ? 0 : 1;
			if(p == 1)
				return lower_tail ? 1 : 0;
		}
		//p = R_DT_qIv(p);
		p = (log_p ? (lower_tail ? exp(p) : - expm1(p)) : (lower_tail ? (p) : (0.5 - (p) + 0.5)));

		/* Invert pnbeta(.) :
		 * 1. finding an upper and lower bound */
		if(p > 1 - DBL_EPSILON) return 1.0;
		pp = min(1 - DBL_EPSILON, p * (1 + Eps));
		for(ux = 0.5; ux < 1 - DBL_EPSILON && cumulative(ux, a, b, ncp, true, false) < pp; ux = 0.5*(1+ux)) ;
		pp = p * (1 - Eps);
		for(lx = 0.5; lx > Double.MIN_VALUE && cumulative(lx, a, b, ncp, true, false) > pp; lx *= 0.5) ;

		/* 2. interval (lx,ux)  halving : */
		do {
			nx = 0.5 * (lx + ux);
			if (cumulative(nx, a, b, ncp, true, false) > p) ux = nx; else lx = nx;
		} while ((ux - lx) / nx > accu);

		return 0.5 * (ux + lx);
	}

	public static final double random(double a, double b, double ncp, QRandomEngine random)
	{
		if (ncp == 0)
			return Beta.random(a, b, random);
		double x = NonCentralChiSquare.random(2 * a, ncp, random);
		x = x / (x + NonCentralChiSquare.random(2 * b, ncp, random));
		return x;
	}
}
