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

public class NonCentralT {
	/**<pre>
	 *    The non-central t density is
	 *
	 *	   f(x, df, ncp) =
	 *		df^(df/2) * exp(-.5*ncp^2) /
	 *		(sqrt(pi)*gamma(df/2)*(df+x^2)^((df+1)/2)) *
	 *		sum_{k=0}^Inf  gamma((df + k + df)/2)*ncp^k /
	 *				prod(1:k)*(2*x^2/(df+x^2))^(k/2)
	 *
	 *    The functional relationship
	 *
	 *	   f(x, df, ncp) = df/x *
	 *			      (F(sqrt((df+2)/df)*x, df+2, ncp) - F(x, df, ncp))
	 *
	 *    is used to evaluate the density at x != 0 and
	 *
	 *	   f(0, df, ncp) = exp(-.5*ncp^2) /
	 *				(sqrt(pi)*sqrt(df)*gamma(df/2))*gamma((df+1)/2)
	 *
	 *    is used for x=0.
	 *
	 *    All calculations are done on log-scale to increase stability.
	 * </pre>
	 */
	public static final double density(double x, double df, double ncp, boolean give_log)
	{
		double u;
		if (Double.isNaN(x) || Double.isNaN(df)) return x + df;

		/* If non-positive df then error */
		if (df <= 0.0) return Double.NaN;

		if(ncp == 0.0) return T.density(x, df, give_log);

		/* If x is infinite then return 0 */
		if(Double.isInfinite(x)) return (give_log ? Double.NEGATIVE_INFINITY : 0.);

		/* If infinite df then the density is identical to a
	       normal distribution with mean = ncp.  However, the formula
	       loses a lot of accuracy around df=1e9
		 */
		if(Double.isInfinite(df) || df > 1e8)
			return Normal.density(x, ncp, 1., give_log);

		/* Do calculations on log scale to stabilize */

		/* Consider two cases: x ~= 0 or not */
		if (abs(x) > sqrt(df * DBL_EPSILON)) {
			u = log(df) - log(abs(x)) +
					log(abs(cumulative(x*sqrt((df+2)/df), df+2, ncp, true, false) -
							cumulative(x, df, ncp, true, false)));
			/* FIXME: the above still suffers from cancellation (but not horribly) */
		}
		else {  /* x ~= 0 : -> same value as for  x = 0 */
			u = lgammafn((df+1)/2) - lgammafn(df/2)
					- .5*(log(M_PI) + log(df) + ncp*ncp);
		}

		return (give_log ? u : exp(u));
	}

	/*  Algorithm AS 243  Lenth,R.V. (1989). Appl. Statist., Vol.38, 185-189.
	 *  ----------------
	 *  Cumulative probability at t of the non-central t-distribution
	 *  with df degrees of freedom (may be fractional) and non-centrality
	 *  parameter delta.
	 *
	 *  NOTE
	 *
	 *    Requires the following auxiliary routines:
	 *
	 *	lgammafn(x)	- log gamma function
	 *	pbeta(x, a, b)	- incomplete beta function
	 *	pnorm(x)	- normal distribution function
	 *
	 *  CONSTANTS
	 *
	 *    M_SQRT_2dPI  = 1/ {gamma(1.5) * sqrt(2)} = sqrt(2 / pi)
	 *    M_LN_SQRT_PI = ln(sqrt(pi)) = ln(pi)/2
	 */
	public static final double cumulative(double t, double df, double ncp, boolean lower_tail, boolean log_p)
	{
		double albeta, a, b, del, errbd, lambda, rxb, tt, x;
		/* long */ double geven, godd, p, q, s, tnc, xeven, xodd; // TODO long double
		int it; boolean negdel;

		/* note - itrmax and errmax may be changed to suit one's needs. */

		final int itrmax = 1000;
		final double errmax = 1.e-12;

		if (df <= 0.0) return Double.NaN;
		if(ncp == 0.0) return T.cumulative(t, df, lower_tail, log_p);

		if(Double.isInfinite(t)) {
			return (t < 0) ?
				(lower_tail ? (log_p ? Double.NEGATIVE_INFINITY : 0.) : (log_p ? 0. : 1.)) :
				(lower_tail ? (log_p ? 0. : 1.) : (log_p ? Double.NEGATIVE_INFINITY : 0.));
		}
		if (t >= 0.) {
			negdel = false; tt = t;	 del = ncp;
		}
		else {
			/* We deal quickly with left tail if extreme,
		   since pt(q, df, ncp) <= pt(0, df, ncp) = \Phi(-ncp) */
			if (ncp > 40 && (!log_p || !lower_tail)) return (lower_tail ? (log_p ? Double.NEGATIVE_INFINITY : 0.) : (log_p ? 0. : 1.));
			negdel = true;	tt = -t; del = -ncp;
		}

		if (df > 4e5 || del*del > 2*M_LN2*(-(DBL_MIN_EXP))) {
			/*-- 2nd part: if del > 37.62, then p=0 below
		  FIXME: test should depend on `df', `tt' AND `del' ! */
			/* Approx. from	 Abramowitz & Stegun 26.7.10 (p.949) */
			s = 1./(4.*df);

			return Normal.cumulative(tt*(1. - s), del, sqrt(1. + tt*tt*2.*s),
					lower_tail != negdel, log_p);
		}

		/* initialize twin series */
		/* Guenther, J. (1978). Statist. Computn. Simuln. vol.6, 199. */

		x = t * t;
		rxb = df/(x + df);/* := (1 - x) {x below} -- but more accurately */
		x = x / (x + df);/* in [0,1) */
		if (x > 0.) {/* <==>  t != 0 */
			lambda = del * del;
			p = .5 * exp(-.5 * lambda);
			if(p == 0.) { /* underflow! */

				/*========== really use an other algorithm for this case !!! */
				//ML_ERROR(ME_UNDERFLOW, "pnt");
				//ML_ERROR(ME_RANGE, "pnt"); /* |ncp| too large */
				System.err.println("Underflow error in NonCentralT.cumulative");
				return (lower_tail ? (log_p ? Double.NEGATIVE_INFINITY : 0.) : (log_p ? 0. : 1.));
			}
			q = M_SQRT_2dPI * p * del;
			s = .5 - p;
			/* s = 0.5 - p = 0.5*(1 - exp(-.5 L)) =  -0.5*expm1(-.5 L)) */
			if(s < 1e-7)
				s = -0.5 * expm1(-0.5 * lambda);
			a = .5;
			b = .5 * df;
			/* rxb = (1 - x) ^ b   [ ~= 1 - b*x for tiny x --> see 'xeven' below]
			 *       where '(1 - x)' =: rxb {accurately!} above */
			rxb = pow(rxb, b);
			albeta = M_LN_SQRT_PI + lgammafn(b) - lgammafn(.5 + b);
			xodd = Beta.cumulative(x, a, b, /*lower*/true, /*log_p*/false);
			godd = 2. * rxb * exp(a * log(x) - albeta);
			tnc = b * x;
			xeven = (tnc < DBL_EPSILON) ? tnc : 1. - rxb;
			geven = tnc * rxb;
			tnc = p * xodd + q * xeven;

			/* repeat until convergence or iteration limit */
			for(it = 1; it <= itrmax; it++) {
				a += 1.;
				xodd  -= godd;
				xeven -= geven;
				godd  *= x * (a + b - 1.) / a;
				geven *= x * (a + b - .5) / (a + .5);
				p *= lambda / (2 * it);
				q *= lambda / (2 * it + 1);
				tnc += p * xodd + q * xeven;
				s -= p;
				/* R 2.4.0 added test for rounding error here. */
				if(s < -1.e-10) { /* happens e.g. for (t,df,ncp)=(40,10,38.5), after 799 it.*/
					//ML_ERROR(ME_PRECISION, "pnt");
					System.err.println("Precision error in NonCentralT.cumulative");
					//goto finis;
					break;
				}
				if(s <= 0 && it > 1) break; // goto finis;
				errbd = 2. * s * (xodd - godd);
				if(abs(errbd) < errmax) break; //goto finis;/*convergence*/
			}
			/* non-convergence:*/
			//ML_ERROR(ME_NOCONV, "pnt");
			System.err.println("Non-convergence error in NonCentralT.cumulative");
		} else { /* x = t = 0 */
			tnc = 0.;
		}
		//finis:
		tnc += Normal.cumulative(- del, 0., 1., /*lower*/true, /*log_p*/false);

		lower_tail = lower_tail != negdel; /* xor */
		if(tnc > 1 - 1e-10 && lower_tail) {
			//ML_ERROR(ME_PRECISION, "pnt{final}");
			System.err.println("Precision error in final section of NonCentralT.cumulative");

		}

		//return R_DT_val(min(tnc, 1.) /* Precaution */);
		tnc = min(tnc, 1.);
		return (lower_tail ? (log_p ? log(tnc) : (tnc))  : (log_p	? log1p(-(tnc)) : (0.5 - (tnc) + 0.5)));
	}

	public static final double quantile(double p, double df, double ncp, boolean lower_tail, boolean log_p)
	{
		final double accu = 1e-13;
		final double Eps = 1e-11; /* must be > accu */

		double ux, lx, nx, pp;

		if (Double.isNaN(p) || Double.isNaN(df) || Double.isNaN(ncp)) return p + df + ncp;
		if (Double.isInfinite(df)) return Double.NaN;

		/* Was
		 * df = floor(df + 0.5);
		 * if (df < 1 || ncp < 0) return Double.NaN;
		 */
		if (df <= 0.0) return Double.NaN;

		if(ncp == 0.0 && df >= 1.0) return T.quantile(p, df, lower_tail, log_p);

		//R_Q_P01_boundaries(p, ML_NEGINF, ML_POSINF);
		if (log_p) {
			if(p > 0)
				return Double.NaN;
			if(p == 0) /* upper bound*/
				return lower_tail ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY;
			if(p == Double.NEGATIVE_INFINITY)
				return lower_tail ? Double.NEGATIVE_INFINITY : Double.POSITIVE_INFINITY;
		}
		else { /* !log_p */
			if(p < 0 || p > 1)
				return Double.NaN;
			if(p == 0)
				return lower_tail ? Double.NEGATIVE_INFINITY : Double.POSITIVE_INFINITY;
			if(p == 1)
				return lower_tail ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY;
		}

		//p = R_DT_qIv(p);
		p = (log_p ? (lower_tail ? exp(p) : - expm1(p)) : (lower_tail ? (p) : (0.5 - (p) + 0.5)));

		/* Invert pnt(.) :
		 * 1. finding an upper and lower bound */
		if(p > 1 - DBL_EPSILON) return Double.POSITIVE_INFINITY;
		pp = min(1 - DBL_EPSILON, p * (1 + Eps));
		for(ux = max(1., ncp);
				ux < Double.MAX_VALUE && cumulative(ux, df, ncp, true, false) < pp;
				ux *= 2);
		pp = p * (1 - Eps);
		for(lx = min(-1., -ncp);
				lx > -Double.MAX_VALUE && cumulative(lx, df, ncp, true, false) > pp;
				lx *= 2);

		/* 2. interval (lx,ux)  halving : */
		do {
			nx = 0.5 * (lx + ux);
			if (cumulative(nx, df, ncp, true, false) > p) ux = nx; else lx = nx;
		}
		while ((ux - lx) / abs(nx) > accu);

		return 0.5 * (lx + ux);
	}

	public static final double random(double df, double ncp, QRandomEngine random)
	{
		if (ncp == 0)
			return T.random(df, random);
		return Normal.random(ncp,1,random) / sqrt(ChiSquare.random(df, random)/df);
	}
}
