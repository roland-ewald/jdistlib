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

public class NonCentralChiSquare {
	/*
	 * The density of the noncentral chi-squared distribution with "df"
	 * degrees of freedom and noncentrality parameter "ncp".
	 */
	public static final double density(double x, double df, double ncp, boolean give_log)
	{
		final double eps = 5e-15;

		double i, ncp2, q, mid, dfmid = 0, imax;
		/* long */ double sum, term; // TODO long double

		if (Double.isNaN(x) || Double.isNaN(df) || Double.isNaN(ncp)) return x + df + ncp;
		if (ncp < 0 || df <= 0) return Double.NaN;

		if (Double.isInfinite(df) || Double.isInfinite(ncp))
			return Double.NaN;

		if(x < 0) return (give_log ? Double.NEGATIVE_INFINITY : 0.);
		if(x == 0 && df < 2.)
			return Double.POSITIVE_INFINITY;
		if(ncp == 0)
			return ChiSquare.density(x, df, give_log);
		if(x == Double.POSITIVE_INFINITY) return (give_log ? Double.NEGATIVE_INFINITY : 0.);

		ncp2 = 0.5 * ncp;

		/* find max element of sum */
		imax = ceil((-(2+df) +sqrt((2-df) * (2-df) + 4 * ncp * x))/4);
		if (imax < 0) imax = 0;
		if(!Double.isInfinite(imax)) {
			dfmid = df + 2 * imax;
			mid = Poisson.density_raw(imax, ncp2, false) * ChiSquare.density(x, dfmid, false);
		} else /* imax = Inf */
			mid = 0;

		if(mid == 0) {
			/* underflow to 0 -- maybe numerically correct; maybe can be more accurate,
			 * particularly when  give_log = TRUE */
			/* Use  central-chisq approximation formula when appropriate;
			 * ((FIXME: the optimal cutoff also depends on (x,df);  use always here? )) */
			if(give_log || ncp > 1000.) {
				double nl = df + ncp, ic = nl/(nl + ncp);/* = "1/(1+b)" Abramowitz & St.*/
				return ChiSquare.density(x*ic, nl*ic, give_log);
			} else
				return (give_log ? Double.NEGATIVE_INFINITY : 0.);
		}

		sum = mid;

		/* errorbound := term * q / (1-q)  now subsumed in while() / if() below: */

		/* upper tail */
		term = mid; df = dfmid; i = imax;
		double x2 = x * ncp2;
		do {
			i++;
			q = x2 / i / df;
			df += 2;
			term *= q;
			sum += term;
		} while (q >= 1 || term * q > (1-q)*eps || term > 1e-10*sum);
		/* lower tail */
		term = mid; df = dfmid; i = imax;
		while (i != 0) {
			df -= 2;
			q = i * df / x2;
			i--;
			term *= q;
			sum += term;
			if (q < 1 && term * q <= (1-q)*eps) break;
		}
		//return R_D_val(sum);
		return (give_log ? log(sum) : (sum));
	}

	@SuppressWarnings("unused")
	public static final double cumulative_raw(double x, double f, double theta, double errmax, double reltol, int itrmax, boolean lower_tail)
	{
		double lam, x2, f2, term, bound, f_x_2n, f_2n;
		double l_lam = -1., l_x = -1.; /* initialized for -Wall */
		int n;
		boolean lamSml, tSml, is_r, is_b, is_it;
		/* long */ double ans, u, v, t, lt, lu =-1; // TODO long double

		final double _dbl_min_exp = M_LN2 * DBL_MIN_EXP;
		/*= -708.3964 for IEEE double precision */

		if (x <= 0.) {
			if(x == 0. && f == 0.)
				return lower_tail ? exp(-0.5*theta) : -expm1(-0.5*theta);
				/* x < 0  or {x==0, f > 0} */
				return lower_tail ? 0. : 1.;
		}
		if(Double.isInfinite(x))	return lower_tail ? 1. : 0.;


		if(theta < 80) { /* use 110 for Inf, as ppois(110, 80/2, lower.tail=FALSE) is 2e-20 */
			/* long */ double sum = 0, sum2 = 0, lambda = 0.5*theta, pr = exp(-lambda); // TODO long double
			//double ans;
			int i;
			/* we need to renormalize here: the result could be very close to 1 */
			for(i = 0; i < 110;  pr *= lambda/++i) {
				sum2 += pr;
				sum += pr * ChiSquare.cumulative(x, f+2*i, lower_tail, false);
				if (sum2 >= 1-1e-15) break;
			}
			ans = sum/sum2;
			return ans;
		}


		lam = .5 * theta;
		lamSml = (-lam < _dbl_min_exp);
		if(lamSml) {
			/* MATHLIB_ERROR(
		   "non centrality parameter (= %g) too large for current algorithm",
		   theta) */
			u = 0;
			lu = -lam;/* == ln(u) */
			l_lam = log(lam);
		} else {
			u = exp(-lam);
		}

		/* evaluate the first term */
		v = u;
		x2 = .5 * x;
		f2 = .5 * f;
		f_x_2n = f - x;

		if(f2 * DBL_EPSILON > 0.125 && /* very large f and x ~= f: probably needs */
				abs(t = x2 - f2) <         /* another algorithm anyway */
				sqrt(DBL_EPSILON) * f2) {
			/* evade cancellation error */
			/* t = exp((1 - t)*(2 - t/(f2 + 1))) / sqrt(2*M_PI*(f2 + 1));*/
			lt = (1 - t)*(2 - t/(f2 + 1)) - 0.5 * log(2*M_PI*(f2 + 1));
		}
		else {
			/* Usual case 2: careful not to overflow .. : */
			lt = f2*log(x2) -x2 - lgammafn(f2 + 1);
		}

		tSml = (lt < _dbl_min_exp);
		if(tSml) {
			if (x > f + theta +  5* sqrt( 2*(f + 2*theta))) {
				/* x > E[X] + 5* sigma(X) */
				return lower_tail ? 1. : 0.; /* FIXME: We could be more accurate than 0. */
			} /* else */
			l_x = log(x);
			ans = term = t = 0.;
		}
		else {
			t = exp(lt);
			ans = term = v * t;
		}

		for (n = 1, f_2n = f + 2., f_x_2n += 2.;  ; n++, f_2n += 2, f_x_2n += 2) {
			/* f_2n    === f + 2*n
			 * f_x_2n  === f - x + 2*n   > 0  <==> (f+2n)  >   x */
			if (f_x_2n > 0) {
				/* find the error bound and check for convergence */

				bound = t * x / f_x_2n;
				is_r = is_it = false;
				/* convergence only if BOTH absolute and relative error < 'bnd' */
				if (((is_b = (bound <= errmax)) &&
						(is_r = (term <= reltol * ans))) || (is_it = (n > itrmax)))
				{
					break; /* out completely */
				}

			}

			/* evaluate the next term of the */
			/* expansion and then the partial sum */

			if(lamSml) {
				lu += l_lam - log(n); /* u = u* lam / n */
				if(lu >= _dbl_min_exp) {
					/* no underflow anymore ==> change regime */
					v = u = exp(lu); /* the first non-0 'u' */
					lamSml = false;
				}
			} else {
				u *= lam / n;
				v += u;
			}
			if(tSml) {
				lt += l_x - log(f_2n);/* t <- t * (x / f2n) */
				if(lt >= _dbl_min_exp) {
					/* no underflow anymore ==> change regime */
					t = exp(lt); /* the first non-0 't' */
					tSml = false;
				}
			} else {
				t *= x / f_2n;
			}
			if(!lamSml && !tSml) {
				term = v * t;
				ans += term;
			}

		} /* for(n ...) */

		if (is_it) {
			// MATHLIB_WARNING2(_("pnchisq(x=%g, ..): not converged in %d iter."), x, itrmax);
			System.err.println("NonCentralChiSquare.density non-convergence error");
		}
		return lower_tail ? ans : 1 - ans;
	}

	/*
	 *  Algorithm AS 275 Appl.Statist. (1992), vol.41, no.2
	 *  original  (C) 1992	     Royal Statistical Society
	 *
	 *  Computes the noncentral chi-squared distribution function with
	 *  positive real degrees of freedom df and nonnegative noncentrality
	 *  parameter ncp.
	 *
	 *    Ding, C. G. (1992)
	 *    Algorithm AS275: Computing the non-central chi-squared
	 *    distribution function. Appl.Statist., 41, 478-482.
	 */
	public static final double cumulative(double x, double df, double ncp, boolean lower_tail, boolean log_p)
	{
		double ans;
		if (Double.isNaN(x) || Double.isNaN(df) || Double.isNaN(ncp)) return x + df + ncp;
		if (Double.isInfinite(df) || Double.isInfinite(ncp)) return Double.NaN;

		if (df < 0. || ncp < 0.) return Double.NaN;

		ans = cumulative_raw(x, df, ncp, 1e-12, 8*DBL_EPSILON, 1000000, lower_tail);
		if(ncp >= 80) {
			if(lower_tail) {
				ans = min(ans, 1.0);  /* e.g., pchisq(555, 1.01, ncp = 80) */
			} else { /* !lower_tail */
				/* since we computed the other tail cancellation is likely */
				if(ans < 1e-10) {
					//ML_ERROR(ME_PRECISION, "pnchisq");
					System.err.println("Precision error NonCentralChiSquare.cumulative");
				}
				ans = max(ans, 0.0);  /* Precaution PR#7099 */
			}
		}
		if (!log_p) return ans;
		/* if ans is near one, we can do better using the other tail */
		if (ncp >= 80 || ans < 1 - 1e-8) return log(ans);
		ans = cumulative_raw(x, df, ncp, 1e-12, 8*DBL_EPSILON, 1000000, !lower_tail);
		return log1p(-ans);
	}

	public static final double quantile(double p, double df, double ncp, boolean lower_tail, boolean log_p)
	{
		final double accu = 1e-13;
		final double racc = 4*DBL_EPSILON;
		/* these two are for the "search" loops, can have less accuracy: */
		final double Eps = 1e-11; /* must be > accu */
		final double rEps= 1e-10; /* relative tolerance ... */

		double ux, lx, ux0, nx, pp;

		if (Double.isNaN(p) || Double.isNaN(df) || Double.isNaN(ncp)) return p + df + ncp;
		if (Double.isInfinite(df)) return Double.NaN;

		/* Was
		 * df = floor(df + 0.5);
		 * if (df < 1 || ncp < 0) return Double.NaN;
		 */
		if (df < 0 || ncp < 0) return Double.NaN;

		//R_Q_P01_boundaries(p, 0, ML_POSINF);
		if (log_p) {
			if(p > 0)
				return Double.NaN;
			if(p == 0) /* upper bound*/
				return lower_tail ? Double.POSITIVE_INFINITY : 0;
			if(p == Double.NEGATIVE_INFINITY)
				return lower_tail ? 0 : Double.POSITIVE_INFINITY;
		}
		else { /* !log_p */
			if(p < 0 || p > 1)
				return Double.NaN;
			if(p == 0)
				return lower_tail ? 0 : Double.POSITIVE_INFINITY;
			if(p == 1)
				return lower_tail ? Double.POSITIVE_INFINITY : 0;
		}


		/* Invert pnchisq(.) :
		 * 1. finding an upper and lower bound */
		{
			/* This is Pearson's (1959) approximation,
	          which is usually good to 4 figs or so.  */
			double b, c, ff;
			b = (ncp*ncp)/(df + 3*ncp);
			c = (df + 3*ncp)/(df + 2*ncp);
			ff = (df + 2 * ncp)/(c*c);
			ux = b + c * ChiSquare.quantile(p, ff, lower_tail, log_p);
			if(ux < 0) ux = 1;
			ux0 = ux;
		}
		//p = R_D_qIv(p);
		p = (log_p ? exp(p) : (p));

		if(!lower_tail && ncp >= 80) {
			/* pnchisq is only for lower.tail = TRUE */
			if(p < 1e-10) {
				//ML_ERROR(ME_PRECISION, "qnchisq");
				System.err.println("Precision error NonCentralChiSquare.quantile");
			}
			p = 1. - p;
			lower_tail = true;
		}

		if(lower_tail) {
			if(p > 1 - DBL_EPSILON) return Double.POSITIVE_INFINITY;
			pp = min(1 - DBL_EPSILON, p * (1 + Eps));
			for(/*ux = ux*/;
					ux < Double.MAX_VALUE &&
					cumulative_raw(ux, df, ncp, Eps, rEps, 10000, true) < pp;
					ux *= 2);
			pp = p * (1 - Eps);
			for(lx = min(ux0, Double.MAX_VALUE);
					lx > Double.MIN_VALUE &&
					cumulative_raw(lx, df, ncp, Eps, rEps, 10000, true) > pp;
					lx *= 0.5);
		}
		else {
			if(p > 1 - DBL_EPSILON) return 0.0;
			pp = min(1 - DBL_EPSILON, p * (1 + Eps));
			for(/*ux = ux*/;
					ux < Double.MAX_VALUE &&
					cumulative_raw(ux, df, ncp, Eps, rEps, 10000, false) > pp;
					ux *= 2);
			pp = p * (1 - Eps);
			for(lx = min(ux0, Double.MAX_VALUE);
					lx > Double.MIN_VALUE &&
					cumulative_raw(lx, df, ncp, Eps, rEps, 10000, false) < pp;
					lx *= 0.5);
		}

		/* 2. interval (lx,ux)  halving : */
		if(lower_tail) {
			do {
				nx = 0.5 * (lx + ux);
				if (cumulative_raw(nx, df, ncp, accu, racc, 100000, true) > p)
					ux = nx;
				else
					lx = nx;
			}
			while ((ux - lx) / nx > accu);
		} else {
			do {
				nx = 0.5 * (lx + ux);
				if (cumulative_raw(nx, df, ncp, accu, racc, 100000, false) < p)
					ux = nx;
				else
					lx = nx;
			}
			while ((ux - lx) / nx > accu);
		}
		return 0.5 * (ux + lx);
	}

	/**<pre>
	 According to Hans R. Kuensch's suggestion (30 sep 2002):

	  It should be easy to do the general case (ncp > 0) by decomposing it
	  as the sum of a central chisquare with df degrees of freedom plus a
	  noncentral chisquare with zero degrees of freedom (which is a Poisson
	  mixture of central chisquares with integer degrees of freedom),
	  see Formula (29.5b-c) in Johnson, Kotz, Balakrishnan (1995).

	  The noncentral chisquare with arbitary degrees of freedom is of interest
	  for simulating the Cox-Ingersoll-Ross model for interest rates in
	  finance.

	  R code that works is

	    rchisq0 <- function(n, ncp) {
		p <- 0 < (K <- rpois(n, lambda = ncp / 2))
		r <- numeric(n)
		r[p] <- rchisq(sum(p), df = 2*K[p])
		r
	    }

	    rchisq <- function(n, df, ncp=0) {
		if(missing(ncp)) .Internal(rchisq(n, df))
		else rchisq0(n, ncp) + .Internal(rchisq(n, df))
	    }</pre>
	 */
	public static final double random(double df, double lambda, QRandomEngine random)
	{
		if (Double.isInfinite(df) || Double.isInfinite(lambda) || df < 0. || lambda < 0.)
			return Double.NaN;

		if(lambda == 0.) {
			if (df == 0.) return Double.NaN;
			return Gamma.random(df / 2., 2., random);
		}
		double r = Poisson.random( lambda / 2., random);
		if (r > 0.)  r = ChiSquare.random(2. * r, random);
		if (df > 0.) r += Gamma.random(df / 2., 2., random);
		return r;
	}
}
