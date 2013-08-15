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

public class Binomial
{
	public class RandomState {
		public double c, fm, npq, p1, p2, p3, p4, qn;
		public double xl, xll, xlr, xm, xr;

		public double psave = -1.0;
		public int nsave = -1;
		public int m;
	}
	private static final Binomial singleton = new Binomial();
	private Binomial() {}

	public static final RandomState create_random_state()
	{	return singleton.new RandomState(); }

	public static final double density_raw(double x, double n, double p, double q, boolean log_p)
	{
		double lf, lc;

		if (p == 0) return((x == 0) ? (log_p ? 0. : 1.) : (log_p ? Double.NEGATIVE_INFINITY : 0.));
		if (q == 0) return((x == n) ? (log_p ? 0. : 1.) : (log_p ? Double.NEGATIVE_INFINITY : 0.));

		if (x == 0) {
			if(n == 0) return (log_p ? 0. : 1.);
			lc = (p < 0.1) ? -bd0(n,n*q) - n*p : n*log(q);
			//return( R_D_exp(lc) );
			return (log_p ? (lc) : exp(lc));
		}
		if (x == n) {
			lc = (q < 0.1) ? -bd0(n,n*p) - n*q : n*log(p);
			return (log_p ? (lc) : exp(lc));
		}
		if (x < 0 || x > n) return (log_p ? Double.NEGATIVE_INFINITY : 0.);

		/* n*p or n*q can underflow to zero if n and p or q are small.  This
	       used to occur in dbeta, and gives NaN as from R 2.3.0.  */
		lc = stirlerr(n) - stirlerr(x) - stirlerr(n-x) - bd0(x,n*p) - bd0(n-x,n*q);

		/* f = (M_2PI*x*(n-x))/n; could overflow or underflow */
		/* Upto R 2.7.1:
		 * lf = log(M_2PI) + log(x) + log(n-x) - log(n);
		 * -- following is much better for  x << n : */
		lf = log(M_2PI) + log(x) + log1p(- x/n);

		return log_p ? (lc - 0.5*lf) : exp(lc - 0.5*lf);
	}

	public static final double density(double x, double n, double p, boolean give_log)
	{
	    /* NaNs propagated correctly */
	    if (Double.isNaN(x) || Double.isNaN(n) || Double.isNaN(p)) return x + n + p;

	    //if (p < 0 || p > 1 || R_D_negInonint(n))
	    if (p < 0 || p > 1 || (n < 0. || (abs((n) - floor((n)+0.5)) > 1e-7)))
			return Double.NaN;
	    //R_D_nonint_check(x);
	    if((abs((x) - floor((x)+0.5)) > 1e-7))
	    	return (give_log ? Double.NEGATIVE_INFINITY : 0.);
	    if (x < 0 || Double.isInfinite(x)) return (give_log ? Double.NEGATIVE_INFINITY : 0.);

	    n = floor((n) + 0.5);
	    x = floor((x) + 0.5);

	    return density_raw(x, n, p, 1-p, give_log);
	}

	public static final double cumulative(double x, double n, double p, boolean lower_tail, boolean log_p)
	{
		if (Double.isNaN(x) || Double.isNaN(n) || Double.isNaN(p)) return x + n + p;
		if (Double.isInfinite(n) || Double.isInfinite(p)) return Double.NaN;
		if((abs((n) - floor((n)+0.5)) > 1e-7)) return Double.NaN;
		n = floor((n) + 0.5);
		/* PR#8560: n=0 is a valid value */
		if(n < 0 || p < 0 || p > 1) return Double.NaN;

		if (x < 0) return 0;
		x = floor(x + 1e-7);
		if (n <= x) return 1;
		return Beta.cumulative(p, x + 1, n - x, !lower_tail, log_p);
	}

	static final double do_search(double y, double[] z, double p, double n, double pr, double incr)
	{
		if(z[0] >= p) {
			/* search to the left */
			for(;;) {
				double newz = cumulative(y - incr, n, pr, /*l._t.*/true, /*log_p*/false);
				if(y == 0 || newz < p)
					return y;
				y = max(0, y - incr);
				z[0] = newz;
			}
		}
		else {		/* search to the right */
			for(;;) {
				y = min(y + incr, n);
				if(y == n || (z[0] = cumulative(y, n, pr, /*l._t.*/true, /*log_p*/false)) >= p)
					return y;
			}
		}
	}

	public static final double quantile(double p, double n, double pr, boolean lower_tail, boolean log_p)
	{
		double q, mu, sigma, gamma, z, y;

		if (Double.isNaN(p) || Double.isNaN(n) || Double.isNaN(pr)) return p + n + pr;
		if(Double.isInfinite(n) || Double.isInfinite(pr)) return Double.NaN;
		/* if log_p is true, p = -Inf is a legitimate value */
		if(Double.isInfinite(p) && !log_p) return Double.NaN;

		if(n != floor(n + 0.5)) return Double.NaN;
		if (pr < 0 || pr > 1 || n < 0) return Double.NaN;

		// R_Q_P01_boundaries(p, 0, n);
		if (log_p) {
			if(p > 0)
				return Double.NaN;
			if(p == 0) /* upper bound*/
				return lower_tail ? n : 0;
			if(p == Double.NEGATIVE_INFINITY)
				return lower_tail ? 0 : n;
		}
		else { /* !log_p */
			if(p < 0 || p > 1)
				return Double.NaN;
			if(p == 0)
				return lower_tail ? 0 : n;
			if(p == 1)
				return lower_tail ? n : 0;
		}

		if (pr == 0. || n == 0) return 0.;

		q = 1 - pr;
		if(q == 0.) return n; /* covers the full range of the distribution */
		mu = n * pr;
		sigma = sqrt(n * pr * q);
		gamma = (q - pr) / sigma;

		/* Note : "same" code in qpois.c, qbinom.c, qnbinom.c --
		 * FIXME: This is far from optimal [cancellation for p ~= 1, etc]: */
		if(!lower_tail || log_p) {
			//p = R_DT_qIv(p); /* need check again (cancellation!): */
			p = (log_p ? (lower_tail ? exp(p) : - expm1(p)) : (lower_tail ? (p) : (0.5 - (p) + 0.5)));
			if (p == 0.) return 0.;
			if (p == 1.) return n;
		}
		/* temporary hack --- FIXME --- */
		if (p + 1.01*DBL_EPSILON >= 1.) return n;

		/* y := approx.value (Cornish-Fisher expansion) :  */
		z = Normal.quantile(p, 0., 1., /*lower_tail*/true, /*log_p*/false);
		y = floor(mu + sigma * (z + gamma * (z*z - 1) / 6) + 0.5);

		if(y > n) /* way off */ y = n;

		z = cumulative(y, n, pr, /*lower_tail*/true, /*log_p*/false);

		/* fuzz to ensure left continuity: */
		p *= 1 - 64*DBL_EPSILON;

		double[] zp = new double[] {z};
		if(n < 1e5) return do_search(y, zp, p, n, pr, 1);
		/* Otherwise be a bit cleverer in the search */
		{
			double incr = floor(n * 0.001), oldincr;
			do {
				oldincr = incr;
				y = do_search(y, zp, p, n, pr, incr);
				incr = max(1, floor(incr/100));
			} while(oldincr > 1 && incr > n*1e-15);
			return y;
		}
	}

	public static final double random(double nin, double pp, QRandomEngine random, RandomState state)
	{
		/* FIXME: These should become THREAD_specific globals : */
		if (state == null) state = singleton.new RandomState();

		double f, f1, f2, u, v, w, w2, x, x1, x2, z, z2;
		double p, q, np, g, r, al, alv, amaxp, ffm, ynorm;
		int i,ix,k, n;

		if (Double.isInfinite(nin)) return Double.NaN;
		r = floor(nin + 0.5);
		if (r != nin) return Double.NaN;
		if (Double.isInfinite(pp) ||
				/* n=0, p=0, p=1 are not errors <TSL>*/
				r < 0 || pp < 0. || pp > 1.)	return Double.NaN;

		if (r == 0 || pp == 0.) return 0;
		if (pp == 1.) return r;

		if (r >= Integer.MAX_VALUE)/* evade integer overflow, and r == INT_MAX gave only even values */
			return quantile(random.nextDouble(), r, pp, /*lower_tail*/ false, /*log_p*/ false);
		/* else */
		n = (int) r;

		p = min(pp, 1. - pp);
		q = 1. - p;
		np = n * p;
		r = p / q;
		g = r * (n + 1);

		/* Setup, perform only when parameters change [using static (globals): */

		/* FIXING: Want this thread safe
	       -- use as little (thread globals) as possible
		 */
		if (pp != state.psave || n != state.nsave) {
			state.psave = pp;
			state.nsave = n;
			if (np < 30.0) {
				/* inverse cdf logic for mean less than 30 */
				state.qn = pow(q, (double) n);
				// goto L_np_small;
				//L_np_small:
				/*---------------------- np = n*p < 30 : ------------------------- */
				for(;;) {
					ix = 0;
					f = state.qn;
					u = random.nextDouble();
					for(;;) {
						if (u < f) {
							//goto finis;
							if (state.psave > 0.5) ix = n - ix;
							return (double)ix;
						}
						if (ix > 110)
							break;
						u -= f;
						ix++;
						f *= (g / ix - r);
					}
				}
			} else {
				ffm = np + p;
				state.m = (int) ffm;
				state.fm = state.m;
				state.npq = np * q;
				state.p1 = (int)(2.195 * sqrt(state.npq) - 4.6 * q) + 0.5;
				state.xm = state.fm + 0.5;
				state.xl = state.xm - state.p1;
				state.xr = state.xm + state.p1;
				state.c = 0.134 + 20.5 / (15.3 + state.fm);
				al = (ffm - state.xl) / (ffm - state.xl * p);
				state.xll = al * (1.0 + 0.5 * al);
				al = (state.xr - ffm) / (state.xr * q);
				state.xlr = al * (1.0 + 0.5 * al);
				state.p2 = state.p1 * (1.0 + state.c + state.c);
				state.p3 = state.p2 + state.c / state.xll;
				state.p4 = state.p3 + state.c / state.xlr;
			}
		} else if (n == state.nsave) {
			if (np < 30.0) {
				// goto L_np_small;
				//L_np_small:
				/*---------------------- np = n*p < 30 : ------------------------- */
				for(;;) {
					ix = 0;
					f = state.qn;
					u = random.nextDouble();
					for(;;) {
						if (u < f) {
							//goto finis;
							if (state.psave > 0.5) ix = n - ix;
							return (double)ix;
						}
						if (ix > 110)
							break;
						u -= f;
						ix++;
						f *= (g / ix - r);
					}
				}
			}
		}

		/*-------------------------- np = n*p >= 30 : ------------------- */
		for(;;) {
			u = random.nextDouble() * state.p4;
			v = random.nextDouble();
			/* triangular region */
			if (u <= state.p1) {
				ix = (int) (state.xm - state.p1 * v + u);
				// goto finis;
				if (state.psave > 0.5) ix = n - ix;
				return (double)ix;
			}
			/* parallelogram region */
			if (u <= state.p2) {
				x = state.xl + (u - state.p1) / state.c;
				v = v * state.c + 1.0 - abs(state.xm - x) / state.p1;
				if (v > 1.0 || v <= 0.)
					continue;
				ix = (int) x;
			} else {
				if (u > state.p3) {	/* right tail */
					ix = (int) (state.xr - log(v) / state.xlr);
					if (ix > n)
						continue;
					v = v * (u - state.p3) * state.xlr;
				} else {/* left tail */
					ix = (int) (state.xl + log(v) / state.xll);
					if (ix < 0)
						continue;
					v = v * (u - state.p2) * state.xll;
				}
			}
			/* determine appropriate way to perform accept/reject test */
			k = abs(ix - state.m);
			if (k <= 20 || k >= state.npq / 2 - 1) {
				/* explicit evaluation */
				f = 1.0;
				if (state.m < ix) {
					for (i = state.m + 1; i <= ix; i++)
						f *= (g / i - r);
				} else if (state.m != ix) {
					for (i = ix + 1; i <= state.m; i++)
						f /= (g / i - r);
				}
				if (v <= f) {
					// goto finis;
					if (state.psave > 0.5) ix = n - ix;
					return (double)ix;
				}
			} else {
				/* squeezing using upper and lower bounds on log(f(x)) */
				amaxp = (k / state.npq) * ((k * (k / 3. + 0.625) + 0.1666666666666) / state.npq + 0.5);
				ynorm = -k * k / (2.0 * state.npq);
				alv = log(v);
				if (alv < ynorm - amaxp) {
					//goto finis;
					if (state.psave > 0.5) ix = n - ix;
					return (double)ix;
				}
				if (alv <= ynorm + amaxp) {
					/* stirling's formula to machine accuracy */
					/* for the final acceptance/rejection test */
					x1 = ix + 1;
					f1 = state.fm + 1.0;
					z = n + 1 - state.fm;
					w = n - ix + 1.0;
					z2 = z * z;
					x2 = x1 * x1;
					f2 = f1 * f1;
					w2 = w * w;
					if (alv <= state.xm * log(f1 / x1) + (n - state.m + 0.5) * log(z / w) + (ix - state.m) * log(w * p / (x1 * q)) + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / f2) / f2) / f2) / f2) / f1 / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / z2) / z2) / z2) / z2) / z / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / x2) / x2) / x2) / x2) / x1 / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / w2) / w2) / w2) / w2) / w / 166320.)
					{
						//goto finis;
						if (state.psave > 0.5) ix = n - ix;
						return (double)ix;
					}
				}
			}
		}
	}

}
