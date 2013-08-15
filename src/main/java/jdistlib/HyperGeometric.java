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

public class HyperGeometric {
	public class RandomState {
		public int ks = -1;
		public int n1s = -1, n2s = -1;

		public int k, m;
		public int minjx, maxjx, n1, n2;

		public double a, d, s, w;
		public double tn, xl, xr, kl, kr, lamdl, lamdr, p1, p2, p3;
	}
	private static final HyperGeometric singleton = new HyperGeometric();
	private HyperGeometric(){}

	public static final RandomState create_random_state()
	{	return singleton.new RandomState(); }

	public static final double density(double x, double r, double b, double n, boolean give_log)
	{
	    double p, q, p1, p2, p3;

	    if (Double.isNaN(x) || Double.isNaN(r) || Double.isNaN(b) || Double.isNaN(n)) return x + r + b + n;

	    if ((r < 0. || (abs((r) - floor((r)+0.5)) > 1e-7)) ||
	    	(b < 0. || (abs((b) - floor((b)+0.5)) > 1e-7)) ||
	    	(n < 0. || (abs((n) - floor((n)+0.5)) > 1e-7)) || n > r+b) return Double.NaN;
	    //if (R_D_negInonint(x))
	    if (x < 0. || (abs((x) - floor((x)+0.5)) > 1e-7))
	    	return(give_log ? Double.NEGATIVE_INFINITY : 0.);

	    x = floor((x) + 0.5); // x = R_D_forceint(x);
	    r = floor((r) + 0.5); // r = R_D_forceint(r);
	    b = floor((b) + 0.5); // b = R_D_forceint(b);
	    n = floor((n) + 0.5); // n = R_D_forceint(n);

	    if (n < x || r < x || n - x > b) return(give_log ? Double.NEGATIVE_INFINITY : 0.);
	    if (n == 0) return((x == 0) ? (give_log ? 0. : 1.) : (give_log ? Double.NEGATIVE_INFINITY : 0.));

	    p = ((double)n)/((double)(r+b));
	    q = ((double)(r+b-n))/((double)(r+b));

	    p1 = Binomial.density_raw(x, r, p,q,give_log);
	    p2 = Binomial.density_raw(n-x,b, p,q,give_log);
	    p3 = Binomial.density_raw(n,r+b, p,q,give_log);

	    return( (give_log) ? p1 + p2 - p3 : p1*p2/p3 );
	}

	/* * Current implementation based on posting
	 * From: Morten Welinder <terra@gnome.org>
	 * Cc: R-bugs@biostat.ku.dk
	 * Subject: [Rd] phyper accuracy and efficiency (PR#6772)
	 * Date: Thu, 15 Apr 2004 18:06:37 +0200 (CEST)
	 ......
	
	 The current version has very serious cancellation issues.  For example,
	 if you ask for a small right-tail you are likely to get total cancellation.
	 For example,  phyper(59, 150, 150, 60, false, false) gives 6.372680161e-14.
	 The right answer is dhyper(0, 150, 150, 60, false) which is 5.111204798e-22.
	
	 phyper is also really slow for large arguments.
	
	 Therefore, I suggest using the code below. This is a sniplet from Gnumeric ...
	 The code isn't perfect.  In fact, if  x*(NR+NB)  is close to	n*NR,
	 then this code can take a while. Not longer than the old code, though.
	
	 -- Thanks to Ian Smith for ideas.
	*/
	public static final double pdhyper(double x, double NR, double NB, double n, boolean log_p)
	{
		/*
		 * Calculate
		 *
		 *	    phyper (x, NR, NB, n, true, false)
		 *   [log]  ----------------------------------
		 *	       dhyper (x, NR, NB, n, false)
		 *
		 * without actually calling phyper.  This assumes that
		 *
		 *     x * (NR + NB) <= n * NR
		 *
		 */
		/* long */ double sum = 0; // TODO long double
		/* long */ double term = 1; // TODO long double

		while (x > 0 && term >= DBL_EPSILON * sum) {
			term *= x * (NB - n + x) / (n + 1 - x) / (NR + 1 - x);
			sum += term;
			x--;
		}

		return log_p ? log1p(sum) : 1 + sum;
	}


	/* FIXME: The old phyper() code was basically used in ./qhyper.c as well
	 * -----  We need to sync this again!
	 */
	public static final double cumulative(double x, double NR, double NB, double n, boolean lower_tail, boolean log_p)
	{
		/* Sample of  n balls from  NR red  and	 NB black ones;	 x are red */
		double d, pd;
		if(Double.isNaN(x) || Double.isNaN(NR) || Double.isNaN(NB) || Double.isNaN(n)) return x + NR + NB + n;

		x = floor (x + 1e-7);
		NR = floor((NR) + 0.5); // NR = R_D_forceint(NR);
		NB = floor((NB) + 0.5); // NB = R_D_forceint(NB);
		n = floor((n) + 0.5);  // n  = R_D_forceint(n);

		if (NR < 0 || NB < 0 || Double.isInfinite(NR + NB) || n < 0 || n > NR + NB)
			return Double.NaN;

		if (x * (NR + NB) > n * NR) {
			/* Swap tails.	*/
			double oldNB = NB;
			NB = NR;
			NR = oldNB;
			x = n - x - 1;
			lower_tail = !lower_tail;
		}

		if (x < 0)
			return (lower_tail ? (log_p ? Double.NEGATIVE_INFINITY : 0.) : (log_p ? 0. : 1.));
		if (x >= NR || x >= n)
			return (lower_tail ? (log_p ? 0. : 1.) : (log_p ? Double.NEGATIVE_INFINITY : 0.));

		d  = density (x, NR, NB, n, log_p);
		pd = pdhyper(x, NR, NB, n, log_p);

		//return log_p ? R_DT_Log(d + pd) : R_D_Lval(d * pd);
		return log_p ? (lower_tail? (log_p	? (d + pd) : log(d + pd)) : (log_p ? ((d + pd) > -M_LN2 ? log(-expm1(d + pd)) : log1p(-exp(d + pd))) : log1p(-(d + pd))))
			: (lower_tail ? (d * pd) : (0.5 - (d * pd) + 0.5));
	}

	static final double lfastchoose(double n, double k)
	{
	    return -log(n + 1.) - lbeta(n - k + 1., k + 1.);
	}

	public static final double quantile (double p, double NR, double NB, double n, boolean lower_tail, boolean log_p)
	{
		/* This is basically the same code as  ./phyper.c  *used* to be --> FIXME! */
		double N, xstart, xend, xr, xb, sum, term;
		boolean small_N;
		if(Double.isNaN(p) || Double.isNaN(NR) || Double.isNaN(NB) || Double.isNaN(n)) return p + NR + NB + n;
		if(Double.isInfinite(p) || Double.isInfinite(NR) || Double.isInfinite(NB) || Double.isInfinite(n))
			return Double.NaN;
		NR = floor(NR + 0.5);
		NB = floor(NB + 0.5);
		N = NR + NB;
		n = floor(n + 0.5);
		if (NR < 0 || NB < 0 || n < 0 || n > N) return Double.NaN;
		/* Goal:  Find  xr (= #{red balls in sample}) such that
		 *   phyper(xr,  NR,NB, n) >= p > phyper(xr - 1,  NR,NB, n)
		 */
		xstart = max(0, n - NB);
		xend = min(n, NR);

		//R_Q_P01_boundaries(p, xstart, xend);
		if (log_p) {
			if(p > 0)
				return Double.NaN;
			if(p == 0) /* upper bound*/
				return lower_tail ? xend : xstart;
			if(p == Double.NEGATIVE_INFINITY)
				return lower_tail ? xstart : xend;
		}
		else { /* !log_p */
			if(p < 0 || p > 1)
				return Double.NaN;
			if(p == 0)
				return lower_tail ? xstart : xend;
			if(p == 1)
				return lower_tail ? xend : xstart;
		}

		xr = xstart;
		xb = n - xr;/* always ( = #{black balls in sample} ) */

		small_N = (N < 1000); /* won't have underflow in product below */
		/* if N is small,  term := product.ratio( bin.coef );
	       otherwise work with its logarithm to protect against underflow */
		term = lfastchoose(NR, xr) + lfastchoose(NB, xb) - lfastchoose(N, n);
		if(small_N) term = exp(term);
		NR -= xr;
		NB -= xb;

		if(!lower_tail || log_p) {
			//p = R_DT_qIv(p);
			p = (log_p ? (lower_tail ? exp(p) : - expm1(p)) : (lower_tail ? (p) : (0.5 - (p) + 0.5)));
		}
		p *= 1 - 1000*DBL_EPSILON; /* was 64, but failed on FreeBSD sometimes */
		sum = small_N ? term : exp(term);

		while(sum < p && xr < xend) {
			xr++;
			NB++;
			if (small_N) term *= (NR / xr) * (xb / NB);
			else term += log((NR / xr) * (xb / NB));
			sum += small_N ? term : exp(term);
			xb--;
			NR--;
		}
		return xr;
	}

	static final double afc(int i)
	{
		final double al[] =
			{
				0.0,
				0.0,/*ln(0!)=ln(1)*/
				0.0,/*ln(1!)=ln(1)*/
				0.69314718055994530941723212145817,/*ln(2) */
				1.79175946922805500081247735838070,/*ln(6) */
				3.17805383034794561964694160129705,/*ln(24)*/
				4.78749174278204599424770093452324,
				6.57925121201010099506017829290394,
				8.52516136106541430016553103634712
				/*, 10.60460290274525022841722740072165*/
			};
		double di, value;

		if (i < 0) {
			//MATHLIB_WARNING(("rhyper.c: afc(i), i=%d < 0 -- SHOULD NOT HAPPEN!\n"), i);
			return -1;/* unreached (Wall) */
		} else if (i <= 7) {
			value = al[i + 1];
		} else {
			di = i;
			value = (di + 0.5) * log(di) - di + 0.08333333333333 / di
					- 0.00277777777777 / di / di / di + 0.9189385332;
		}
		return value;
	}
	public static final double random(double nn1in, double nn2in, double kkin, QRandomEngine random)
	{	return random(nn1in, nn2in, kkin, random, null); }

	public static final double random(double nn1in, double nn2in, double kkin, QRandomEngine random, RandomState state)
	{
		final double con = 57.56462733;
		final double deltal = 0.0078;
		final double deltau = 0.0034;
		final double scale = 1e25;

		/* extern double afc(int); */

		int nn1, nn2, kk;
		int i, ix;
		boolean reject, setup1, setup2;

		double e, f, g, p, r, t, u, v, y;
		double de, dg, dr, ds, dt, gl, gu, nk, nm, ub;
		double xk, xm, xn, y1, ym, yn, yk, alv;

		/* These should become `thread_local globals' : */
		if (state == null)
			state = singleton.new RandomState();

		/* check parameter validity */

		if(Double.isInfinite(nn1in) || Double.isInfinite(nn2in) || Double.isInfinite(kkin))
			return Double.NaN;

		nn1 = (int) floor(nn1in+0.5);
		nn2 = (int) floor(nn2in+0.5);
		kk	= (int) floor(kkin +0.5);

		if (nn1 < 0 || nn2 < 0 || kk < 0 || kk > nn1 + nn2)
			return Double.NaN;

		/* if new parameter values, initialize */
		reject = true;
		if (nn1 != state.n1s || nn2 != state.n2s) {
			setup1 = true;	setup2 = true;
		} else if (kk != state.ks) {
			setup1 = false;	setup2 = true;
		} else {
			setup1 = false;	setup2 = false;
		}
		if (setup1) {
			state.n1s = nn1;
			state.n2s = nn2;
			state.tn = nn1 + nn2;
			if (nn1 <= nn2) {
				state.n1 = nn1;
				state.n2 = nn2;
			} else {
				state.n1 = nn2;
				state.n2 = nn1;
			}
		}
		if (setup2) {
			state.ks = kk;
			if (kk + kk >= state.tn) {
				state.k = (int) (state.tn - kk);
			} else {
				state.k = kk;
			}
		}
		if (setup1 || setup2) {
			state.m = (int) ((state.k + 1.0) * (state.n1 + 1.0) / (state.tn + 2.0));
			state.minjx = max(0,state.k - state.n2);
			state.maxjx = min(state.n1, state.k);
		}
		/* generate random variate --- Three basic cases */

		if (state.minjx == state.maxjx) { /* I: degenerate distribution ---------------- */
			ix = state.maxjx;
			/* return ix;
		   No, need to unmangle <TSL>*/
			/* return appropriate variate */

			if (kk + kk >= state.tn) {
				if (nn1 > nn2) {
					ix = kk - nn2 + ix;
				} else {
					ix = nn1 - ix;
				}
			} else {
				if (nn1 > nn2)
					ix = kk - ix;
			}
			return ix;

		} else if (state.m - state.minjx < 10) { /* II: inverse transformation ---------- */
			if (setup1 || setup2) {
				if (state.k < state.n2) {
					state.w = exp(con + afc(state.n2) + afc(state.n1 + state.n2 - state.k)
							- afc(state.n2 - state.k) - afc(state.n1 + state.n2));
				} else {
					state.w = exp(con + afc(state.n1) + afc(state.k)
							- afc(state.k - state.n2) - afc(state.n1 + state.n2));
				}
			}
			L10: while (true) {
				p = state.w;
				ix = state.minjx;
				u = random.nextDouble() * scale;
				L20: while (true) {
					if (u > p) {
						u -= p;
						p *= (state.n1 - ix) * (state.k - ix);
						ix++;
						p = p / ix / (state.n2 - state.k + ix);
						if (ix > state.maxjx)
							continue L10;
						continue L20;
					}
					break L10;
				}
			}
		} else { /* III : h2pe --------------------------------------------- */

			if (setup1 || setup2) {
				state.s = sqrt((state.tn - state.k) * state.k * state.n1 * state.n2 / (state.tn - 1) / state.tn / state.tn);

				/* remark: d is defined in reference without int. */
				/* the truncation centers the cell boundaries at 0.5 */

				state.d = (int) (1.5 * state.s) + .5;
				state.xl = state.m - state.d + .5;
				state.xr = state.m + state.d + .5;
				state.a = afc(state.m) + afc(state.n1 - state.m) + afc(state.k - state.m) + afc(state.n2 - state.k + state.m);
				state.kl = exp(state.a - afc((int) (state.xl)) - afc((int) (state.n1 - state.xl))
						- afc((int) (state.k - state.xl))
						- afc((int) (state.n2 - state.k + state.xl)));
				state.kr = exp(state.a - afc((int) (state.xr - 1))
						- afc((int) (state.n1 - state.xr + 1))
						- afc((int) (state.k - state.xr + 1))
						- afc((int) (state.n2 - state.k + state.xr - 1)));
				state.lamdl = -log(state.xl * (state.n2 - state.k + state.xl) / (state.n1 - state.xl + 1) / (state.k - state.xl + 1));
				state.lamdr = -log((state.n1 - state.xr + 1) * (state.k - state.xr + 1) / state.xr / (state.n2 - state.k + state.xr));
				state.p1 = state.d + state.d;
				state.p2 = state.p1 + state.kl / state.lamdl;
				state.p3 = state.p2 + state.kr / state.lamdr;
			}
			L30: while(true) {
				u = random.nextDouble() * state.p3;
				v = random.nextDouble();
				if (u < state.p1) {		/* rectangular region */
					ix = (int) (state.xl + u);
				} else if (u <= state.p2) {	/* left tail */
					ix = (int) (state.xl + log(v) / state.lamdl);
					if (ix < state.minjx)
						continue L30;
					v = v * (u - state.p1) * state.lamdl;
				} else {		/* right tail */
					ix = (int) (state.xr - log(v) / state.lamdr);
					if (ix > state.maxjx)
						continue L30;
					v = v * (u - state.p2) * state.lamdr;
				}

				/* acceptance/rejection test */

				if (state.m < 100 || ix <= 50) {
					/* explicit evaluation */
					/* The original algorithm (and TOMS 668) have
			   f = f * i * (n2 - k + i) / (n1 - i) / (k - i);
		       in the (m > ix) case, but the definition of the
		       recurrence relation on p134 shows that the +1 is
		       needed. */
					f = 1.0;
					if (state.m < ix) {
						for (i = state.m + 1; i <= ix; i++)
							f = f * (state.n1 - i + 1) * (state.k - i + 1) / (state.n2 - state.k + i) / i;
					} else if (state.m > ix) {
						for (i = ix + 1; i <= state.m; i++)
							f = f * i * (state.n2 - state.k + i) / (state.n1 - i + 1) / (state.k - i + 1);
					}
					if (v <= f) {
						reject = false;
					}
				} else {
					/* squeeze using upper and lower bounds */
					y = ix;
					y1 = y + 1.0;
					ym = y - state.m;
					yn = state.n1 - y + 1.0;
					yk = state.k - y + 1.0;
					nk = state.n2 - state.k + y1;
					r = -ym / y1;
					state.s = ym / yn;
					t = ym / yk;
					e = -ym / nk;
					g = yn * yk / (y1 * nk) - 1.0;
					dg = 1.0;
					if (g < 0.0)
						dg = 1.0 + g;
					gu = g * (1.0 + g * (-0.5 + g / 3.0));
					gl = gu - .25 * (g * g * g * g) / dg;
					xm = state.m + 0.5;
					xn = state.n1 - state.m + 0.5;
					xk = state.k - state.m + 0.5;
					nm = state.n2 - state.k + xm;
					ub = y * gu - state.m * gl + deltau
							+ xm * r * (1. + r * (-0.5 + r / 3.0))
							+ xn * state.s * (1. + state.s * (-0.5 + state.s / 3.0))
							+ xk * t * (1. + t * (-0.5 + t / 3.0))
							+ nm * e * (1. + e * (-0.5 + e / 3.0));
					/* test against upper bound */
					alv = log(v);
					if (alv > ub) {
						reject = true;
					} else {
						/* test against lower bound */
						dr = xm * (r * r * r * r);
						if (r < 0.0)
							dr /= (1.0 + r);
						ds = xn * (state.s * state.s * state.s * state.s);
						if (state.s < 0.0)
							ds /= (1.0 + state.s);
						dt = xk * (t * t * t * t);
						if (t < 0.0)
							dt /= (1.0 + t);
						de = nm * (e * e * e * e);
						if (e < 0.0)
							de /= (1.0 + e);
						if (alv < ub - 0.25 * (dr + ds + dt + de)
								+ (y + state.m) * (gl - gu) - deltal) {
							reject = false;
						}
						else {
							/* * Stirling's formula to machine accuracy
							 */
							if (alv <= (state.a - afc(ix) - afc(state.n1 - ix)
									- afc(state.k - ix) - afc(state.n2 - state.k + ix))) {
								reject = false;
							} else {
								reject = true;
							}
						}
					}
				}
				if (reject)
					continue L30;
				break L30;
			}
		}

		/* return appropriate variate */

		if (kk + kk >= state.tn) {
			if (nn1 > nn2) {
				ix = kk - nn2 + ix;
			} else {
				ix = nn1 - ix;
			}
		} else {
			if (nn1 > nn2)
				ix = kk - ix;
		}
		return ix;
	}

}
