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

public class Beta {

	public static final double density(double x, double a, double b, boolean log_p)
	{
	    double lval;

	    /* NaNs propagated correctly */
	    if (Double.isNaN(x) || Double.isNaN(a) || Double.isNaN(b)) return x + a + b;

	    if (a <= 0 || b <= 0) return Double.NaN;;
	    if (x < 0 || x > 1) return(log_p ? Double.NEGATIVE_INFINITY : 0.);
	    if (x == 0) {
		if(a > 1) return(log_p ? Double.NEGATIVE_INFINITY : 0.);
		if(a < 1) return(Double.POSITIVE_INFINITY);
		/* a == 1 : */ return((log_p	? log(b) : (b)));
	    }
	    if (x == 1) {
		if(b > 1) return(log_p ? Double.NEGATIVE_INFINITY : 0.);
		if(b < 1) return(Double.POSITIVE_INFINITY);
		/* b == 1 : */ return((log_p	? log(a) : (a)));
	    }
	    if (a <= 2 || b <= 2)
		lval = (a-1)*log(x) + (b-1)*log1p(-x) - lbeta(a, b);
	    else
		lval = log(a+b-1) + Binomial.density_raw(a-1, a+b-2, x, 1-x, true);

	    return (log_p ? (lval) : exp(lval));
	}

	public static final double cumulative_raw(double x, double pin, double qin, boolean lower_tail, boolean log_p)
	{
	    double x1 = 0.5 - x + 0.5, w, wc;
	    int ierr;
	    double[] temp = MathFunctions.bratio(pin, qin, x, x1, log_p); /* -> ./toms708.c */
	    w = temp[0]; wc = temp[1]; ierr = (int) temp[2];
	    /* ierr = 8 is about inaccuracy in extreme cases */
	    if(ierr > 0 && (ierr != 8 || log_p) )
	    	//MATHLIB_WARNING(_("pbeta_raw() -> bratio() gave error code %d"), ierr);
	    	return Double.NaN;
	    return lower_tail ? w : wc;
	} /* pbeta_raw() */

	public static final double cumulative(double x, double pin, double qin, boolean lower_tail, boolean log_p)
	{
	    if (Double.isNaN(x) || Double.isNaN(pin) || Double.isNaN(qin)) return x + pin + qin;
	    if (pin <= 0 || qin <= 0) return Double.NaN;;

	    if (x <= 0)	return (lower_tail ? (log_p ? Double.NEGATIVE_INFINITY : 0.) : (log_p ? 0. : 1.));
	    if (x >= 1)	return (lower_tail ? (log_p ? 0. : 1.) : (log_p ? Double.NEGATIVE_INFINITY : 0.));
	    return cumulative_raw(x, pin, qin, lower_tail, log_p);
	}

	public static final double quantile (double alpha, double p, double q, boolean lower_tail, boolean log_p)
	{
		final double fpu = 3e-308;
		/* acu_min:  Minimal value for accuracy 'acu' which will depend on (a,p); acu_min >= fpu ! */
		final double acu_min = 1e-300;
		final double lower=fpu;
		final double upper=1-2.22e-16;

		final double const1 = 2.30753;
		final double const2 = 0.27061;
		final double const3 = 0.99229;
		final double const4 = 0.04481;
		int i_pb, i_inn;
		double a, adj, logbeta, g, h, pp, p_, prev, qq, r, s, t, tx, w, y, yprev;
		double acu;
		double xinbta;
		boolean swap_tail;

		/* test for admissibility of parameters */
		if (Double.isNaN(p) || Double.isNaN(q) || Double.isNaN(alpha)) return p + q + alpha;
		if(p < 0. || q < 0.) return Double.NaN;

		//R_Q_P01_boundaries(alpha, 0, 1);
	    if (log_p) {
	    	if(alpha > 0)
	    	    return Double.NaN;
	    	if(alpha == 0) /* upper bound*/
	    	    return lower_tail ? 1 : 0;
	    	if(alpha == Double.NEGATIVE_INFINITY)
	    	    return lower_tail ? 0 : 1;
	        }
	        else { /* !log_p */
	    	if(alpha < 0 || alpha > 1)
	    	    return Double.NaN;
	    	if(alpha == 0)
	    	    return lower_tail ? 0 : 1;
	    	if(alpha == 1)
	    	    return lower_tail ? 1 : 0;
	        }

		//p_ = R_DT_qIv(alpha);/* lower_tail prob (in any case) */
	    p_ = (log_p ? (lower_tail ? exp(alpha) : - expm1(alpha)) : (lower_tail ? (alpha) : (0.5 - (alpha) + 0.5)));

		if(log_p && (p_ == 0. || p_ == 1.))
			return p_; /* better than NaN or infinite loop;
			      FIXME: suboptimal, since -Inf < alpha ! */

		/* initialize */
		logbeta = lbeta(p, q);

		/* change tail if necessary;  afterwards   0 < a <= 1/2	 */
		if (p_ <= 0.5) {
			a = p_;	pp = p; qq = q; swap_tail = false;
		} else { /* change tail, swap  p <-> q :*/
			a = (!lower_tail && !log_p)? alpha : 1 - p_;
			pp = q; qq = p; swap_tail = true;
		}

		/* calculate the initial approximation */

		/* y := {fast approximation of} qnorm(1 - a) :*/
		r = sqrt(-2 * log(a));
		y = r - (const1 + const2 * r) / (1. + (const3 + const4 * r) * r);
		if (pp > 1 && qq > 1) {
			r = (y * y - 3.) / 6.;
			s = 1. / (pp + pp - 1.);
			t = 1. / (qq + qq - 1.);
			h = 2. / (s + t);
			w = y * sqrt(h + r) / h - (t - s) * (r + 5. / 6. - 2. / (3. * h));
			xinbta = pp / (pp + qq * exp(w + w));
		} else {
			r = qq + qq;
			t = 1. / (9. * qq);
			t = r * pow(1. - t + y * sqrt(t), 3.0);
			if (t <= 0.)
				xinbta = 1. - exp((log1p(-a)+ log(qq) + logbeta) / qq);
			else {
				t = (4. * pp + r - 2.) / t;
				if (t <= 1.)
					xinbta = exp((log(a * pp) + logbeta) / pp);
				else
					xinbta = 1. - 2. / (t + 1.);
			}
		}

		/* solve for x by a modified newton-raphson method, */
		/* using the function pbeta_raw */

		r = 1 - pp;
		t = 1 - qq;
		yprev = 0.;
		adj = 1;
		/* Sometimes the approximation is negative! */
		if (xinbta < lower)
			xinbta = 0.5;
		else if (xinbta > upper)
			xinbta = 0.5;

		/* Desired accuracy should depend on  (a,p)
		 * This is from Remark .. on AS 109, adapted.
		 * However, it's not clear if this is "optimal" for IEEE double prec.

		 * acu = fmax2(acu_min, pow(10., -25. - 5./(pp * pp) - 1./(a * a)));

		 * NEW: 'acu' accuracy NOT for squared adjustment, but simple;
		 * ---- i.e.,  "new acu" = sqrt(old acu)

		 */
		acu = max(acu_min, pow(10., -13 - 2.5/(pp * pp) - 0.5/(a * a)));
		tx = prev = 0.;	/* keep -Wall happy */

		for (i_pb=0; i_pb < 1000; i_pb++) {
			y = cumulative_raw(xinbta, pp, qq, /*lower_tail = */ true, false);
			if(Double.isInfinite(y))
				return Double.NaN;

			y = (y - a) *
					exp(logbeta + r * log(xinbta) + t * log1p(-xinbta));
			if (y * yprev <= 0.)
				prev = max(abs(adj),fpu);
			g = 1;
			for (i_inn=0; i_inn < 1000;i_inn++) {
				adj = g * y;
				if (abs(adj) < prev) {
					tx = xinbta - adj; /* trial new x */
					if (tx >= 0. && tx <= 1) {
						if (prev <= acu || abs(y) <= acu) //goto L_converged;
							return swap_tail ? 1 - xinbta : xinbta;
						if (tx != 0. && tx != 1)
							break;
					}
				}
				g /= 3;
			}
			if (abs(tx - xinbta) < 1e-15*xinbta) // goto L_converged;
				return swap_tail ? 1 - xinbta : xinbta;
			xinbta = tx;
			yprev = y;
		}
		/*-- NOT converged: Iteration count --*/
		//ML_ERROR(ME_PRECISION, "qbeta");
		return Double.NaN;
	}

	public static final double random(double aa, double bb, QRandomEngine random)
	{
		final double expmax = (DBL_MAX_EXP * M_LN2);
		double a, b, alpha;
		double r, s, t, u1, u2, v, w, y, z;

		//boolean qsame;
		/* FIX-ME:  Keep Globals (properly) for threading */
		/* Uses these GLOBALS to save time when many rv's are generated : */
		// RJ's modification: The paltry saving isn't worth it since we much prefer threading
		double beta, gamma, delta, k1, k2;
		//double olda = -1.0;
		//double oldb = -1.0;

		if (aa <= 0. || bb <= 0. || (Double.isInfinite(aa) && Double.isInfinite(bb)))
			return Double.NaN;

		if (Double.isInfinite(aa))
			return 1.0;

		if (Double.isInfinite(bb))
			return 0.0;

		/* Test if we need new "initializing" */
		//qsame = (olda == aa) && (oldb == bb);
		//if (!qsame) { olda = aa; oldb = bb; }

		a = min(aa, bb);
		b = max(aa, bb); /* a <= b */
		alpha = a + b;

		//#define v_w_from__u1_bet(AA)
		//	    v = beta * log(u1 / (1.0 - u1));
		//	    if (v <= expmax) {
		//		w = AA * exp(v);
		//		if(Double.isInfinite(w)) w = Double.MAX_VALUE;
		//	    } else
		//		w = Double.MAX_VALUE;


		if (a <= 1.0) {	/* --- Algorithm BC --- */

			/* changed notation, now also a <= b (was reversed) */

			//if (!qsame) { /* initialize */
			beta = 1.0 / a;
			delta = 1.0 + b - a;
			k1 = delta * (0.0138889 + 0.0416667 * a) / (b * beta - 0.777778);
			k2 = 0.25 + (0.5 + 0.25 / delta) * a;
			//}
			/* FIXME: "do { } while()", but not trivially because of "continue"s:*/
			for(;;) {
				u1 = random.nextDouble();
				u2 = random.nextDouble();
				if (u1 < 0.5) {
					y = u1 * u2;
					z = u1 * y;
					if (0.25 * u2 + z - y >= k1)
						continue;
				} else {
					z = u1 * u1 * u2;
					if (z <= 0.25) {
						//v_w_from__u1_bet(b);
						v = beta * log(u1 / (1.0 - u1));
						if (v <= expmax) {
							w = b * exp(v);
							if(Double.isInfinite(w)) w = Double.MAX_VALUE;
						} else
							w = Double.MAX_VALUE;
						break;
					}
					if (z >= k2)
						continue;
				}

				//v_w_from__u1_bet(b);
				v = beta * log(u1 / (1.0 - u1));
				if (v <= expmax) {
					w = b * exp(v);
					if(Double.isInfinite(w)) w = Double.MAX_VALUE;
				} else
					w = Double.MAX_VALUE;

				if (alpha * (log(alpha / (a + w)) + v) - 1.3862944 >= log(z))
					break;
			}
			return (aa == a) ? a / (a + w) : w / (a + w);

		}
		else {		/* Algorithm BB */

			//if (!qsame) { /* initialize */
			beta = sqrt((alpha - 2.0) / (2.0 * a * b - alpha));
			gamma = a + 1.0 / beta;
			//}
			do {
				u1 = random.nextDouble();
				u2 = random.nextDouble();

				//v_w_from__u1_bet(a);
				v = beta * log(u1 / (1.0 - u1));
				if (v <= expmax) {
					w = a * exp(v);
					if(Double.isInfinite(w)) w = Double.MAX_VALUE;
				} else
					w = Double.MAX_VALUE;

				z = u1 * u1 * u2;
				r = gamma * v - 1.3862944;
				s = a + r - w;
				if (s + 2.609438 >= 5.0 * z)
					break;
				t = log(z);
				if (s > t)
					break;
			}
			while (r + alpha * log(alpha / (b + w)) < t);

			return (aa != a) ? b / (b + w) : w / (b + w);
		}
	}

}
