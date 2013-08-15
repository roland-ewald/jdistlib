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

public class T {
	public static final double density(double x, double n, boolean give_log)
	{
		if (Double.isNaN(x) || Double.isNaN(n))
			return x + n;

		if (n <= 0) return Double.NaN;
		if(Double.isInfinite(x))
			return (give_log ? Double.NEGATIVE_INFINITY : 0.);
		if(Double.isInfinite(n))
			return Normal.density(x, 0., 1., give_log);

		double u, ax = abs(x),
				t = -bd0(n/2.,(n+1)/2.) + stirlerr((n+1)/2.) - stirlerr(n/2.),
				x2n = x*x/n, // in  [0, Inf]
				l_x2n; // := log(sqrt(1 + x2n)) = log(1 + x2n)/2
		boolean lrg_x2n =  (x2n > 1./DBL_EPSILON);
		if (lrg_x2n) { // large x^2/n :
			l_x2n = log(ax) - log(n)/2.; // = log(x2n)/2 = 1/2 * log(x^2 / n)
			u = //  log(1 + x2n) * n/2 =  n * log(1 + x2n)/2 =
				n * l_x2n;
		}
		else if (x2n > 0.2) {
			l_x2n = log(1 + x2n)/2.;
			u = n * l_x2n;
		} else {
			l_x2n = log1p(x2n)/2.;
			u = -bd0(n/2.,(n+x*x)/2.) + x*x/2.;
		}

		if(give_log)
			return t-u - (M_LN_SQRT_2PI + l_x2n);

		// else :  if(lrg_x2n) : sqrt(1 + 1/x2n) ='= sqrt(1) = 1
		double I_sqrt_ = (lrg_x2n ? sqrt(n)/ax : exp(-l_x2n));
		return exp(t-u) * M_1_SQRT_2PI * I_sqrt_;
	}

	public static final double cumulative(double x, double n, boolean lower_tail, boolean log_p)
	{
		/* return  P[ T <= x ]	where
		 * T ~ t_{n}  (t distrib. with n degrees of freedom).

		 *	--> ./pnt.c for NON-central
		 */
		double val, nx;
		if (Double.isNaN(x) || Double.isNaN(n)) return x + n;
		if (n <= 0.0) return Double.NaN;

		if(Double.isInfinite(x))
			return (x < 0) ? (lower_tail ? (log_p ? Double.NEGATIVE_INFINITY : 0.) : (log_p ? 0. : 1.)) : (lower_tail ? (log_p ? 0. : 1.) : (log_p ? Double.NEGATIVE_INFINITY : 0.));
			if(Double.isInfinite(n))
				return Normal.cumulative(x, 0.0, 1.0, lower_tail, log_p);

			nx = 1 + (x/n)*x;
			/* FIXME: This test is probably losing rather than gaining precision,
			 * now that pbeta(*, log_p = TRUE) is much better.
			 * Note however that a version of this test *is* needed for x*x > D_MAX */
			if(nx > 1e100) { /* <==>  x*x > 1e100 * n  */
				/* Danger of underflow. So use Abramowitz & Stegun 26.5.4
		   pbeta(z, a, b) ~ z^a(1-z)^b / aB(a,b) ~ z^a / aB(a,b),
		   with z = 1/nx,  a = n/2,  b= 1/2 :
				 */
				double lval;
				lval = -0.5*n*(2*log(abs(x)) - log(n)) - lbeta(0.5*n, 0.5) - log(0.5*n);
				val = log_p ? lval : exp(lval);
			} else {
				val = (n > x * x)
					? Beta.cumulative (x * x / (n + x * x), 0.5, n / 2., false, log_p)
					: Beta.cumulative (1. / nx,             n / 2., 0.5, true, log_p);
			}
			/* Use "1 - v"  if	lower_tail  and	 x > 0 (but not both):*/
			if(x <= 0.)
				lower_tail = !lower_tail;
			if(log_p) {
				if(lower_tail) return log1p(-0.5*exp(val));
				else return val - M_LN2; /* = log(.5* pbeta(....)) */
			}
			else {
				val /= 2.;
				//return R_D_Cval(val);
				return (lower_tail ? (0.5 - (val) + 0.5) : (val));
			}
	}

	public static final double quantile(double p, double ndf, boolean lower_tail, boolean log_p)
	{
		final double eps = 1.e-12;

		double P, q;

		if (Double.isNaN(p) || Double.isNaN(ndf)) return p + ndf;

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

		if (ndf <= 0) return Double.NaN;

		if (ndf < 1) { /* based on qnt */
			final double accu = 1e-13;
			final double Eps = 1e-11; /* must be > accu */

			double ux, lx, nx, pp;

			int iter = 0;

			//p = R_DT_qIv(p);
			p = (log_p ? (lower_tail ? exp(p) : - expm1(p)) : (lower_tail ? (p) : (0.5 - (p) + 0.5)));

			/* Invert pt(.) :
			 * 1. finding an upper and lower bound */
			if(p > 1 - DBL_EPSILON) return Double.POSITIVE_INFINITY;
			pp = min(1 - DBL_EPSILON, p * (1 + Eps));
			for(ux = 1.; ux < Double.MAX_VALUE && cumulative(ux, ndf, true, false) < pp; ux *= 2);
			pp = p * (1 - Eps);
			for(lx =-1.; lx > -Double.MAX_VALUE && cumulative(lx, ndf, true, false) > pp; lx *= 2);

			/* 2. interval (lx,ux)  halving
		   regula falsi failed on qt(0.1, 0.1)
			 */
			do {
				nx = 0.5 * (lx + ux);
				if (cumulative(nx, ndf, true, false) > p) ux = nx; else lx = nx;
			} while ((ux - lx) / abs(nx) > accu && ++iter < 1000);

			if(iter >= 1000) return Double.NaN; //ML_ERROR(ME_PRECISION, "qt");

			return 0.5 * (lx + ux);
		}

		/* Old comment:
		 * FIXME: "This test should depend on  ndf  AND p  !!
		 * -----  and in fact should be replaced by
		 * something like Abramowitz & Stegun 26.7.5 (p.949)"
		 *
		 * That would say that if the qnorm value is x then
		 * the result is about x + (x^3+x)/4df + (5x^5+16x^3+3x)/96df^2
		 * The differences are tiny even if x ~ 1e5, and qnorm is not
		 * that accurate in the extreme tails.
		 */
		if (ndf > 1e20) return Normal.quantile(p, 0., 1., lower_tail, log_p);

		//P = R_D_qIv(p); /* if exp(p) underflows, we fix below */
		P = (log_p ? exp(p) : (p));

		boolean neg = (!lower_tail || P < 0.5) && (lower_tail || P > 0.5),
			is_neg_lower = (lower_tail == neg); /* both TRUE or FALSE == !xor */;
		if(neg)
			P = 2 * (log_p ? (lower_tail ? P : -expm1(p)) : (lower_tail ? (p) : (0.5 - (p) + 0.5))); //R_D_Lval(p));
		else
			P = 2 * (log_p ? (lower_tail ? -expm1(p) : P) : (lower_tail ? (0.5 - (p) + 0.5) : (p))); // R_D_Cval(p));
		/* 0 <= P <= 1 ; P = 2*min(P', 1 - P')  in all cases */

		/* Use this if(log_p) only : */
		//#define P_is_exp_2p (lower_tail == neg) /* both TRUE or FALSE == !xor */

		if (abs(ndf - 2) < eps) {	/* df ~= 2 */
			if(P > Double.MIN_VALUE) {
				if(3* P < DBL_EPSILON) /* P ~= 0 */
					q = 1 / sqrt(P);
				else if (P > 0.9)	   /* P ~= 1 */
					q = (1 - P) * sqrt(2 /(P * (2 - P)));
				else /* eps/3 <= P <= 0.9 */
					q = sqrt(2 / (P * (2 - P)) - 2);
			}
			else { /* P << 1, q = 1/sqrt(P) = ... */
				if(log_p)
					q = is_neg_lower ? exp(- p/2) / M_SQRT_2 : 1/sqrt(-expm1(p));
					else
						q = Double.POSITIVE_INFINITY;
			}
		}
		else if (ndf < 1 + eps) { /* df ~= 1  (df < 1 excluded above): Cauchy */
			if(P > 0)
				q = 1/tan(P * M_PI_2);/* == - tan((P+1) * M_PI_2) -- suffers for P ~= 0 */

			else { /* P = 0, but maybe = 2*exp(p) ! */
				if(log_p) /* 1/tan(e) ~ 1/e */
					q = is_neg_lower ? M_1_PI * exp(-p) : -1./(M_PI * expm1(p));
					else
						q = Double.POSITIVE_INFINITY;
			}
		}
		else {		/*-- usual case;  including, e.g.,  df = 1.1 */
			double x = 0., y = 0, log_P2 = 0./* -Wall */,
					a = 1 / (ndf - 0.5),
					b = 48 / (a * a),
					c = ((20700 * a / b - 98) * a - 16) * a + 96.36,
					d = ((94.5 / (b + c) - 3) / b + 1) * sqrt(a * M_PI_2) * ndf;

			boolean P_ok1 = P > Double.MIN_VALUE || !log_p,  P_ok = P_ok1;
			if(P_ok1) {
				y = pow(d * P, 2 / ndf);
				P_ok = (y >= DBL_EPSILON);
			}
			if(!P_ok) { /* log_p && P very small */
				log_P2 = is_neg_lower ? p : (p > -M_LN2 ? log(-expm1(p)) : log1p(-exp(p))); //R_Log1_Exp(p); /* == log(P / 2) */
				x = (log(d) + M_LN2 + log_P2) / ndf;
				y = exp(2 * x);
			}

			if ((ndf < 2.1 && P > 0.5) || y > 0.05 + a) { /* P > P0(df) */
				/* Asymptotic inverse expansion about normal */
				if(P_ok)
					x = Normal.quantile(0.5 * P, 0., 1., /*lower_tail*/true,  /*log_p*/false);
				else /* log_p && P underflowed */
					x = Normal.quantile(log_P2,  0., 1., lower_tail,	        /*log_p*/ true);

				y = x * x;
				if (ndf < 5)
					c += 0.3 * (ndf - 4.5) * (x + 0.6);
				c = (((0.05 * d * x - 5) * x - 7) * x - 2) * x + b + c;
				y = (((((0.4 * y + 6.3) * y + 36) * y + 94.5) / c
						- y - 3) / b + 1) * x;
				y = expm1(a * y * y);
				q = sqrt(ndf * y);
			} else if(!P_ok && x < - M_LN2 * DBL_MANT_DIG) {/* 0.5* log(DBL_EPSILON) */
				/* y above might have underflown */
				q = sqrt(ndf) * exp(-x);
			}
			else {
				y = ((1 / (((ndf + 6) / (ndf * y) - 0.089 * d - 0.822)
						* (ndf + 2) * 3) + 0.5 / (ndf + 4))
						* y - 1) * (ndf + 1) / (ndf + 2) + 1 / y;
				q = sqrt(ndf * y);
			}


			/* Now apply 2-term Taylor expansion improvement (1-term = Newton):
			 * as by Hill (1981) [ref.above] */

			/* FIXME: This can be far from optimal when log_p = TRUE
			 *      but is still needed, e.g. for qt(-2, df=1.01, log=TRUE).
			 *	Probably also improvable when  lower_tail = FALSE */

			if(P_ok1) {
				int it=0;
				while(it++ < 10 && (y = density(q, ndf, false)) > 0 &&
						!Double.isInfinite(x = (cumulative(q, ndf, false, false) - P/2) / y) &&
						abs(x) > 1e-14*abs(q))
					/* Newton (=Taylor 1 term):
					 *  q += x;
					 * Taylor 2-term : */
					q += x * (1. + x * q * (ndf + 1) / (2 * (q * q + ndf)));
			}
		}
		if(neg) q = -q;
		return q;
	}

	public static final double random(double df, QRandomEngine random)
	{
		if (Double.isNaN(df) || df <= 0.0) return Double.NaN;

		if(Double.isInfinite(df))
			return Normal.random_standard(random);
		else {
			/* Some compilers (including MW6) evaluated this from right to left
			return norm_rand() / sqrt(rchisq(df) / df); */
			double num = Normal.random_standard(random);
			return num / sqrt(ChiSquare.random(df, random) / df);
		}
	}
}
