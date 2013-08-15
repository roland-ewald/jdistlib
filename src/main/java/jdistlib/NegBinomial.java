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

public class NegBinomial {
	public static final double density(double x, double size, double prob, boolean give_log)
	{
		double ans, p;
		if (Double.isNaN(x) || Double.isNaN(size) || Double.isNaN(prob)) return x + size + prob;

		if (prob <= 0 || prob > 1 || size < 0) return Double.NaN;
		//R_D_nonint_check(x);
		if((abs((x) - floor((x)+0.5)) > 1e-7)) {
			//MATHLIB_WARNING("non-integer x = %f", x);
			return (give_log ? Double.NEGATIVE_INFINITY : 0.);
		}

		if (x < 0 || Double.isInfinite(x)) return (give_log ? Double.NEGATIVE_INFINITY : 0.);
		//x = R_D_forceint(x);
		x = floor((x) + 0.5);

		ans = Binomial.density_raw(size, x+size, prob, 1-prob, give_log);
		p = ((double)size)/(size+x);
		return((give_log) ? log(p) + ans : p * ans);
	}

	public static final double density_mu(double x, double size, double mu, boolean give_log)
	{
		/* originally, just set  prob :=  size / (size + mu)  and called dbinom_raw(),
		 * but that suffers from cancellation when   mu << size  */
		double ans, p;

		if (Double.isNaN(x) || Double.isNaN(size) || Double.isNaN(size)) return x + size + mu;

		if (mu < 0 || size < 0) return Double.NaN;
		// R_D_nonint_check(x);
		if((abs((x) - floor((x)+0.5)) > 1e-7)) {
			//MATHLIB_WARNING("non-integer x = %f", x);
			return (give_log ? Double.NEGATIVE_INFINITY : 0.);
		}

		if (x < 0 || Double.isInfinite(x)) return (give_log ? Double.NEGATIVE_INFINITY : 0.);
		//x = R_D_forceint(x);
		x = floor((x) + 0.5);
		if(x == 0) { /* be accurate, both for n << mu, and n >> mu :*/
			x = size * (size < mu ? log(size/(size+mu)) : log1p(- mu/(size+mu)));
			return (give_log ? (x) : exp(x));
		}
		if(x < 1e-10 * size) { /* don't use dbinom_raw() but MM's formula: */
			/* FIXME --- 1e-8 shows problem; rather use algdiv() from ./toms708.c */
			x = x * log(size*mu / (size+mu)) - mu - lgammafn(x+1) + log1p(x*(x-1)/(2*size));
			return (give_log ? (x) : exp(x));
		}
		/* else: no unnecessary cancellation inside dbinom_raw, when
		 * x_ = size and n_ = x+size are so close that n_ - x_ loses accuracy
		 */
		ans = Binomial.density_raw(size, x+size, size/(size+mu), mu/(size+mu), give_log);
		p = ((double)size)/(size+x);
		return((give_log) ? log(p) + ans : p * ans);
	}

	public static final double cumulative(double x, double size, double prob, boolean lower_tail, boolean log_p)
	{
		if (Double.isNaN(x) || Double.isNaN(size) || Double.isNaN(prob)) return x + size + prob;
		if(Double.isInfinite(size) || Double.isInfinite(prob)) return Double.NaN;
		if (size <= 0 || prob <= 0 || prob > 1)	return Double.NaN;

		if (x < 0) return (lower_tail ? (log_p ? Double.NEGATIVE_INFINITY : 0.) : (log_p ? 0. : 1.));
		if (Double.isInfinite(x)) return (lower_tail ? (log_p ? 0. : 1.) : (log_p ? Double.NEGATIVE_INFINITY : 0.));
		x = floor(x + 1e-7);
		return Beta.cumulative(prob, size, x + 1, lower_tail, log_p);
	}

	public static final double cumulative_mu(double x, double size, double mu, boolean lower_tail, boolean log_p)
	{
		if (Double.isNaN(x) || Double.isNaN(size) || Double.isNaN(size)) return x + size + mu;
		if(Double.isInfinite(size) || Double.isInfinite(mu)) return Double.NaN;
		if (size <= 0 || mu < 0) return Double.NaN;

		if (x < 0) return (lower_tail ? (log_p ? Double.NEGATIVE_INFINITY : 0.) : (log_p ? 0. : 1.));
		if (Double.isInfinite(x)) return (lower_tail ? (log_p ? 0. : 1.) : (log_p ? Double.NEGATIVE_INFINITY : 0.));
		x = floor(x + 1e-7);
		/* return
		 * pbeta(pr, size, x + 1, lower_tail, log_p);  pr = size/(size + mu), 1-pr = mu/(size+mu)
		 *
		 *= pbeta_raw(pr, size, x + 1, lower_tail, log_p)
		 *            x.  pin   qin
		 *=  bratio (pin,  qin, x., 1-x., &w, &wc, &ierr, log_p),  and return w or wc ..
		 *=  bratio (size, x+1, pr, 1-pr, &w, &wc, &ierr, log_p) */
		{
			double w, wc;
			double[] temp = bratio(size, x+1, size/(size+mu), mu/(size+mu), log_p);
			w = temp[0]; wc = temp[1];
			//if(temp[2] > 0)
			//	MATHLIB_WARNING(_("pnbinom_mu() -> bratio() gave error code %d"), ierr);
			return lower_tail ? w : wc;
		}
	}

	static final double do_search(double y, double []z, double p, double n, double pr, double incr)
	{
		if(z[0] >= p) {
			/* search to the left */
			for(;;) {
				if(y == 0 ||
						(z[0] = cumulative(y - incr, n, pr, /*l._t.*/true, /*log_p*/false)) < p)
					return y;
				y = max(0, y - incr);
			}
		}
		else {		/* search to the right */
			for(;;) {
				y = y + incr;
				if((z[0] = cumulative(y, n, pr, /*l._t.*/true, /*log_p*/false)) >= p)
					return y;
			}
		}
	}

	public static final double quantile(double p, double size, double prob, boolean lower_tail, boolean log_p)
	{
		double P, Q, mu, sigma, gamma, z, y;

		if (Double.isNaN(p) || Double.isNaN(size) || Double.isNaN(prob)) return p + size + prob;
		if (prob <= 0 || prob > 1 || size <= 0) return Double.NaN;
		/* FIXME: size = 0 is well defined ! */
		if (prob == 1) return 0;

		// R_Q_P01_boundaries(p, 0, ML_POSINF);
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

		Q = 1.0 / prob;
		P = (1.0 - prob) * Q;
		mu = size * P;
		sigma = sqrt(size * P * Q);
		gamma = (Q + P)/sigma;

		/* Note : "same" code in qpois.c, qbinom.c, qnbinom.c --
		 * FIXME: This is far from optimal [cancellation for p ~= 1, etc]: */
		if(!lower_tail || log_p) {
			//p = R_DT_qIv(p); /* need check again (cancellation!): */
			p = (log_p ? (lower_tail ? exp(p) : - expm1(p)) : (lower_tail ? (p) : (0.5 - (p) + 0.5)));
			if (p == (lower_tail ? (log_p ? Double.NEGATIVE_INFINITY : 0.) : (log_p ? 0. : 1.))) return 0;
			if (p == (lower_tail ? (log_p ? 0. : 1.) : (log_p ? Double.NEGATIVE_INFINITY : 0.))) return Double.POSITIVE_INFINITY;
		}
		/* temporary hack --- FIXME --- */
		if (p + 1.01*DBL_EPSILON >= 1.) return Double.POSITIVE_INFINITY;

		/* y := approx.value (Cornish-Fisher expansion) :  */
		z = Normal.quantile(p, 0., 1., /*lower_tail*/true, /*log_p*/false);
		y = floor(mu + sigma * (z + gamma * (z*z - 1) / 6) + 0.5);

		z = cumulative(y, size, prob, /*lower_tail*/true, /*log_p*/false);

		/* fuzz to ensure left continuity: */
		p *= 1 - 64*DBL_EPSILON;

		/* If the C-F value is not too large a simple search is OK */
		double[] zp = new double[] {z};
		if(y < 1e5) return do_search(y, zp, p, size, prob, 1);
		/* Otherwise be a bit cleverer in the search */
		{
			double incr = floor(y * 0.001), oldincr;
			do {
				oldincr = incr;
				y = do_search(y, zp, p, size, prob, incr);
				incr = max(1, floor(incr/100));
			} while(oldincr > 1 && incr > y*1e-15);
			return y;
		}
	}

	public static final double quantile_mu(double p, double size, double mu, boolean lower_tail, boolean log_p)
	{
		/* FIXME!  Implement properly!! (not losing accuracy for very large size (prob ~= 1)*/
		return quantile(p, size, /* prob = */ size/(size+mu), lower_tail, log_p);
	}

	public static final double random(double size, double prob, QRandomEngine random)
	{
	    if(Double.isInfinite(size) || Double.isInfinite(prob) || size <= 0 || prob <= 0 || prob > 1)
	    	/* prob = 1 is ok, PR#1218 */
	    	return Double.NaN;
	    return (prob == 1) ? 0 : Poisson.random(Gamma.random(size, (1 - prob) / prob, random), random);
	}
}
