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
import jdistlib.rng.QRandomEngine;

public class Logistic {
	public static final double density(double x, double location, double scale, boolean give_log)
	{
	    double e, f;
	    if (Double.isNaN(x) || Double.isNaN(location) || Double.isNaN(scale)) return x + location + scale;
	    if (scale <= 0.0) return Double.NaN;

	    x = abs((x - location) / scale);
	    e = exp(-x);
	    f = 1.0 + e;
	    return give_log ? -(x + log(scale * f * f)) : e / (scale * f * f);
	}

	/* Compute  log(1 + exp(x))  without overflow (and fast for x > 18)
	   For the two cutoffs, consider
	   curve(log1p(exp(x)) - x,       33.1, 33.5, n=2^10)
	   curve(x+exp(-x) - log1p(exp(x)), 15, 25,   n=2^11)
	*/
	private static final double log1pexp(double x) {
	    if(x <= 18.) return log1p(exp(x));
	    if(x > 33.3) return x;
	    // else: 18.0 < x <= 33.3 :
	    return x + exp(-x);
	}

	public static final double cumulative(double x, double location, double scale, boolean lower_tail, boolean log_p)
	{
		if (Double.isNaN(x) || Double.isNaN(location) || Double.isNaN(scale)) return x + location + scale;
		if (scale <= 0.0) return Double.NaN;

		x = (x - location) / scale;
		if (Double.isNaN(x)) return Double.NaN;
		//R_P_bounds_Inf_01(x);
		if(Double.isInfinite(x)) {
			if (x > 0) return (lower_tail ? (log_p ? 0. : 1.) : (log_p ? Double.NEGATIVE_INFINITY : 0.));
			/* x < 0 */return (lower_tail ? (log_p ? Double.NEGATIVE_INFINITY : 0.) : (log_p ? 0. : 1.));
		}

		return log_p ?
			-log1pexp(lower_tail ? -x : x) :
			1 / (1 + exp(lower_tail ? -x : x));
	}

	public static final double quantile(double p, double location, double scale, boolean lower_tail, boolean log_p)
	{
		if (Double.isNaN(p) || Double.isNaN(location) || Double.isNaN(scale)) return p + location + scale;
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

		if (scale <	 0.) return Double.NaN;
		if (scale == 0.) return location;

		/* p := logit(p) = log( p / (1-p) )	 : */
		if(log_p) {
			if(lower_tail)
				//p = p - R_Log1_Exp(p);
				p = p - ((p) > -M_LN2 ? log(-expm1(p)) : log1p(-exp(p)));
			else
				//p = R_Log1_Exp(p) - p;
				p = ((p) > -M_LN2 ? log(-expm1(p)) : log1p(-exp(p))) - p;
		}
		else
			p = log(lower_tail ? (p / (1. - p)) : ((1. - p) / p));

		return location + scale * p;
	}

	public static final double random(double location, double scale, QRandomEngine random)
	{
		if (Double.isNaN(location) || Double.isInfinite(scale)) return Double.NaN;
		if (scale == 0. || Double.isInfinite(location)) return location;
		double u = random.nextDouble();
		return location + scale * log(u / (1. - u));
	}
}
