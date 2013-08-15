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

public class Geometric {
	public static final double density(double x, double p, boolean give_log)
	{ 
		double prob;
		if (Double.isNaN(x) || Double.isNaN(p)) return x + p;
		if (p <= 0 || p > 1) return Double.NaN;

		//R_D_nonint_check(x);
		if((abs((x) - floor((x)+0.5)) > 1e-7)) {
			//MATHLIB_WARNING("non-integer x = %f", x);
			return (give_log ? Double.NEGATIVE_INFINITY : 0.);
		}
		if (x < 0 || Double.isInfinite(x) || p == 0) return (give_log ? Double.NEGATIVE_INFINITY : 0.);
		//x = R_D_forceint(x);
		x = floor((x) + 0.5);
		/* prob = (1-p)^x, stable for small p */
		prob = Binomial.density_raw(0.,x, p,1-p, give_log);
		return((give_log) ? log(p) + prob : p*prob);
	}

	public static final double cumulative(double x, double p, boolean lower_tail, boolean log_p)
	{
		if (Double.isNaN(x) || Double.isNaN(p)) return x + p;
		if (p <= 0 || p > 1) return Double.NaN;

		if (x < 0.) return (lower_tail ? (log_p ? Double.NEGATIVE_INFINITY : 0.) : (log_p ? 0. : 1.));
		if (Double.isInfinite(x)) return (lower_tail ? (log_p ? 0. : 1.) : (log_p ? Double.NEGATIVE_INFINITY : 0.));
		x = floor(x +1e-7);

		if(p == 1.) { /* we cannot assume IEEE */
			x = lower_tail ? 1: 0;
			return log_p ? log(x) : x;
		}
		x = log1p(-p) * (x + 1);
		if (log_p)
			//return R_DT_Clog(x);
			return (lower_tail? (log_p ? ((x) > -M_LN2 ? log(-expm1(x)) : log1p(-exp(x))) : log1p(-x)) : (log_p ? (x) : log(x)));
		return lower_tail ? -expm1(x) : exp(x);
	}

	public static final double quantile(double p, double prob, boolean lower_tail, boolean log_p)
	{
		if (prob <= 0 || prob > 1) return Double.NaN;
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
		if (Double.isNaN(p) || Double.isNaN(prob)) return p + prob;
		if (prob == 1) return(0);
		/* add a fuzz to ensure left continuity */
		//return ceil(R_DT_Clog(p) / log1p(- prob) - 1 - 1e-7);
		p = (lower_tail? (log_p ? ((p) > -M_LN2 ? log(-expm1(p)) : log1p(-exp(p))) : log1p(-p)) : (log_p ? (p) : log(p)));
		return max(0, ceil(p / log1p(- prob) - 1 - 1e-12));
	}

	public static final double random(double p, QRandomEngine random)
	{
	    if (Double.isInfinite(p) || p <= 0 || p > 1) return Double.NaN;
	    return Poisson.random(Exponential.random_standard(random) * ((1 - p) / p), random);
	}
}
