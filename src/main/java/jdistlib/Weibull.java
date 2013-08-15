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

public class Weibull {
	public static final double density(double x, double shape, double scale, boolean give_log)
	{
		double tmp1, tmp2;
		if (Double.isNaN(x) || Double.isNaN(shape) || Double.isNaN(scale)) return x + shape + scale;
		if (shape <= 0 || scale <= 0) return Double.NaN;

		if (x < 0) return (give_log ? Double.NEGATIVE_INFINITY : 0.);
		if (Double.isInfinite(x)) return (give_log ? Double.NEGATIVE_INFINITY : 0.);
		/* need to handle x == 0 separately */
		if(x == 0 && shape < 1) return Double.POSITIVE_INFINITY;
		tmp1 = pow(x / scale, shape - 1);
		tmp2 = tmp1 * (x / scale);
		/* These are incorrect if tmp1 == 0 */
		return  give_log ?
				-tmp2 + log(shape * tmp1 / scale) :
				shape * tmp1 * exp(-tmp2) / scale;
	}

	public static final double cumulative(double x, double shape, double scale, boolean lower_tail, boolean log_p)
	{
		if (Double.isNaN(x) || Double.isNaN(shape) || Double.isNaN(scale)) return x + shape + scale;
		if (shape <= 0 || scale <= 0) return Double.NaN;

		if (x <= 0)	return (lower_tail ? (log_p ? Double.NEGATIVE_INFINITY : 0.) : (log_p ? 0. : 1.));
		x = -pow(x / scale, shape);
		if (lower_tail)
			return (log_p
					/* log(1 - exp(x))  for x < 0 : */
					//? R_Log1_Exp(x) : -expm1(x));
					? ((x) > -M_LN2 ? log(-expm1(x)) : log1p(-exp(x))) : -expm1(x));
		/* else:  !lower_tail */
		//return R_D_exp(x);
		return (log_p ? (x) : exp(x));
	}

	public static final double quantile(double p, double shape, double scale, boolean lower_tail, boolean log_p)
	{
		if (Double.isNaN(p) || Double.isNaN(shape) || Double.isNaN(scale)) return p + shape + scale;
		if (shape <= 0 || scale <= 0) return Double.NaN;

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

		//return scale * pow(- R_DT_Clog(p), 1./shape) ;
		p = (lower_tail? (log_p ? ((p) > -M_LN2 ? log(-expm1(p)) : log1p(-exp(p))) : log1p(-p)) : (log_p ? (p) : log(p)));
		return scale * pow(-p, 1./shape) ;
	}

	public static final double random(double shape, double scale, QRandomEngine random)
	{
		if (Double.isInfinite(shape) || Double.isInfinite(scale) || shape <= 0. || scale <= 0.) {
			if(scale == 0.) return 0.;
			/* else */
			return Double.NaN;
		}
		return scale * pow(-log(random.nextDouble()), 1.0 / shape);
	}
}
