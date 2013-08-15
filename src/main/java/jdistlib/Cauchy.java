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

public class Cauchy {
	public static final double density(double x, double location, double scale, boolean give_log)
	{
		double y;
		/* NaNs propagated correctly */
		if (Double.isNaN(x) || Double.isNaN(location) || Double.isNaN(scale)) return x + location + scale;
		if (scale <= 0) return Double.NaN;

		y = (x - location) / scale;
		return give_log ?
			- log(M_PI * scale * (1. + y * y)) :
			1. / (M_PI * scale * (1. + y * y));
	}

	public static final double cumulative(double x, double location, double scale, boolean lower_tail, boolean log_p)
	{
		if (Double.isNaN(x) || Double.isNaN(location) || Double.isNaN(scale)) return x + location + scale;
		if (scale <= 0) return Double.NaN;

		x = (x - location) / scale;
		if (Double.isNaN(x)) return Double.NaN;
		if(Double.isInfinite(x)) {
			if(x < 0) return (lower_tail ? (log_p ? Double.NEGATIVE_INFINITY : 0.) : (log_p ? 0. : 1.));
			else return (lower_tail ? (log_p ? 0. : 1.) : (log_p ? Double.NEGATIVE_INFINITY : 0.));
		}
		if (!lower_tail)
			x = -x;
		/* for large x, the standard formula suffers from cancellation.
		 * This is from Morten Welinder thanks to  Ian Smith's  atan(1/x) : */
		if (abs(x) > 1) {
			double y = atan(1/x) / M_PI;
			//return (x > 0) ? R_D_Clog(y) : R_D_val(-y);
			return (x > 0) ? (log_p	? log1p(-(y)) : (0.5 - (y) + 0.5)) : (log_p	? log(-y) : (-y));
		}
		//return R_D_val(0.5 + atan(x) / M_PI);
		x = 0.5 + atan(x) / M_PI;
		return (log_p ? log(x) : (x));
	}

	public static final double quantile(double p, double location, double scale, boolean lower_tail, boolean log_p)
	{
		if (Double.isNaN(p) || Double.isNaN(location) || Double.isNaN(scale))
			return p + location + scale;
		//R_Q_P01_check(p);
		if ((log_p	&& p > 0) || (!log_p && (p < 0 || p > 1)) ) return Double.NaN;
		if (scale <= 0 || Double.isInfinite(scale)) {
			if (scale == 0) return location;
			/* else */ return Double.NaN;
		}

		if (log_p) {
			if (p > -1) {
				/* when ep := exp(p),
				 * tan(pi*ep)= -tan(pi*(-ep))= -tan(pi*(-ep)+pi) = -tan(pi*(1-ep)) =
				 *		 = -tan(pi*(-expm1(p))
				 * for p ~ 0, exp(p) ~ 1, tan(~0) may be better than tan(~pi).
				 */
				if (p == 0.) /* needed, since 1/tan(-0) = -Inf  for some arch. */
					return location + (lower_tail ? scale : -scale) * Double.POSITIVE_INFINITY;
				lower_tail = !lower_tail;
				p = -expm1(p);
			} else
				p = exp(p);
		} else if (p == 1.)
			return location + (lower_tail ? scale : -scale) * Double.POSITIVE_INFINITY;

		return location + (lower_tail ? -scale : scale) / tan(M_PI * p);
		/*	-1/tan(pi * p) = -cot(pi * p) = tan(pi * (p - 1/2))  */
	}

	public static final double random(double location, double scale, QRandomEngine random)
	{
		if (Double.isNaN(location) || Double.isInfinite(scale) || scale < 0) return Double.NaN;
		if (scale == 0. || Double.isInfinite(location)) return location;
		return location + scale * tan(M_PI * random.nextDouble());
	}
}
