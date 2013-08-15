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
import jdistlib.rng.QRandomEngine;

public class Uniform {
	public static final double density(double x, double a, double b, boolean give_log)
	{
		if (Double.isNaN(x) || Double.isNaN(a) || Double.isNaN(b)) return x + a + b;
		if (b <= a) return Double.NaN;

		if (a <= x && x <= b) return give_log ? -log(b - a) : 1. / (b - a);
		return (give_log ? Double.NEGATIVE_INFINITY : 0.);
	}

	public static final double cumulative(double x, double a, double b, boolean lower_tail, boolean log_p)
	{
		if (Double.isNaN(x) || Double.isNaN(a) || Double.isNaN(b)) return x + a + b;
		if (b < a || Double.isInfinite(a) || Double.isInfinite(b)) return Double.NaN;
		if (x >= b) return (lower_tail ? (log_p ? 0. : 1.) : (log_p ? Double.NEGATIVE_INFINITY : 0.));
		if (x <= a) return (lower_tail ? (log_p ? Double.NEGATIVE_INFINITY : 0.) : (log_p ? 0. : 1.));
		if (lower_tail) {
			x = (x - a) / (b - a);
			return (log_p ? log(x) : (x));
		}
		x = (b - x) / (b - a);
		return (log_p ? log(x) : (x));
	}

	public static final double quantile(double p, double a, double b, boolean lower_tail, boolean log_p)
	{
		if (Double.isNaN(p) || Double.isNaN(a) || Double.isNaN(b)) return p + a + b;
		//R_Q_P01_check(p);
		if ((log_p	&& p > 0) || (!log_p && (p < 0 || p > 1)) ) return Double.NaN;
		if (Double.isInfinite(a) || Double.isInfinite(b) || b < a) return Double.NaN;
		if (b == a) return a;

		//return a + R_DT_qIv(p) * (b - a);
		p = (log_p ? (lower_tail ? exp(p) : - expm1(p)) : (lower_tail ? (p) : (0.5 - (p) + 0.5)));
		return a + p * (b - a);
	}

	public static final double random(double a, double b, QRandomEngine random)
	{
		if (Double.isInfinite(a) || Double.isInfinite(b) || b < a) return Double.NaN;
		if (a == b)
			return a;
		else {
			double u;
			/* This is true of all builtin generators, but protect against user-supplied ones */
			do {u = random.nextDouble();} while (u <= 0 || u >= 1);
			return a + (b - a) * u;
		}
	}
}
