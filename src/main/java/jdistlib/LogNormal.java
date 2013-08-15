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

public class LogNormal {
	public static final double density(double x, double meanlog, double sdlog, boolean give_log)
	{
		double y;
		if (Double.isNaN(x) || Double.isNaN(meanlog) || Double.isNaN(sdlog)) return x + meanlog + sdlog;
		if(sdlog <= 0) return Double.NaN;
		if(x <= 0) return (give_log ? Double.NEGATIVE_INFINITY : 0.);
		y = (log(x) - meanlog) / sdlog;
		return (give_log ?
				-(M_LN_SQRT_2PI   + 0.5 * y * y + log(x * sdlog)) :
					M_1_SQRT_2PI * exp(-0.5 * y * y)  /	 (x * sdlog));
		/* M_1_SQRT_2PI = 1 / sqrt(2 * pi) */
	}

	public static final double cumulative(double x, double meanlog, double sdlog, boolean lower_tail, boolean log_p)
	{
		if (Double.isNaN(x) || Double.isNaN(meanlog) || Double.isNaN(sdlog)) return x + meanlog + sdlog;
		if(sdlog <= 0) return Double.NaN;
		if (x > 0) return Normal.cumulative(log(x), meanlog, sdlog, lower_tail, log_p);
		return (log_p ? Double.NEGATIVE_INFINITY : 0.);
	}

	public static final double quantile(double p, double meanlog, double sdlog, boolean lower_tail, boolean log_p)
	{
		if (Double.isNaN(p) || Double.isNaN(meanlog) || Double.isNaN(sdlog)) return p + meanlog + sdlog;
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
		return exp(Normal.quantile(p, meanlog, sdlog, lower_tail, log_p));
	}

	public static final double random(double meanlog, double sdlog, QRandomEngine random)
	{
		if(Double.isNaN(meanlog) || Double.isInfinite(sdlog) || sdlog < 0.) return Double.NaN;
		return exp(Normal.random(meanlog, sdlog, random));
	}
}
