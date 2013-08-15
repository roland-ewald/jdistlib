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

import jdistlib.rng.QRandomEngine;

public class ChiSquare {
	public static final double density(double x, double df, boolean give_log)
	{
	    return Gamma.density(x, df / 2., 2., give_log);
	}

	public static final double cumulative(double x, double df, boolean lower_tail, boolean log_p)
	{
	    return Gamma.cumulative(x, df/2., 2., lower_tail, log_p);
	}

	public static final double quantile(double p, double df, boolean lower_tail, boolean log_p)
	{
	    return Gamma.quantile(p, 0.5 * df, 2.0, lower_tail, log_p);
	}

	public static final double random(double df, QRandomEngine random)
	{
	    if (Double.isInfinite(df) || df < 0.0) return Double.NaN;
	    return Gamma.random(df / 2.0, 2.0, random);
	}
}
