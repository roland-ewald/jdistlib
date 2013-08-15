/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; version 3 of the License.
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
package jdistlib.evd;

import static java.lang.Math.*;
import jdistlib.Exponential;
import jdistlib.rng.QRandomEngine;

/**
 * Generalized Pareto Distribution
 * Taken from EVD package of R
 *
 */
public class GeneralizedPareto {
	public static final double density(double x, double loc, double scale, double shape, boolean log)
	{
		if (scale <= 0)
			return Double.NaN;
		x = (x - loc) / scale;
		if ((1 + shape * x) <= 0)
			return Double.NEGATIVE_INFINITY;
		x = shape == 0 ? -log(scale) - x : -log(scale) - (1.0 / shape + 1) * log(1 + shape * x);
		return !log ? exp(x) : x;
	}

	public static final double cumulative(double q, double loc, double scale, double shape, boolean lower_tail)
	{
		if (scale <= 0)
			return Double.NaN;
		q = max(q - loc, 0) / scale;
		q = shape == 0 ? 1 - exp(-q) : 1 - pow(max(1 + shape * q, 0), -1.0 / shape);
		return !lower_tail ? 1 - q : q;
	}

	public static final double quantile(double p, double loc, double scale, double shape, boolean lower_tail)
	{
		if (p <= 0 || p >= 1 || scale < 0)
			return Double.NaN;
		if (!lower_tail)
			p = 1 - p;
		p = shape == 0 ? 1 - exp(-p) : 1 - pow(max(1 + shape * p, 0), -1.0 / shape);
		return shape == 0 ? (loc - scale * log(p)) :
			(loc + scale * (pow(p, -shape) - 1) / shape);
	}

	public static final double random(double loc, double scale, double shape, QRandomEngine random)
	{
		if (scale < 0)
			return Double.NaN;
		return shape == 0 ?
			(loc + scale * log(Exponential.random_standard(random))) :
			(loc + scale * ((pow(random.nextDouble(), -shape) - 1) / shape) );
	}
}
