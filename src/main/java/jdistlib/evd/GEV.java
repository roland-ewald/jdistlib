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

import jdistlib.Exponential;
import jdistlib.rng.QRandomEngine;

import static java.lang.Math.*;

/**
 * Generalized extreme value distribution.
 * Taken from EVD package of R
 */
public class GEV {
	public static final double density(double x, double loc, double scale, double shape, boolean log)
	{
		if (scale <= 0)
			return Double.NaN;
		x = (x - loc) / scale;
		if (shape == 0)
			x = - log(scale) - x - exp(-x);
		else {
			x = 1 + shape * x;
			x = x <= 0 ? Double.NEGATIVE_INFINITY : -log(scale) - pow(x, -1.0/shape) - (1.0/shape + 1) * log(x);
		}
		return !log ? exp(x) : x;
	}

	public static final double cumulative(double q, double loc, double scale, double shape, boolean lower_tail)
	{
		if (scale <= 0)
			return Double.NaN;
		q = (q - loc) / scale;
		q = shape == 0 ? exp(-exp(-q)) : exp(-pow(max(1 + shape * q, 0), -1.0/shape));
		return !lower_tail ? 1 - q : q;
	}

	public static final double quantile(double p, double loc, double scale, double shape, boolean lower_tail)
	{
		if (p <= 0 || p >= 1 || scale < 0)
			return Double.NaN;
		if (!lower_tail)
			p = 1 - p;
		return shape == 0 ? loc - scale * log(-log(p)) : loc + scale * ((pow(-log(p), -shape) - 1) / shape);
	}

	public static final double random(double loc, double scale, double shape, QRandomEngine random)
	{
		if (scale < 0)
			return Double.NaN;
		return shape == 0 ?
			(loc - scale * log(Exponential.random_standard(random))) :
			(loc + scale * ((pow(Exponential.random_standard(random), -shape) - 1) / shape) );
	}
}
