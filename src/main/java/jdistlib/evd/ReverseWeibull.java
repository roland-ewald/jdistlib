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
 * Reverse Weibull distribution.
 * Taken from EVD package of R
 */
public class ReverseWeibull {
	public static final double density(double x, double loc, double scale, double shape, boolean log)
	{
		if (scale <= 0 || shape <= 0)
			return Double.NaN;
		x = (x - loc) / scale;
		if (x >= 0)
			return Double.NEGATIVE_INFINITY;
		x = log(shape / scale) + (shape - 1.0) * log(-x) - pow(-x, shape);
		return !log ? exp(x) : x;
	}

	public static final double cumulative(double q, double loc, double scale, double shape, boolean lower_tail)
	{
		if (scale <= 0 || shape <= 0)
			return Double.NaN;
		q = exp(-pow(-min((q - loc) / scale, 0), shape));
		return !lower_tail ? 1 - q : q;
	}

	public static final double quantile(double p, double loc, double scale, double shape, boolean lower_tail)
	{
		if (p <= 0 || p >= 1 || scale < 0 || shape <= 0)
			return Double.NaN;
		if (!lower_tail)
			p = 1 - p;
		return loc - scale * pow(-log(p), 1.0 / shape);
	}

	public static final double random(double loc, double scale, double shape, QRandomEngine random)
	{
		if (scale < 0 || shape <= 0)
			return Double.NaN;
		return loc - scale * pow(Exponential.random_standard(random), 1.0 / shape);
	}
}
