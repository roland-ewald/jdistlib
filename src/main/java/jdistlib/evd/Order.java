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

import jdistlib.generic.GenericDistribution;

import static java.lang.Math.*;
import static jdistlib.MathFunctions.*;

/**
 * Order distribution.
 * Taken from EVD package of R
 *
 */
public class Order {
	public static final double density(double x, GenericDistribution dist, int mlen, int j, boolean largest, boolean log)
	{
		if (mlen <= 0 || j <= 0 || j > mlen)
			return Double.NaN;
		if (!largest)
			j = mlen + 1 - j;
		double dens = dist.density(x, true);
		if (Double.isInfinite(dens))
			return Double.NEGATIVE_INFINITY;
		double cum = dist.cumulative(x, true);
		cum = (mlen - j) * log(cum) + (j - 1) * log (1 - cum);
		x = lgammafn(mlen + 1) - lgammafn(j) - lgammafn(mlen - j + 1) + dens + cum;
		return !log ? exp(x) : x;
	}

	public static final double cumulative(double q, GenericDistribution dist, int mlen, int j, boolean largest, boolean lower_tail)
	{
		if (mlen <= 0 || j <= 0 || j > mlen)
			return Double.NaN;
		return 0;
	}
}
