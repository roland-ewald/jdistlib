/*
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

/**
 * @author Roby Joehanes
 */
public class Logarithmic {
	public static final double density(double x, double mu, boolean give_log)
	{
		if (Double.isNaN(x) || Double.isNaN(mu))
			return x + mu;
		if (mu <= 0 || mu >= 1)
			return Double.NaN;
		double logfy = x * log(mu) - log(x) -log(-log(1-mu));
		return give_log ? logfy : exp(logfy);
	}

	public static final double cumulative(double q, double mu, boolean lower_tail, boolean log_p)
	{
		if (Double.isNaN(q) || Double.isNaN(mu))
			return q + mu;
		if (mu <= 0 || mu >= 1 || q <= 0)
			return Double.NaN;
		double sum = 0;
		for (int i = 0; i < q; i++)
			sum += density(i, mu, false);
		sum = lower_tail ? sum : 1 - sum;
		return log_p ? log(sum) : sum;
	}

	public static final double quantile(double p, double mu, boolean lower_tail, boolean log_p)
	{	return quantile(p, mu, lower_tail, log_p, 10000); }

	public static final double quantile(double p, double mu, boolean lower_tail, boolean log_p, int max_value)
	{
		if (Double.isNaN(p) || Double.isNaN(mu))
			return p + mu;
		if (mu <= 0 || mu >= 1)
			return Double.NaN;
		if (p < 0)
			return Double.NEGATIVE_INFINITY;
		if (p > 1)
			return Double.POSITIVE_INFINITY;
		p = log_p ? exp(p) : p;
		p = lower_tail ? p : 1 - p;
		for (int i = 1; i <= max_value; i++)
			if (p <= cumulative(i, mu, true, false))
				return i;
		return Double.POSITIVE_INFINITY;
	}

	public static final double random(double mu, QRandomEngine random)
	{
		return quantile(random.nextDouble(), mu, false, false);
	}
}
