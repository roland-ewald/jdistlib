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
package jdistlib.disttest;

/**
 * 
 * @author Roby Joehanes
 *
 */
public class Utils {
	/**
	 * Calculate the CDF of empirical distribution
	 * @param sorted
	 * @return
	 */
	public static final double[] calculate_ecdf(double[] sorted)
	{
		int n = sorted.length;
		double
			cdf[] = new double[n],
			last = sorted[n - 1];
	
		for (int i = 0; i < n; i++)
			cdf[i] = (i + 1.0) / n;
	
		// Detect duplicate
		for (int i = n - 2; i >= 0; i--)
		{
			double val = sorted[i];
			if (val == last)
				cdf[i] = cdf[i + 1];
			last = val;
		}
		return cdf;
	}
}
