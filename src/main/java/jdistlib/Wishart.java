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
package jdistlib;

import jdistlib.rng.QRandomEngine;

import static java.lang.Math.sqrt;
import static jdistlib.matrix.QMatrixUtils.*;

/**
 * Wishart distribution
 * By: Roby Joehanes
 */
public class Wishart {

	/**
	 * <P>Sample Wishart distribution with n degree of freedom and LL' scale parameter.
	 * Wishart distribution is multivariate generalization of Gamma distribution.
	 * Remember that Gamma distribution is a generalization of both Chi-square and
	 * exponential distribution.
	 * 
	 * <P>I use L instead of the usual V to denote the scale parameter. Note that V = LL'.
	 * Compute L using Cholesky decomposition. L is a lower triangular matrix.
	 * 
	 * <P>The algorithm described here is due to:<br>
	 * Smith, W. B. and Hocking, R. R.; Algorithm AS 53: Wishart Variate Generator, Applied Statistician, v.21(341--345), 1972
	 * http://www.jstor.org/stable/2346290
	 * 
	 * <P>See also:
	 * Gentle, J. E., Random Number Generation and Monte Carlo Methods, 2nd ed., Springer-Verlag, 2003, pp. 199
	 * 
	 * <P>Note: If you want to generate an inverse-wishart distribution, simply generate
	 * a matrix from Wishart distribution and invert it.
	 * 
	 * @param n degree of freedom
	 * @param L square matrix of d x d, lower triangular Cholesky factor of sigma matrix
	 * @param rand Randomizer
	 * @return
	 */
	public static final double[][] random(double n, double[][] L, QRandomEngine rand)
	{
		int d = L.length;
		assert (d == L[0].length); // Crude check for square-ness

		double
			z[][] = new double[d][d],
			y[] = new double[d],
			B[][] = new double[d][d];

		for (int i = 1; i < d; i++)
			for (int j = 0; j < i; j++)
				z[i][j] = rand.nextGaussian();

		for (int i = 0; i < d; i++)
		{
			double sum = y[i] = ChiSquare.random(n - i + 1, rand);
			for (int j = 0; j < i; j++)
			{
				double val = z[j][i];
				sum += val * val;
			}
			B[i][i] = sum / n;
			if (i > 0)
				B[0][i] = B[i][0] = z[i][0] * sqrt(y[0]) / n;
		}

		for (int i = 2; i < d; i++)
			for (int j = 1; j < i; j++)
			{
				double sum = z[i][j] * sqrt(y[j]);
				for (int k = 0; k < j; k++)
					sum += z[i][k] * z[j][k];
				B[i][j] = B[j][i] = sum / n;
			}

		return calculateXYt(calculateXY(L,B),L);
	}
}
