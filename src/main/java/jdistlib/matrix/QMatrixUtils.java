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
package jdistlib.matrix;

/**
 * Matrix utilities
 * @author Roby Joehanes
 *
 */
public class QMatrixUtils {
	/**
	 * Calculate XY. Eliminate the need to instantiate a matrix.
	 * @param X (n x m)
	 * @param Y (m x p)
	 * @return
	 */
	public static final double[][] calculateXY(double[][] X, double[][] Y)
	{
		int
			rows = X.length,
			m = Y.length,
			cols = Y[0].length;
		assert (m == X.length);
		double[][] result = new double[rows][cols];
	
		for (int i = 0; i < rows; i++)
		{
			double[] Xrow = X[i];
			for (int j = 0; j < cols; j++)
			{
				double sum = 0;
				for (int k = 0; k < m; k++)
					sum += Xrow[k] * Y[k][j];
				result[i][j] = sum;
			}
		}
		return result;
	}

	/**
	 * Calculate XY'. Eliminate the need to transpose.
	 * @param X
	 * @param Y
	 * @return
	 */
	public static final double[][] calculateXYt(double[][] X, double[][] Y)
	{
		int
			rows = X.length,
			m = Y[0].length,
			cols = Y.length;
		assert (m == X.length);
		double[][] result = new double[rows][cols];
	
		for (int i = 0; i < rows; i++)
		{
			double[] Xrow = X[i];
			for (int j = 0; j < cols; j++)
			{
				double
					sum = 0,
					Yrow[] = Y[j];
				for (int k = 0; k < m; k++)
					sum += Xrow[k] * Yrow[k];
				result[i][j] = sum;
			}
		}
		return result;
	}

	/**
	 * Calculate X'Y. Eliminate the need to transpose.
	 * @param X
	 * @param Y
	 * @return
	 */
	public static final double[][] calculateXtY(double[][] X, double[][] Y)
	{
		int
			rows = X[0].length,
			m = Y.length,
			cols = Y[0].length;
		assert (m == X.length);
		double[][] result = new double[rows][cols];
	
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
			{
				double sum = 0;
				for (int k = 0; k < m; k++)
					sum += X[k][i] * Y[k][j];
				result[i][j] = sum;
			}
		return result;
	}
}
