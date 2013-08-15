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

import static java.lang.Math.abs;
import static java.lang.Math.floor;
import static java.util.Arrays.sort;
import static jdistlib.disttest.Utils.calculate_ecdf;


/**
 * Comparing two distributions
 * @author Roby Joehanes
 *
 */
public class DistributionTest {
	/**
	 * Compute the Kolmogorov-Smirnov test to test between two distribution.
	 * Note: I don't multiply the D score with sqrt(nX*nY / (nX + nY)), which
	 * is needed for P-value computation
	 * 
	 * @param X an array with length of nX
	 * @param Y an array with length of nY
	 * @return
	 */
	public static final double kolmogorov_smirnov_statistic(double[] X, double[] Y)
	{
		int
			nX = X.length,
			nY = Y.length,
			idxX = 0,
			idxY = 0;
		double
			sortedX[] = X.clone(),
			sortedY[] = Y.clone(),
			maxDiv = 0;
		sort(sortedX);
		sort(sortedY);
	
		// Pathological case
		if (sortedX[nX - 1] < sortedY[0] || sortedY[nY - 1] < sortedX[0])
			return 1.0;
	
		// Scan for duplicate values
		double
			cdfX[] = calculate_ecdf(sortedX),
			cdfY[] = calculate_ecdf(sortedY),
			pX = 0,
			pY = 0,
			div = 0;
	
		while (idxX < nX && idxY < nY)
		{
			double
				x = sortedX[idxX],
				y = sortedY[idxY];
			if (y < x)
			{
				pY = cdfY[idxY];
				idxY++;
			}
			else if (y > x)
			{
				pX = cdfX[idxX];
				idxX++;
			}
			else
			{
				pX = cdfX[idxX];
				pY = cdfY[idxY];
				idxX++; idxY++;
			}
			div = abs(pX - pY);
			//div = abs(idxX / ((double) nX) - (idxY - 1.0) / nY);
			if (div > maxDiv)
				maxDiv = div;
		}
	
		return maxDiv;
	}

	/**
	 * Compute the P-value out of the D-score produced by <tt>kolmogorov_smirnov_statistic</tt>.
	 * 
	 * @param maxDiv
	 * @param lengthX
	 * @param lengthY
	 * @return
	 */
	public static final double kolmogorov_smirnov_pvalue(double maxDiv, int lengthX, int lengthY)
	{
		/*
		Set<Double> set = new HashSet<Double>();
		for (double x: X)
			set.add(x);
		m = set.size();
		set.clear();
		for (double y: Y)
			set.add(y);
		n = set.size();
		set.clear();
		set = null;
		//*/
	
		if (lengthX > lengthY)
		{
			int temp = lengthY;
			lengthY = lengthX;
			lengthX = temp;
		}
	
		double
			q = floor(maxDiv * lengthX * lengthY - 1e-7) / (lengthX * lengthY),
			u[] = new double[lengthY + 1],
			md = lengthX,
			nd = lengthY;
	
		for (int j = 0; j <= lengthY; j++)
			u[j] = (j / nd) > q ? 0: 1;
		for(int i = 1; i <= lengthX; i++)
		{
			double w = (double)(i) / ((double)(i + lengthY));
			u[0] = (i / md) > q ? 0 : w * u[0];
			for(int j = 1; j <= lengthY; j++)
				u[j] = abs(i / md - j / nd) > q ? 0 : w * u[j] + u[j - 1];
		}
		return 1 - u[lengthY];
	}

}
