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

import static java.util.Arrays.*;
import static jdistlib.MathFunctions.*;
import jdistlib.rng.QRandomEngine;

/**
 * Ansari-Bradley test statistic
 * 
 *
 */
public class Ansari {
	static final double cansari(int k, int m, int n, double w[][][])
	{
		int l, u;

		l = (m + 1) * (m + 1) / 4;
		u = l + m * n / 2;

		if ((k < l) || (k > u))
			return(0);

		if (w[m][n] == null) {
			w[m][n] = new double[u+1];// (double *) R_alloc(u + 1, sizeof(double));
			//memset(w[m][n], '\0', (u + 1) * sizeof(double));
			//for (int i = 0; i <= u; i++) w[m][n][i] = -1;
			fill(w[m][n], -1);
		}
		if (w[m][n][k] < 0) {
			if (m == 0)
				w[m][n][k] = (k == 0) ? 1 : 0;
			else if (n == 0)
				w[m][n][k] = (k == l) ? 1 : 0;
			else
				w[m][n][k] = cansari(k, m, n - 1, w) + cansari(k - (m + n + 1) / 2, m - 1, n, w);
		}
		return(w[m][n][k]);
	}

	public static final double[] density(int[] x, int m, int n)
	{
		int i, len = x.length;
		double[][][] w = new double[m+1][n+1][];
		double[] result = new double[len];
		for (i = 0; i < len; i++)
			result[i] = cansari((int) x[i], m, n, w) / choose(m + n, m);
		return result;
	}

	public static final double density(int x, int m, int n)
	{
		return cansari(x, m, n, new double[m+1][n+1][]) / choose(m + n, m);
	}

	public static final double[] cumulative(int[] x, int m, int n)
	{
		int i, j, l, u, len = x.length;
		double c, p;
		double[][][] w = new double[m+1][n+1][];
		double[] result = new double[len];

		l = (m + 1) * (m + 1) / 4;
		u = l + m * n / 2;
		c = choose(m + n, m);
		for (i = 0; i < len; i++) {
			//q = floor(x[i] + 1e-7);
			if (x[i] < l)
				result[i] = 0;
			else if (x[i] > u)
				result[i] = 1;
			else {
				int qq = x[i];
				p = 0;
				for (j = l; j <= qq; j++) {
					p += cansari(j, m, n, w);
				}
				result[i] = p / c;
			}
		}
		return result;
	}

	public static final double cumulative(int x, int m, int n)
	{
		int j, l, u;
		double c, p;
		double[][][] w = new double[m+1][n+1][];

		l = (m + 1) * (m + 1) / 4;
		u = l + m * n / 2;
		c = choose(m + n, m);
		if (x < l)
			return 0;
		else if (x > u)
			return 1;
		p = 0;
		for (j = l; j <= x; j++) {
			p += cansari(j, m, n, w);
		}
		return p / c;
	}

	public static final int[] quantile(double[] x, int m, int n)
	{
		int i, l, u, q, len = x.length;
		double c, p, xi;
		double[][][] w = new double[m+1][n+1][];
		int[] result = new int[len];

		l = (m + 1) * (m + 1) / 4;
		u = l + m * n / 2;
		c = choose(m + n, m);
		for (i = 0; i < len; i++) {
			xi = x[i];
			if(xi < 0 || xi > 1) {
				//error(_("probabilities outside [0,1] in qansari()"));
				System.err.println("probabilities outside [0,1] in Ansari.quantile");
				result[i] = Integer.MIN_VALUE;
			}
			if(xi == 0)
				result[i] = l;
			else if(xi == 1)
				result[i] = u;
			else {
				p = 0;
				q = 0;
				for(;;) {
					p += cansari(q, m, n, w) / c;
					if (p >= xi)
						break;
					q++;
				}
				result[i] = q;
			}
		}
		return result;
	}

	public static final int quantile(double xi, int m, int n)
	{
		int l, u, q;
		double c, p;
		double[][][] w = new double[m+1][n+1][];

		l = (m + 1) * (m + 1) / 4;
		u = l + m * n / 2;
		c = choose(m + n, m);
		if(xi < 0 || xi > 1) {
			//error(_("probabilities outside [0,1] in qansari()"));
			System.err.println("probabilities outside [0,1] in Ansari.quantile");
			return Integer.MIN_VALUE;
		}
		if(xi == 0)
			return l;
		else if(xi == 1)
			return u;
		p = 0;
		q = 0;
		for(;;) {
			p += cansari(q, m, n, w) / c;
			if (p >= xi)
				break;
			q++;
		}
		return q;
	}

	/**
	 * Ansari RNG by inversion -- WARNING: Untested
	 * @param m
	 * @param n
	 * @param random
	 * @return
	 */
	public static final double random(int m, int n, QRandomEngine random)
	{
		double u1 = random.nextDouble();
		u1 = (int) (134217728 * u1) + random.nextDouble();
		u1 = quantile(u1 / 134217728, m, n);
		return u1;
	}
}
