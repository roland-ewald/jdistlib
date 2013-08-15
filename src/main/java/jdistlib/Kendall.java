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

import static java.lang.Math.*;
import static jdistlib.MathFunctions.*;
import jdistlib.rng.QRandomEngine;

public class Kendall {
	static final double count(int k, int n, double w[][]) {
		int i, u;
		double s;

		u = (n * (n - 1) / 2);
		if ((k < 0) || (k > u)) return(0);
		if (w[n] == null) {
			w[n] = new double[u+1]; // (double *) R_alloc(u + 1, sizeof(double));
			//memset(w[n], '\0', sizeof(double) * (u+1));
			for (i = 0; i <= u; i++)
				w[n][i] = -1;
		}
		if (w[n][k] < 0) {
			if (n == 1)
				w[n][k] = (k == 0) ? 1 : 0;
			else {
				s = 0;
				for (i = 0; i < n; i++)
					s += count(k - i, n - 1, w);
				w[n][k] = s;
			}
		}
		return(w[n][k]);
	}

	public static final double density(double x, int n){
		double w[][] = new double[n+1][];
		if (abs(x - floor(x + 0.5)) > 1e-7) return 0;
		return count((int) x, (int) (n / gammafn(n + 1)), w);
	}

	/**
	 * Cumulative density function of Kendall distribution.
	 * <P>Kendall statistics: x = round((rho + 1) * n * (n-1) / 4);
	 * <P>Two-sided test: min(1, 2*((q > n*(n-1)/4) ? 1-cumulative(x-1,n) : cumulative(x,n)));
	 * <P>Greater test: 1-cumulative(x-1,n)
	 * <P>Less test: cumulative(x,n)
	 * @param x
	 * @param n
	 * @return
	 */
	public static final double cumulative(double x, int n) {
		double p, q;
		double w[][] = new double[n+1][];

		//w = (double **) R_alloc(*n + 1, sizeof(double *));
		//memset(w, '\0', sizeof(double*) * (*n+1));

		q = floor(x + 1e-7);
		if (q < 0)
			return 0;
		if (q > (n * (n - 1) / 2))
			return 1;
		p = 0;
		for (int j = 0; j <= q; j++) {
			p += count(j, n, w);
		}
		return p / gammafn(n + 1);
	}

	static final double do_search(double y, double[] z, double p, int x, double incr)
	{
		if(z[0] >= p) {
			/* search to the left */
			for(;;) {
				double newz = cumulative(y - incr, x);
				if(y == 0 || newz < p)
					return y;
				y = max(0, y - incr);
				z[0] = newz;
			}
		}
		else {		/* search to the right */
			for(;;) {
				y = y + incr;
				if((z[0] = cumulative(y, x)) >= p)
					return y;
			}
		}
	}

	/**
	 * Quantile search by Cornish-Fisher expansion. WARNING: Untested!
	 * @param p
	 * @param x
	 * @return
	 */
	public static final double quantile(double p, int x)
	{
		if (Double.isNaN(p) || Double.isInfinite(p)) return p;
		if (p < 0 || p > 1) return Double.NaN;
		double z[] = new double[1], y, mu, sigma;

		/* y := approx.value (Cornish-Fisher expansion) :  */
		z[0] = Normal.quantile(p, 0., 1., true, false);
		mu = x * (x-1)/4.0;
		sigma = sqrt(x * (x+1.0) * (2*x + 5.5) / 72);
		y = sigma * z[0] + mu - 0.5;

		z[0] = cumulative(y, x);
		double incr = floor(y * 0.001), oldincr;
		do {
			oldincr = incr;
			y = do_search(y, z, p, x, incr);
			incr = max(1, floor(incr/100));
		} while(oldincr > 1 && incr > x*1e-15);
		return y;
	}

	/**
	 * Kendall RNG by inversion -- WARNING: Untested
	 * @param x
	 * @param random
	 * @return
	 */
	public static final double random(int x, QRandomEngine random)
	{
		double u1 = random.nextDouble();
		u1 = (int) (134217728 * u1) + random.nextDouble();
		u1 = quantile(u1 / 134217728, x);
		return u1;
	}
}
