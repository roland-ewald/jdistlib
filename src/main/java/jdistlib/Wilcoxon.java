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
import static jdistlib.Constants.*;
import static jdistlib.MathFunctions.*;
import jdistlib.rng.QRandomEngine;

/**<pre>
  SYNOPSIS

    #include <Rmath.h>
    double dwilcox(double x, double m, double n, int give_log)
    double pwilcox(double x, double m, double n, int lower_tail, int log_p)
    double qwilcox(double x, double m, double n, int lower_tail, int log_p);
    double rwilcox(double m, double n)

  DESCRIPTION

    dwilcox	The density of the Wilcoxon distribution.
    pwilcox	The distribution function of the Wilcoxon distribution.
    qwilcox	The quantile function of the Wilcoxon distribution.
    rwilcox	Random variates from the Wilcoxon distribution. </pre>
 * 
 * <P>NOTE: Since the computation of Wilcoxon distribution is highly dependent on
 * the matrix <tt>w</tt>, the dimensions of which depends on m and n (parameters of the 
 * Wilcoxon distribution), I decided to make this class as dynamic. -- Roby Joehanes
 *  
 */
public class Wilcoxon {
	protected QRandomEngine random;
	protected double[][][] w;

	public Wilcoxon(int m, int n) {
		init(m,n);
	}

	protected void init(int m, int n) {
	    int i;

	    if (m > n) {
	    	i = n; n = m; m = i;
	    }
	    w = new double[m+1][n+1][];
	}

	public int getM()
	{	return w.length - 1; }

	public int getN()
	{	return w[0].length - 1; }

	protected double count(int k, int m, int n)
	{
		int c, u, i, j, l;

		u = m * n;
		if (k < 0 || k > u)
			return(0);
		c = (int)(u / 2);
		if (k > c)
			k = u - k; /* hence  k <= floor(u / 2) */
		if (m < n) {
			i = m; j = n;
		} else {
			i = n; j = m;
		} /* hence  i <= j */

		if (j == 0) /* and hence i == 0 */
			return (k == 0) ? 1: 0;

		/* We can simplify things if k is small.  Consider the Mann-Whitney 
	           definition, and sort y.  Then if the statistic is k, no more 
	           than k of the y's can be <= any x[i], and since they are sorted 
	           these can only be in the first k.  So the count is the same as
	           if there were just k y's. 
		 */
		if (j > 0 && k < j) return count(k, i, k);    

		if (w[i][j] == null) {
			w[i][j] = new double[c+1];
			for (l = 0; l <= c; l++)
				w[i][j][l] = -1;
		}
		if (w[i][j][k] < 0) {
			if (j == 0) /* and hence i == 0 */
				w[i][j][k] = (k == 0) ? 1 : 0;
			else
				w[i][j][k] = count(k - j, i - 1, j) + count(k, i, j - 1);
		}
		return(w[i][j][k]);
	}

	public double density(int x, boolean give_log)
	{
		int m = w.length - 1, n = w[0].length - 1;
	    double d;

	    /* NaNs propagated correctly */
	    if (Double.isNaN(x) || Double.isNaN(m) || Double.isNaN(n)) return(x + m + n);
	    //m = floor(m + 0.5);
	    //n = floor(n + 0.5);
	    if (m <= 0 || n <= 0) return Double.NaN;

	    //if (abs(x - floor(x + 0.5)) > 1e-7) return (give_log ? Double.NEGATIVE_INFINITY : 0.);
	    // x = floor(x + 0.5);
	    if ((x < 0) || (x > m * n))
		return (give_log ? Double.NEGATIVE_INFINITY : 0.);

	    // w_init_maybe(m, n);
	    d = give_log ?
		log(count((int) x, (int) m, (int) n)) - lchoose(m + n, n) :
			count((int) x, (int) m, (int) n)  /	 choose(m + n, n);

	    return(d);
	}

	public double cumulative(int q, boolean lower_tail, boolean log_p)
	{
		int m = w.length - 1, n = w[0].length - 1;
		int i;
		double c, p;

		if (Double.isNaN(q) || Double.isNaN(m) || Double.isNaN(n)) return(q + m + n);
		if (Double.isInfinite(m) || Double.isInfinite(n)) return Double.NaN;
		//m = floor(m + 0.5);
		//n = floor(n + 0.5);
		if (m <= 0 || n <= 0) return Double.NaN;

		// q = floor(q + 1e-7);

		if (q < 0.0) return(lower_tail ? (log_p ? Double.NEGATIVE_INFINITY : 0.) : (log_p ? 0. : 1.));
		if (q >= m * n) return(lower_tail ? (log_p ? 0. : 1.) : (log_p ? Double.NEGATIVE_INFINITY : 0.));

		//w_init_maybe(m, n);
		c = choose(m + n, n);
		p = 0;
		/* Use summation of probs over the shorter range */
		if (q <= (m * n / 2)) {
			for (i = 0; i <= q; i++)
				p += count(i, m, n) / c;
		}
		else {
			q = m * n - q;
			for (i = 0; i < q; i++)
				p += count(i, m, n) / c;
			lower_tail = !lower_tail; /* p = 1 - p; */
		}
		//return(R_DT_val(p));
		return (lower_tail ? (log_p ? log(p) : (p)) : (log_p ? log1p(-(p)) : (0.5 - (p) + 0.5)));
	}

	public double quantile(double x, boolean lower_tail, boolean log_p)
	{
		int m = w.length - 1, n = w[0].length - 1, q;
		double c, p;

		if (Double.isNaN(x) || Double.isNaN(m) || Double.isNaN(n)) return(x + m + n);
		if(Double.isInfinite(x) || Double.isInfinite(m) || Double.isInfinite(n)) return Double.NaN;
		//R_Q_P01_check(x);
		if ((log_p && x > 0) || (!log_p && (x < 0 || x > 1)) ) return Double.NaN;

		//m = floor(m + 0.5);
		//n = floor(n + 0.5);
		if (m <= 0 || n <= 0) return Double.NaN;

		if (x == (lower_tail ? (log_p ? Double.NEGATIVE_INFINITY : 0.) : (log_p ? 0. : 1.))) return(0);
		if (x == (lower_tail ? (log_p ? 0. : 1.) : (log_p ? Double.NEGATIVE_INFINITY : 0.)))
			return(m * n);

		if(log_p || !lower_tail)
			//x = R_DT_qIv(x); /* lower_tail,non-log "p" */
			x = (log_p ? (lower_tail ? exp(x) : - expm1(x)) : (lower_tail ? (x) : (0.5 - (x) + 0.5)));

		//w_init_maybe(m, n);
		c = choose(m + n, n);
		p = 0;
		q = 0;
		if (x <= 0.5) {
			x = x - 10 * DBL_EPSILON;
			for (;;) {
				p += count(q, m, n) / c;
				if (p >= x)
					break;
				q++;
			}
		}
		else {
			x = 1 - x + 10 * DBL_EPSILON;
			for (;;) {
				p += count(q, m, n) / c;
				if (p > x) {
					q = m * n - q;
					break;
				}
				q++;
			}
		}
		return(q);
	}

	public void setRandomEngine(QRandomEngine rand)
	{	random = rand; }

	public QRandomEngine getRandomEngine()
	{	return random; }

	public double random()
	{	return random(random); }

	public double random(QRandomEngine random)
	{
		int m = w.length - 1, n = w[0].length - 1;
		int i, j, k, x[];
		double r;

		/* NaNs propagated correctly */
		if (Double.isNaN(m) || Double.isNaN(n)) return(m + n);
		//m = floor(m + 0.5);
		//n = floor(n + 0.5);
		if ((m < 0) || (n < 0))
			return Double.NaN;

		if ((m == 0) || (n == 0))
			return(0);

		r = 0.0;
		k = (int) (m + n);
		x = new int[k]; // (int *) calloc((size_t) k, sizeof(int));
		for (i = 0; i < k; i++)
			x[i] = i;
		for (i = 0; i < n; i++) {
			j = (int) floor(k * random.nextDouble());
			r += x[j];
			x[j] = x[--k];
		}
		x = null; //free(x);
		return(r - n * (n - 1) / 2);
	}
}
