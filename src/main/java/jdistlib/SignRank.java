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
import jdistlib.rng.QRandomEngine;

public class SignRank {
	protected double[] w;
	protected int n;
	protected QRandomEngine random;

	public SignRank(int n) {
		
	}

	protected void init(int _n)
	{
		n = _n;
	    int c = (n * (n + 1) / 4);
		w = new double[c+1];
	}

	protected double count(int k, int n)
	{
		int c, u, j;

		u = n * (n + 1) / 2;
		c = (u / 2);

		if (k < 0 || k > u)
			return 0;
		if (k > c)
			k = u - k;

		if (n == 1)
			return 1.;
		if (w[0] == 1.)
			return w[k];

		w[0] = w[1] = 1.;
		for(j=2; j < n+1; ++j) {
			int i, end = min(j*(j+1)/2, c);
			for(i=end; i >= j; --i) {
				w[i] += w[i-j];
			}
		}

		return w[k];
	}

	public double density(int x, boolean give_log)
	{
	    double d;

	    /* NaNs propagated correctly */
	    if (Double.isNaN(x)) return(x + n);
	    //n = floor(n + 0.5);
	    //if (n <= 0) return Double.NaN;

	    // if (abs(x - floor(x + 0.5)) > 1e-7) return(R_D__0);
	    // x = floor(x + 0.5);
	    if ((x < 0) || (x > (n * (n + 1) / 2)))
	    	return(give_log? Double.NEGATIVE_INFINITY : 0.);

	    //d = R_D_exp(log(csignrank(x, n)) - n * M_LN2);
	    d = log(count(x, n)) - n * M_LN2;

	    return (give_log ? (d) : exp(d));
	}

	public double cumulative(int x, boolean lower_tail, boolean log_p)
	{
		int n = w.length - 1;
		int i;
		double f, p;

		if (Double.isNaN(x)) return(x + n);
		//if (Double.isInfinite(n)) return Double.NaN;
		//n = floor(n + 0.5);
		//if (n <= 0) return Double.NaN;

		//x = floor(x + 1e-7);
		if (x < 0.0) return(lower_tail ? (log_p ? Double.NEGATIVE_INFINITY : 0.) : (log_p ? 0. : 1.));
		if (x >= n * (n + 1) / 2)
			return(lower_tail ? (log_p ? 0. : 1.) : (log_p ? Double.NEGATIVE_INFINITY : 0.));

		f = exp(- n * M_LN2);
		p = 0;
		if (x <= (n * (n + 1) / 4)) {
			for (i = 0; i <= x; i++)
				p += count(i, n) * f;
		}
		else {
			x = n * (n + 1) / 2 - x;
			for (i = 0; i < x; i++)
				p += count(i, n) * f;
			lower_tail = !lower_tail; /* p = 1 - p; */
		}

		//return(R_DT_val(p));
		return (lower_tail ? (log_p ? log(p) : (p))  : (log_p ? log1p(-(p)) : (0.5 - (p) + 0.5)));
	}

	public double quantile(double x, boolean lower_tail, boolean log_p)
	{
		int n = w.length - 1;
		double f, p;//, q;
		int q;

		if (Double.isNaN(x)) return(x + n);
		if (Double.isInfinite(x)) return Double.NaN;
		//R_Q_P01_check(x);
		if ((log_p	&& x > 0) || (!log_p && (x < 0 || x > 1)) ) return Double.NaN;

		//n = floor(n + 0.5);
		//if (n <= 0) return Double.NaN;

		if (x == (lower_tail ? (log_p ? Double.NEGATIVE_INFINITY : 0.) : (log_p ? 0. : 1.)))
			return(0);
		if (x == (lower_tail ? (log_p ? 0. : 1.) : (log_p ? Double.NEGATIVE_INFINITY : 0.)))
			return(n * (n + 1) / 2);

		if(log_p || !lower_tail)
			//x = R_DT_qIv(x); /* lower_tail,non-log "p" */
			x = (log_p ? (lower_tail ? exp(x) : - expm1(x)) : (lower_tail ? (x) : (0.5 - (x) + 0.5)));

		//w_init_maybe(n);
		f = exp(- n * M_LN2);
		p = 0;
		q = 0;
		if (x <= 0.5) {
			x = x - 10 * DBL_EPSILON;
			for (;;) {
				p += count(q, n) * f;
				if (p >= x)
					break;
				q++;
			}
		}
		else {
			x = 1 - x + 10 * DBL_EPSILON;
			for (;;) {
				p += count(q, n) * f;
				if (p > x) {
					q = n * (n + 1) / 2 - q;
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
		int i, k;
		double r;

		//n = floor(n + 0.5);
		//if (n < 0) ML_ERR_return_NAN;

		if (n == 0)
			return(0);
		r = 0.0;
		k = (int) n;
		for (i = 0; i < k; ) {
			r += (++i) * floor(random.nextDouble() + 0.5);
		}
		return(r);
	}
}
