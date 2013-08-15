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

public class F {
	public static final double density(double x, double m, double n, boolean give_log)
	{
		double p, q, f, dens;

		if (Double.isNaN(x) || Double.isNaN(m) || Double.isNaN(n)) return x + m + n;
		if (m <= 0 || n <= 0) return Double.NaN;
		if (x < 0.)  return((give_log ? Double.NEGATIVE_INFINITY : 0.));
		if (x == 0.) return(m > 2 ? (give_log ? Double.NEGATIVE_INFINITY : 0.) : (m == 2 ? (give_log ? 0. : 1.) : Double.POSITIVE_INFINITY));
		if (Double.isInfinite(m) && Double.isInfinite(n)) { /* both +Inf */
			if(x == 1.) return Double.POSITIVE_INFINITY;
			/* else */  return (give_log ? Double.NEGATIVE_INFINITY : 0.);
		}
		if (Double.isInfinite(n)) /* must be +Inf by now */
			return(Gamma.density(x, m/2, 2./m, give_log));
		if (m > 1e14) {/* includes +Inf: code below is inaccurate there */
			dens = Gamma.density(1./x, n/2, 2./n, give_log);
			return give_log ? dens - 2*log(x): dens/(x*x);
		}

		f = 1./(n+x*m);
		q = n*f;
		p = x*m*f;

		if (m >= 2) {
			f = m*q/2;
			dens = Binomial.density_raw((m-2)/2, (m+n-2)/2, p, q, give_log);
		}
		else {
			f = m*m*q / (2*p*(m+n));
			dens = Binomial.density_raw(m/2, (m+n)/2, p, q, give_log);
		}
		return(give_log ? log(f)+dens : f*dens);
	}

	public static final double cumulative(double x, double df1, double df2, boolean lower_tail, boolean log_p)
	{
		if (Double.isNaN(x) || Double.isNaN(df1) || Double.isNaN(df2)) return x + df2 + df1;
		if (df1 <= 0. || df2 <= 0.) return Double.NaN;

		//R_P_bounds_01(x, 0., ML_POSINF);
	    if(x <= 0) return (lower_tail ? (log_p ? Double.NEGATIVE_INFINITY : 0.) : (log_p ? 0. : 1.));
	    if(x >= Double.POSITIVE_INFINITY) return (lower_tail ? (log_p ? 0. : 1.) : (log_p ? Double.NEGATIVE_INFINITY : 0.));

		/* move to pchisq for very large values - was 'df1 > 4e5' in 2.0.x,
		       now only needed for df1 = Inf or df2 = Inf {since pbeta(0,*)=0} : */
		if (df2 == Double.POSITIVE_INFINITY) {
			if (df1 == Double.POSITIVE_INFINITY) {
				if(x <  1.) return (lower_tail ? (log_p ? Double.NEGATIVE_INFINITY : 0.) : (log_p ? 0. : 1.));
				if(x == 1.) return (log_p ? -M_LN2 : 0.5);
				if(x >  1.) return (lower_tail ? (log_p ? 0. : 1.) : (log_p ? Double.NEGATIVE_INFINITY : 0.));
			}
			return ChiSquare.cumulative(x * df1, df1, lower_tail, log_p);
		}

		if (df1 == Double.POSITIVE_INFINITY)/* was "fudge"	'df1 > 4e5' in 2.0.x */
			return ChiSquare.cumulative(df2 / x , df2, !lower_tail, log_p);

		/* Avoid squeezing pbeta's first parameter against 1 :  */
		if (df1 * x > df2)
			x = Beta.cumulative(df2 / (df2 + df1 * x), df2 / 2., df1 / 2., !lower_tail, log_p);
		else
			x = Beta.cumulative(df1 * x / (df2 + df1 * x), df1 / 2., df2 / 2., lower_tail, log_p);

		//return !Double.isNaN(x) ? x : Double.NaN;
		return x;
	}

	public static final double quantile(double p, double df1, double df2, boolean lower_tail, boolean log_p)
	{
		if (Double.isNaN(p) || Double.isNaN(df1) || Double.isNaN(df2)) return p + df2 + df1;
		if (df1 <= 0. || df2 <= 0.) return Double.NaN;

		//R_Q_P01_boundaries(p, 0, ML_POSINF);
		if (log_p) {
			if(p > 0)
				return Double.NaN;
			if(p == 0) /* upper bound*/
				return lower_tail ? Double.POSITIVE_INFINITY : 0;
			if(p == Double.NEGATIVE_INFINITY)
				return lower_tail ? 0 : Double.POSITIVE_INFINITY;
		}
		else { /* !log_p */
			if(p < 0 || p > 1)
				return Double.NaN;
			if(p == 0)
				return lower_tail ? 0 : Double.POSITIVE_INFINITY;
			if(p == 1)
				return lower_tail ? Double.POSITIVE_INFINITY : 0;
		}

		/* fudge the extreme DF cases -- qbeta doesn't do this well.
		       But we still need to fudge the infinite ones.
		 */

		if (df1 <= df2 && df2 > 4e5) {
			if(Double.isInfinite(df1)) /* df1 == df2 == Inf : */
				return 1.;
			/* else */
			return ChiSquare.quantile(p, df1, lower_tail, log_p) / df1;
		}
		if (df1 > 4e5) { /* and so  df2 < df1 */
			return df2 / ChiSquare.quantile(p, df2, !lower_tail, log_p);
		}

		p = (1. / Beta.quantile(p, df2/2, df1/2, !lower_tail, log_p) - 1.) * (df2 / df1);
		//return !Double.isNaN(p) ? p : Double.NaN;
		return p;
	}

	public static final double random(double n1, double n2, QRandomEngine random)
	{
	    double v1, v2;
	    if (Double.isNaN(n1) || Double.isNaN(n2) || n1 <= 0. || n2 <= 0.) return Double.NaN;
	    v1 = !Double.isInfinite(n1) ? (ChiSquare.random(n1, random) / n1) : 1;
	    v2 = !Double.isInfinite(n2) ? (ChiSquare.random(n2, random) / n2) : 1;
	    return v1 / v2;
	}
}
