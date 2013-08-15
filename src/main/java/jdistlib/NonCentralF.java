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
import jdistlib.rng.QRandomEngine;

public class NonCentralF {
	/*
	 * The density function of the non-central F distribution ---
	 * obtained by differentiating the corresp. cumulative distribution function
	 * using dnbeta.
	 * For df1 < 2, since the F density has a singularity as x -> Inf.
	 */
	public static final double density(double x, double df1, double df2, double ncp, boolean give_log)
	{
		double y, z, f;

		if (Double.isNaN(x) || Double.isNaN(df1) || Double.isNaN(df2) || Double.isNaN(ncp))
			return x + df2 + df1 + ncp;

		/* want to compare dnf(ncp=0) behavior with df() one, hence *NOT* :
		 * if (ncp == 0)
		 *   return df(x, df1, df2, give_log); */

		if (df1 <= 0. || df2 <= 0. || ncp < 0) return Double.NaN;
		if (x < 0.)	 return (give_log ? Double.NEGATIVE_INFINITY : 0.);
		if (Double.isInfinite(ncp)) /* ncp = +Inf -- FIXME?: in some cases, limit exists */
			return Double.NaN;

		/* This is not correct for  df1 == 2, ncp > 0 - and seems unneeded:
		 *  if (x == 0.) return(df1 > 2 ? R_D__0 : (df1 == 2 ? R_D__1 : ML_POSINF));
		 */
		if (Double.isInfinite(df1) && Double.isInfinite(df2)) { /* both +Inf */
			/* PR: not sure about this (taken from  ncp==0)  -- FIXME ? */
			if(x == 1.) return Double.POSITIVE_INFINITY;
			/* else */  return (give_log ? Double.NEGATIVE_INFINITY : 0.);
		}
		if (Double.isInfinite(df2)) /* i.e.  = +Inf */
			return df1* NonCentralChiSquare.density(x*df1, df1, ncp, give_log);
		/*	 ==  dngamma(x, df1/2, 2./df1, ncp, give_log)  -- but that does not exist */
		if (df1 > 1e14 && ncp < 1e7) {
			/* includes df1 == +Inf: code below is inaccurate there */
			f = 1 + ncp/df1; /* assumes  ncp << df1 [ignores 2*ncp^(1/2)/df1*x term] */
			z = Gamma.density(1./x/f, df2/2, 2./df2, give_log);
			return give_log ? z - 2*log(x) - log(f) : z / (x*x) / f;
		}

		y = (df1 / df2) * x;
		z = NonCentralBeta.density(y/(1 + y), df1 / 2., df2 / 2., ncp, give_log);
		return  give_log ?
				z + log(df1) - log(df2) - 2 * log1p(y) :
				z * (df1 / df2) /(1 + y) / (1 + y);
	}

	public static final double cumulative(double x, double df1, double df2, double ncp, boolean lower_tail, boolean log_p)
	{
		double y;
		if (Double.isNaN(x) || Double.isNaN(df1) || Double.isNaN(df2) || Double.isNaN(ncp)) return x + df2 + df1 + ncp;
		if (df1 <= 0. || df2 <= 0. || ncp < 0) return Double.NaN;
		if (Double.isInfinite(ncp)) return Double.NaN;
		if (Double.isInfinite(df1) && Double.isInfinite(df2)) /* both +Inf */
			return Double.NaN;

		//R_P_bounds_01(x, 0., ML_POSINF);
		if(x <= 0) return (lower_tail ? (log_p ? Double.NEGATIVE_INFINITY : 0.) : (log_p ? 0. : 1.));
		if(x >= Double.POSITIVE_INFINITY) return (lower_tail ? (log_p ? 0. : 1.) : (log_p ? Double.NEGATIVE_INFINITY : 0.));

		if (df2 > 1e8) /* avoid problems with +Inf and loss of accuracy */
			return NonCentralChiSquare.cumulative(x * df1, df1, ncp, lower_tail, log_p);

		y = (df1 / df2) * x;
		return NonCentralBeta.pnbeta2(y/(1. + y), 1./(1. + y), df1 / 2., df2 / 2., ncp, lower_tail, log_p);
	}

	public static final double quantile(double p, double df1, double df2, double ncp, boolean lower_tail,  boolean log_p)
	{
		double y;
		if (Double.isNaN(p) || Double.isNaN(df1) || Double.isNaN(df2) || Double.isNaN(ncp)) return p + df1 + df2 + ncp;
		if (df1 <= 0. || df2 <= 0. || ncp < 0) return Double.NaN;
		if (Double.isInfinite(ncp)) return Double.NaN;
		if (Double.isInfinite(df1) && Double.isInfinite(df2)) return Double.NaN;
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

		if (df2 > 1e8) /* avoid problems with +Inf and loss of accuracy */
			return NonCentralChiSquare.quantile(p, df1, ncp, lower_tail, log_p)/df1;
		y = NonCentralBeta.quantile(p, df1 / 2., df2 / 2., ncp, lower_tail, log_p);
		return y/(1-y) * (df2/df1);
	}

	public static final double random(double df1, double df2, double ncp, QRandomEngine random)
	{
		if (ncp == 0)
			return F.random(df1, df2, random);
		return (NonCentralChiSquare.random(df1, ncp, random) / df1) / (NonCentralChiSquare.random(df2, ncp, random) / df2);
	}
}
