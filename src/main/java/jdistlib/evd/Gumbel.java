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
package jdistlib.evd;

import jdistlib.rng.QRandomEngine;

/**
 * Gumbel distribution. It is a special case of Extreme Value Distribution with shape == 0.
 * Taken from EVD package of R
 *
 */
public class Gumbel {
	public static final double density(double x, double loc, double scale, boolean log)
	{	return GEV.density(x, loc, scale, 0, log); }

	public static final double cumulative(double q, double loc, double scale, boolean lower_tail)
	{	return GEV.cumulative(q, loc, scale, 0, lower_tail); }

	public static final double quantile(double p, double loc, double scale, boolean lower_tail)
	{	return GEV.quantile(p, loc, scale, 0, lower_tail); }

	public static final double random(double loc, double scale, QRandomEngine random)
	{	return GEV.random(loc, scale, 0, random); }
}
