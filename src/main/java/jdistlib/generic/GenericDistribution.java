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
package jdistlib.generic;

import jdistlib.rng.QRandomEngine;

/**
 * An interface for a generic distribution. All parameters have to be encoded (either as fields or otherwise).
 * Treat this interface as an adapter to the other distributions.
 * 
 * @author Roby Joehanes
 *
 */
public interface GenericDistribution {
	public double density(double x, boolean log);
	public double cumulative(double p, boolean lower_tail);
	public double quantile(double q, boolean lower_tail);
	public double random(QRandomEngine random);
}
