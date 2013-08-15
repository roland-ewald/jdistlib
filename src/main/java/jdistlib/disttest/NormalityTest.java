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

import java.util.HashSet;
import java.util.Set;

import jdistlib.ChiSquare;
import jdistlib.Normal;

import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.asin;
import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.signum;
import static java.lang.Math.sqrt;
import static java.util.Arrays.sort;
import static jdistlib.disttest.Utils.calculate_ecdf;

/** 
 * A package about normality testing. Taken from various places.
 * @author Roby Joehanes
 */
public class NormalityTest {
	public static final double shapiro_wilk_statistic(double[] X)
	{
		final double kSqrtHalf = sqrt(0.5);

		// constant for Shapiro-wilk
		final double[]
			c1 = {0, 0.221157, -0.147981, -2.07119, 4.434685, -2.706056},
			c2 = {0, 0.042981, -0.293762, -1.752461, 5.682633, -3.582633};
		int	n = X.length, n2 = n/2; // yes, integer division
		if (n < 3)
			return 0;
		double[] a = new double[n2];
		if (n == 3)
			a[0] = kSqrtHalf;
		else
		{
			double an25 = n + 0.25, sum2 = 0, sqrtSum2;
			for (int i = 0; i < n2; i++)
			{
				double val = a[i] = Normal.quantile((i + 0.625) / an25, 0, 1, true, false);
				sum2 += val * val;
			}
			sum2 *= 2;
			sqrtSum2 = sqrt(sum2);
			double rsn = 1 / sqrt(n), a1 = poly(c1, rsn) - a[0] / sqrtSum2, fac, a2;
			int i1 = 1;
			if (n > 5)
			{
				i1 = 2;
				a2 = -a[1] / sqrtSum2 + poly(c2, rsn);
				fac = -sqrt( (sum2 - 2 * a[0] * a[0] - 2 * a[1] * a[1]) / (1 - 2 * a1 * a1 - 2 * a2 * a2) );
				a[1] = a2;
			}
			else
				fac = -sqrt( (sum2 - 2 * a[0] * a[0]) / (1 - 2 * a1 * a1) );
			a[0] = a1;
			for (int i = i1; i < n2; i++)
				a[i] /= fac;
		}

		double range = X[n - 1] - X[0], xx = X[0] / range, sx = xx, sa = -a[0], xi;
		int j = n - 1;
		for (int i = 1; i < n; j--)
		{
			xi = X[i] / range;
			sx += xi;
			i++;
			if (i != j)
				sa += signum(i - j) * a[min(i,j) - 1];
			xx = xi;
		}

		// W statistic is a squared correlation between data and coefficients
		double ssa = 0, ssx = 0, sax = 0;
		sa /= n; sx /= n; j = n - 1;
		for (int i = 0; i < n; ++i, --j)
		{
			double
			asa = i != j ? signum(i - j) * a[min(i,j)] - sa : -sa,
				xsx = X[i] / range - sx;
			ssa += asa * asa;
			ssx += xsx * xsx;
			sax += asa * xsx;
		}

		// W1 = (1-W) calculated to avoid excessive rounding error for W
		// very near 1 (a potential problem in very large samples)

		double
			ssassx = sqrt(ssa * ssx),
			w1 = (ssassx - sax) * (ssassx + sax) / (ssa * ssx);

		return 1 - w1;
	}

	/**
	 * 
	 * @param w
	 * @param n The length of the array
	 * @return
	 */
	public static final double shapiro_wilk_pvalue(double w, int n)
	{
		final double kVerySmallValue = 1e-19;

		// constant for Shapiro-wilk
		final double[]
			c3 = {0.544, -0.39978, 0.025054, -6.714e-4},
			c4 = {1.3822, -0.77857, 0.062767, -0.0020322},
			c5 = {-1.5861, -0.31082, -0.083751, 0.0038915},
			c6 = {-0.4803, -0.082676f, 0.0030302},
			g = {-2.273, 0.459};

		if (n < 3)
			return 1;
		if (n == 3) // exact P value :
			return max(0, 1.90985931710274 * (asin(sqrt(w)) - 1.04719755119660));

		double
		y = log(1 - w),
		xx = log(n),
		gamma,
		m, s;

		if (n <= 11)
		{
			gamma = poly(g, n);
			if (y >= gamma)
				return kVerySmallValue; // rather use an even smaller value, or NA ?
			y = -log(gamma - y);
			m = poly(c3, n);
			s = exp(poly(c4, n));
		}
		else
		{
			m = poly(c5, xx);
			s = exp(poly(c6, xx));
		}

		return Normal.cumulative(y, m, s*s, false, false);
	}

	public static final double anderson_darling_statistic(double[] X)
	{
		int n = X.length;
		double
			sum = 0,
			sumSq = 0;
	
		for (int i = 0; i < n; i++)
		{
			double value = X[i];
			sum += value;
			sumSq += value * value;
		}
	
		double
			mean = sum / n,
			sd = sqrt((sumSq - sum * mean) / (n - 1)),
			Y[] = new double[n];
	
		// Standardize X
		for (int i = 0; i < n; i++)
			Y[i] = (X[i] - mean) / sd;
		sort(Y);
	
		// Get corresponding normal CDF
		for (int i = 0; i < n; i++)
			Y[i] = Normal.cumulative(Y[i], 0, 1, true, false);
	
		sum = 0;
		for (int i = 1; i <= n; i++)
			sum += (2 * i - 1) * (log(Y[i-1]) + log(1 - Y[n-i]));
	
		return - n - sum / n;
	}

	/**
	 * 
	 * @param value
	 * @param n The length of the array
	 * @return
	 */
	public static final double anderson_darling_pvalue(double value, int n)
	{
		double
			aa = value * (1 + 0.75/n + 2.25 / (n*n)),
			aasq = aa * aa;
		if (aa < 0.2)
			return 1 - exp(-13.436 + 101.14 * aa - 223.73 * aasq);
		else if (aa < 0.34)
			return 1 - exp(-8.318 + 42.796 * aa - 59.938 * aasq);
		else if (aa < 0.6)
			return exp(0.9177 - 4.279 * aa - 1.38 * aasq);
		return exp(1.2937 - 5.709 * aa + 0.0186 * aasq);
	}

	public static final double cramer_vonmises_statistic(double[] X)
	{
		int
			n = X.length,
			n2 = n * 2;
		double
			w = 0,
			sum = 0,
			sumSq = 0;

		for (int i = 0; i < n; i++)
		{
			double val = X[i];
			sum += val;
			sumSq += val * val;
		}

		double
			mean = sum / n,
			sd = sqrt((sumSq - sum * mean) / (n - 1)),
			sortedZ[] = new double[n];

		// Standardize X
		for (int i = 0; i < n; i++)
			sortedZ[i] = Normal.cumulative(X[i], mean, sd, true, false);
		sort(sortedZ);

		w = 1.0 / (12 * n);
		for (int i = 0; i < n; i++)
		{
			double val = (2 * i + 1.0) / n2 - sortedZ[i];
			w += val * val;
		}
		return w;
	}

	/**
	 * 
	 * @param w
	 * @param n The length of the array
	 * @return
	 */
	public static final double cramer_vonmises_pvalue(double w, int n)
	{
		double
			ww = (1 + 0.5/n) * w,
			ww2 = ww * ww;
	    if (ww < 0.0275)
	        return 1 - exp(-13.953 + 775.5 * ww - 12542.61 * ww2);
	    else if (ww < 0.051)
	        return 1 - exp(-5.903 + 179.546 * ww - 1515.29 * ww2);
	    else if (ww < 0.092)
	        return exp(0.886 - 31.62 * ww + 10.897 * ww2);
	    ww = exp(1.111 - 34.242 * ww + 12.832 * ww2);
		return ww > 1 ? 0 : ww;
	}

	/**
	 * Calculate D'Agostino-Pearson test for normality. Follows Chi^2 distribution with df = 2.
	 * This is the infamous omnibus test. Note: It breaks down for n <= 7. See:<br>
	 * Doornik, Hansen, An Omnibus Test for Univariate and Multivariate Normality, 1994
	 * 
	 * @param X
	 * @return
	 */
	public static final double dagostino_pearson_statistic(double[] X)
	{
		// These are all magic numbers I took from:
		// Handbook of Parametric and Non-Parametric Statistical Procedures by David Sheskin (3rd ed.)
		int
			n = X.length,
			nSq = n * n,
			nCube = n * nSq;
		double
			nMin1 = n - 1,
			nn1 = n * nMin1,
			n1n3 = (n + 1) * (n + 3),
			nMin2 = n - 2,
			n2n3 = nMin2 * (n - 3),
			n3n5 = (n + 3) * (n + 5),
			n7n9 = (n + 7) * (n + 9),
			sum = 0,
			sumSq = 0,
			sumCube = 0,
			sumQuad = 0;

		for (int i = 0; i < n; i++)
		{
			double
				val = X[i],
				valsq = val * val;
			sum += val;
			sumSq += valsq;
			sumCube += valsq * val; 
			sumQuad += valsq * valsq;
		}

		double
			sumsum = sum * sum,
			kurtosis = ((nCube + nSq) * sumQuad - 4 * (nSq + n) * sumCube * sum - 3 * (nSq - n) * sumSq * sumSq
				+ 12 * n * sumSq * sumsum - 6 * sumsum * sumsum) / (nn1 * n2n3),
			skewness = (n * sumCube - 3 * sum * sumSq + (2 * sumsum * sum / n)) / (nMin1 * nMin2),
			variance = (sumSq - sumsum / n) / nMin1,
			c = sqrt(2 * ((3 * (nSq + 27 * n - 70) * n1n3) / (nMin2 * (n + 5) * n7n9)) - 1) - 1,
			f = (nMin2 * skewness / (variance * sqrt(variance) * sqrt(nn1))) * sqrt(n1n3 * (c - 1) / (12 * nMin2)),
			jinv = n7n9 / ((6 * nSq - 30 * n + 12) * sqrt((6 * n3n5) / (n * n2n3))),
			k = 6 + 8 * jinv * (2 * jinv + sqrt(1 + 4 * jinv * jinv)),
			l = (1 - 2/k) / (1 + (sqrt( n2n3 * n3n5 / (24 * n)) * abs(kurtosis / (variance * variance)) / nMin1) * sqrt(2 / (k - 4))),
			k2 = 2 / (9 * k),
			z1 = log(f + sqrt(f * f + 1)) / sqrt(0.5 * log(c)),
			z2 = (1 - k2 - pow(l, 1/3.0)) / sqrt(k2);

		double value = z1 * z1 + z2 * z2;
		return value > 50 ? 50 : value;
	}

	public static final double dagostino_pearson_pvalue(double value)
	{	return ChiSquare.cumulative(value, 2, true, false); }

	/**
	 * Calculate Jarque-Bera Normality Test. Follows Chi^2 distribution with df = 2<br>
	 * http://en.wikipedia.org/wiki/Jarque-Bera_test
	 * 
	 * @param X
	 * @return
	 */
	public static final double jarque_bera_statistic(double[] X)
	{
		// These are all magic numbers I took from:
		// http://en.wikipedia.org/wiki/Jarque-Bera_test
		int n = X.length;
		double
			sum = 0,
			sumSq = 0,
			sumCube = 0,
			sumQuad = 0;
	
		for (int i = 0; i < n; i++)
		{
			double val = X[i];
			sum += val;
			double valsq = val * val;
			sumSq += valsq;
			sumCube += valsq * val; 
			sumQuad += valsq * valsq;
		}
	
		int
			nSq = n * n,
			nCube = n * nSq;
		double
			sumsum = sum * sum,
			variance = (sumSq - sum * sum / n) / n,
			skewness = sumCube / n - 3 * sumSq * sum / nSq + 2 * sumsum * sum / nCube,
			kurtosis = sumQuad / n - 4 * sumCube * sum / nSq + 6 * sumSq * sumsum / nCube
				- 3 * sumsum * sumsum / (nSq * nSq);
		skewness = skewness / (variance * sqrt(variance));
		kurtosis = kurtosis / (variance * variance);
		double kMin3 = kurtosis - 3;
	
		double value = n * (skewness * skewness + kMin3 * kMin3 / 4) / 6;
		return value > 50 ? 50 : value; // Cap it 
	}

	public static final double jarque_bera_pvalue(double value)
	{	return ChiSquare.cumulative(value, 2, true, false); }

	public static final double kolmogorov_smirnov_statistic(double[] X)
	{
		int n = X.length;
		double
			sum = 0,
			sumSq = 0;
	
		for (int i = 0; i < n; i++)
		{
			double val = X[i];
			sum += val;
			sumSq += val * val;
		}
	
		double
			mean = sum / n,
			sd = sqrt((sumSq - sum * mean) / (n - 1)),
			sortedZ[] = new double[n];
	
		// Standardize X
		for (int i = 0; i < n; i++)
			sortedZ[i] = (X[i] - mean) / sd;
		sort(sortedZ);
	
		double
			cdfZ[] = calculate_ecdf(sortedZ),
			max = 0;
	
		for (int i = 0; i < n; i++)
		{
			double
				cdfNormal = Normal.cumulative(sortedZ[i], 0, 1, true, false),
				M = abs(cdfNormal - cdfZ[i]),
				MPrime = i > 0 ? abs(cdfNormal - cdfZ[i-1]) : abs(cdfNormal),
				curMax = max(M, MPrime);
			if (curMax > max)
				max = curMax;
		}
		return max;
	}

	/**
	 * 
	 * @param d
	 * @param X The original array with which you invoked kolmogorov_smirnov_statistic
	 * @return
	 */
	public static final double kolmogorov_smirnov_pvalue(double d, double[] X)
	{
		int n = X.length;
		Set<Double> set = new HashSet<Double>(n);
		for (double x: X)
			set.add(x);
		if (set.size() < n)
			return kolmogorov_smirnov_pvalue_with_ties(d, X);
		int
			k = (int) (n * d) + 1,
			m = 2 * k - 1,
			mm = m * m;
		double
			h = k - n * d,
			H[] = new double[mm],
			Q[] = new double[mm];
		for (int i = 0; i < m; i++)
			for (int j = 0; j < m; j++)
				H[i * m + j] = i - j + 1 < 0 ? 0 : 1;

		for(int i = 0; i < m; i++)
		{
			H[i * m] -= pow(h, i + 1);
			H[(m - 1) * m + i] -= pow(h, (m - i));
		}
		H[(m - 1) * m] += ((2 * h - 1 > 0) ? pow(2 * h - 1, m) : 0);
		for (int i = 0; i < m; i++)
			for (int j = 0; j < m; j++)
				if(i - j + 1 > 0)
					for(int g = 1; g <= i - j + 1; g++)
						H[i * m + j] /= g;
		int eH = 0;
		double eQ = m_power(H, eH, Q, 0, m, n);
		double s = Q[(k - 1) * m + k - 1];
		for(int i = 1; i <= n; i++)
		{
			s = s * i / n;
			if(s < 1e-140)
			{
				s *= 1e140;
				eQ -= 140;
			}
		}
		s *= pow(10., eQ);
		return 1 - s;
	}

	/**
	 * Calculate an approximation of P-Value of Kolmogorov-Smirnov Normality test.
	 * Don't call this routine unless you're absolutely sure that some elements in X
	 * tie with some others.
	 * 
	 * @param d
	 * @param X
	 * @return
	 */
	private static final double kolmogorov_smirnov_pvalue_with_ties(double d, double[] X)
	{
		final double
			kKSTolerance = 1e-6,
			kPiSqDiv8 = PI * PI / 8,
			k1Sqrt2Pi = 1 / sqrt(2 * PI);
		int
			n = X.length,
			k_max = (int) sqrt(2 - log(kKSTolerance));
		d = sqrt(n) * d;

		double z, w, s;
		if (d < 1)
		{
			z = -kPiSqDiv8 / (d * d);
			w = log(d);
			s = 0;
			for (int k = 1; k < k_max; k += 2)
				s += exp(k * k * z - w);
			return 1 - s / k1Sqrt2Pi;
		}

		z = -2 * d * d;
		s = -1;
		int k = 1;
		double
			old = 0,
			neww = 1;
		while (abs(old - neww) > kKSTolerance)
		{
			old = neww;
			neww += 2 * s * exp(z * k * k);
			s *= -1;
			k++;
		}
		return 1 - neww;
	}

	/**
	 * Exactly identical as kolmogorov_smirnov_statistic
	 * @param X
	 * @return
	 */
	public static final double kolmogorov_lilliefors_statistic(double[] X)
	{	return kolmogorov_smirnov_statistic(X); }

	/**
	 * 
	 * @param k
	 * @param n The length of the array
	 * @return
	 */
	public static final double kolmogorov_lilliefors_pvalue(double k, int n)
	{
		double
			Kd = k,
			nd = n;
		if (n > 100)
		{
			Kd = k * pow((n/100.0),0.49);
			nd = 100;
		}

		double pvalue = exp(-7.01256 * Kd * Kd * (nd + 2.78019) + 2.99587 * Kd * sqrt(nd + 2.78019) - 0.122119 + 0.974598/sqrt(nd) + 1.67997/nd);
		if (pvalue > 0.1)
		{
			double
				KK = (sqrt(n) - 0.01 + 0.85/sqrt(n)) * k,
				KK2 = KK * KK,
				KK3 = KK2 * KK,
				KK4 = KK2 * KK2;
			if (KK <= 0.302)
				return 1;
			else if (KK <= 0.5)
				return 2.76773 - 19.828315 * KK + 80.709644 * KK2 - 138.55152 * KK3 + 81.218052 * KK4;
			else if (KK <= 0.9)
				return -4.901232 + 40.662806 * KK - 97.490286 * KK2 + 94.029866 * KK3 - 32.355711 * KK4;
			else if (KK <= 1.31)
				return 6.198765 - 19.558097 * KK + 23.186922 * KK2 - 12.234627 * KK3 + 2.423045 * KK4;
			return 0;
		}
		return pvalue; 
	}

	/**
	 * Shapiro-Francia normality test
	 * @param X a sorted array of values
	 * @return
	 */
	public static final double shapiro_francia_statistic(double[] X)
	{
		int n = X.length;
		double
			denum = n + 0.25,
			sumX = 0,
			sumY = 0,
			sumXSq = 0,
			sumYSq = 0,
			sumXY = 0,
			x,
			y;

		for (int i = 0; i < n; i++)
		{
			y = Normal.quantile((i + 0.625) / denum, 0, 1, true, false);
			x = X[i];
			sumX += x;
			sumY += y;
			sumXSq += x * x;
			sumYSq += y * y;
			sumXY += x * y;
		}
		double cor = (n * sumXY - sumX * sumY) / (sqrt(n * sumXSq - sumX * sumX) * sqrt(n * sumYSq - sumY * sumY));
		return cor * cor;
	}

	/**
	 * P-value of Shapiro-Francia normality test
	 * @param w the result from ShapiroFrancia's statistic
	 * @param n the length of the original array
	 * @return
	 */
	public static final double shapiro_francia_pvalue(double w, int n)
	{
		double
			a = log(n),
			b = log(a),
			mu = -1.2725 + 1.0521 * (b - a),
			sigma = 1.0308 - 0.26758 * (b + 2/a);
		return Normal.cumulative(log(1 - w), mu, sigma, false, false);
	}

	/**
	 * Helper function to calculate polynomials (for Shapiro-Wilk)
	 * @param coeff
	 * @param x
	 * @return
	 */
	private static final double poly(double[] coeff, double x)
	{
		int n = coeff.length;
		double result = coeff[0];
		if (n > 1)
		{
			double p = x * coeff[n - 1];
			for (int i = n - 2; i > 0; i--)
				p = (p + coeff[i]) * x;
			result += p;
		}
		return result;
	}

	/**
	 * Helper function for Kolmogorov-Smirnov
	 * @param A
	 * @param eA
	 * @param V
	 * @param eV
	 * @param m
	 * @param n
	 * @return
	 */
	private static final int m_power(double[] A, int eA, double[] V, int eV, int m, int n)
	{
		double[] B = new double[m * m];
		int eB;

		if(n == 1)
		{
			for (int i = 0; i < m * m; i++)
				V[i] = A[i];
			return eA;
		}
		eV = m_power(A, eA, V, eV, m, n / 2);
		m_multiply(V, V, B, m);
		eB = 2 * eV;

		if((n & 1) == 0)
		{
			for (int i = 0; i < m * m; i++)
				V[i] = B[i];
			eV = eB;
		}
		else
	    {
			m_multiply(A, B, V, m);
			eV = eA + eB;
	    }

		int mdiv2 = m / 2;
		if(V[mdiv2 * m + mdiv2] > 1e140)
		{
			for (int i = 0; i < m * m; i++)
				V[i] = V[i] * 1e-140;
			eV += 140;
		}
		return eV;
	}

	/**
	 * Helper function for Kolmogorov-Smirnov
	 * @param A
	 * @param B
	 * @param C
	 * @param m
	 */
	private static final void m_multiply(double[] A, double[] B, double[] C, int m)
	{
		for (int i = 0; i < m; i++)
			for (int j = 0; j < m; j++)
			{
				double s = 0.0;
				for (int k = 0; k < m; k++)
					s+= A[i * m + k] * B[k * m + j];
				C[i * m + j] = s;
			}
	}

}
