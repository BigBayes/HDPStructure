/* libnp
 * Copyright (c) 2013, Lloyd T. Elliott and Yee Whye Teh
 */

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static java.lang.Math.pow;

public class crpkprb {
	public static double[] crp_sizes(double alpha, double d, int n) {

		if (n == 0) {
			return new double[] {};
		} else if (n == 1) {
			return new double[] { 1.0 };
		}
		double[] pdf = new double[n];
		double[] T = new double[n - 1];
		for (int i = 1; i < n; i++) {
			T[i - 1] = 0.0;
			for (int j = 1; j < n; j++) {
				T[i - 1] += pow(alpha / j, i);
			}
		}

		pdf[0] = 1.0;
		for (int i = 1; i < n; i++) {
			pdf[0] *= 1.0 - alpha / (i + alpha);
		}
		for (int k = 1; k < n; k++) {
			double p = 0.0;
			int sign = 1;
			for (int i = 1; i <= k; i++) {
				p += sign * pdf[k - i] * T[i - 1];
				sign *= -1;
			}
			pdf[k] = p / k;
		}

		for (int i = 1; i < n; i++) {
			if (Double.isNaN(pdf[i]) || pdf[i] < 0.0) {
				pdf[i] = 0.0;
			}
		}
		return pdf;
	}
}
