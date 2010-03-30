#pragma once

namespace Common
{
	static void SolveTridiagonal(double *a, double *b, double *c, double *d, double *x, int num)
	{
		c[num-1] = 0.0;
		
		c[0] = c[0] / b[0];
		d[0] = d[0] / b[0];

		for (int i = 1; i < num; i++)
		{
			c[i] = c[i] / (b[i] - a[i] * c[i-1]);
			d[i] = (d[i] - d[i-1] * a[i]) / (b[i] - a[i] * c[i-1]);  
		}

		x[num-1] = d[num-1];
	
		for (int i = num-2; i >= 0; i--)
			x[i] = d[i] - c[i] * x[i+1];
	}
}