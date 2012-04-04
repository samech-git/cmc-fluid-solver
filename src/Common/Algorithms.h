/*
 *  Copyright 2010-2011 Nikolai Sakharnykh
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

#pragma once

namespace Common
{
	static void SolveTridiagonal( FTYPE *a, FTYPE *b, FTYPE *c, FTYPE *d, FTYPE *x, int num )
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
