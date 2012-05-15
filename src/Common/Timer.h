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

#ifdef _WIN32
#include <windows.h>
#elif linux
#include <time.h>
#endif

namespace Common
{
	struct cpu_timer
	{
#ifdef _WIN32 
		LARGE_INTEGER	m_startTime, m_endTime;
#elif linux
		timespec startTime, endTime;

		timespec diff(timespec& start, timespec& end)
		{
			timespec temp;
			if ((end.tv_nsec - start.tv_nsec) < 0)
			{
				temp.tv_sec = end.tv_sec - start.tv_sec - 1;
				temp.tv_nsec = 1000000000 + end.tv_nsec - start.tv_nsec;
			} else 
			{
				temp.tv_sec = end.tv_sec-start.tv_sec;
				temp.tv_nsec = end.tv_nsec-start.tv_nsec;
			}
			return temp;
		}
#endif

		cpu_timer() { }
		~cpu_timer() { }

#ifdef _WIN32 
		void start() { QueryPerformanceCounter(&m_startTime); }
		void stop()  { QueryPerformanceCounter(&m_endTime); }
#elif linux
		void start() { clock_gettime(CLOCK_REALTIME, &startTime); }
		void stop()  { clock_gettime(CLOCK_REALTIME, &endTime); }
#endif

		float elapsed_sec() {
#ifdef _WIN32 
			LARGE_INTEGER diff;
			LARGE_INTEGER freq;
			QueryPerformanceFrequency(&freq);
			diff.QuadPart = m_endTime.QuadPart - m_startTime.QuadPart;
			return (float)diff.QuadPart / (float)freq.QuadPart;
#elif linux
			timespec temp = diff(startTime, endTime);
			return (float)temp.tv_sec + float(temp.tv_nsec) / 1000000000.0f; 
#endif
		}

		float elapsed_ms()
		{
			return elapsed_sec() * 1000.0f;
		}

	};
}
