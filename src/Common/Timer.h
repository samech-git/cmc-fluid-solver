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

#include <windows.h>

namespace Common
{
	struct cpu_timer
	{
		LARGE_INTEGER	m_startTime, m_endTime;

		cpu_timer() { }
		~cpu_timer() { }

		void start() { QueryPerformanceCounter(&m_startTime); }
		void stop()  { QueryPerformanceCounter(&m_endTime); }

		float elapsed_sec() { 
			LARGE_INTEGER diff;
			LARGE_INTEGER freq;
			QueryPerformanceFrequency(&freq);
			diff.QuadPart = m_endTime.QuadPart - m_startTime.QuadPart;
			return (float)diff.QuadPart / (float)freq.QuadPart;
		}

		float elapsed_ms()
		{
			return elapsed_sec() * 1000.0f;
		}

	};
}