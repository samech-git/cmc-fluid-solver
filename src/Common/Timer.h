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