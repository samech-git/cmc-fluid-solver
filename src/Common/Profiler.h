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
#include "..\Common\Timer.h"
#include "..\Common\PARAplan.h"
#elif __unix__
#include "../Common/Timer.h"
#include "../Common/PARAplan.h"
#endif

#include <stdio.h>
#include <string>
#include <map>
#include <algorithm>
#include <vector>

using namespace std;

namespace Common
{
	struct EventInfo
	{
		int count;
		double total_ms;
		double avg_ms;
	};

	class Profiler
	{
	private:
		cpu_timer timer;
		map<string, EventInfo> events;

	public:
		Profiler() { events.clear(); }
		~Profiler() { events.clear(); }

		void StartEvent()
		{
#if PROFILE_ENABLE
#ifdef __PARA
			MPI_Barrier(MPI_COMM_WORLD);
#endif
			timer.start();
#endif
		}

		void StopEvent(const char *name)
		{
#if PROFILE_ENABLE
			timer.stop();
			if( events.find(name) == events.end() )
			{
				// create new
				events[name].count = 1;
				events[name].total_ms = timer.elapsed_ms();
				events[name].avg_ms = events[name].total_ms;
			}
			else
			{
				// update 
				events[name].count++;
				events[name].total_ms += timer.elapsed_ms();
			}
#endif
		}

		struct ValueCmp {
			bool operator()(const pair<string,EventInfo> &lhs, const pair<string,EventInfo> &rhs) {
				return lhs.second.total_ms > rhs.second.total_ms;
			}
		};
		
		void PrintTimings(bool csv)
		{
#if PROFILE_ENABLE
			PARAplan *pplan = PARAplan::Instance();
			if( csv )
			{
				double total_time = 0.0;
				printf("%s,%s,%s,%s,\n", "Event Name", "Total (ms)", "Avg (ms)", "Count");
				
				// copy to vector and sort by values
				vector<pair<string,EventInfo> > v(events.begin(), events.end());
				sort(v.begin(), v.end(), ValueCmp());
				
				for( vector<pair<string,EventInfo> >::iterator it = v.begin(); it != v.end(); it++ )
				{
					EventInfo &e = it->second;
					e.avg_ms = e.total_ms / e.count;
					printf("%s,%.2f,%.2f,%i,\n", it->first.c_str(), e.total_ms, e.avg_ms, e.count); 
					total_time += e.total_ms;
				}
				printf("%s,%.2f,sec\n", "Overall", total_time / 1000);
			}
			else
			{
				printf("Profiling data node(%d):\n", pplan->rank());
				double total_time = 0.0;
				printf("%16s%16s%16s%16s\n", "Event Name", "Total (ms)", "Avg (ms)", "Count");
				
				// copy to vector and sort by values
				vector<pair<string,EventInfo> > v(events.begin(), events.end());
				sort(v.begin(), v.end(), ValueCmp());
				
				for( vector<pair<string,EventInfo> >::iterator it = v.begin(); it != v.end(); it++ )
				{
					EventInfo &e = it->second;
					e.avg_ms = e.total_ms / e.count;
					printf("%16s%16.2f%16.2f%16i\n", it->first.c_str(), e.total_ms, e.avg_ms, e.count); 
					total_time += e.total_ms;
				}
				printf("%16s%16.2f sec\n", "Overall", total_time / 1000);
			}
			fflush(stdout);
#endif
		}
	};
}
