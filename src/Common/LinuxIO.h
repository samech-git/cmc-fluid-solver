#pragma once

#include <stdio.h>
#include <string.h>
#include <cstdarg>
#include <sys/types.h>
#include <sys/stat.h>

using namespace std;

#define fscanf_s  fscanf
#define sprintf_s snprintf
#define strcpy_s  strcpy

namespace LinuxIO
{
	extern int  fopen_s( FILE**, const char*, const char*);
	extern int  sprintf_s(char*, const char*, ...);
	extern char* _strrev(char*);
	extern int   _mkdir(const char*);
}
