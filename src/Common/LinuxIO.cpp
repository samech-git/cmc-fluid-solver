#ifdef linux
#include "LinuxIO.h"

int LinuxIO::fopen_s( FILE** pFile, const char *filename, const char *mode )
{
	*pFile = fopen(filename, mode);
	
	if (*pFile != NULL or *pFile != 0)
		return 0;
	else   
		return 1;

int LinuxIO::sprintf_s(char *buffer, const char *fmt, ...)
{	
	va_list argptr;
	va_start(argptr,fmt);
	int iout = vsprintf(buffer, fmt, argptr);
	va_end(argptr);
	return iout;
}

char* LinuxIO::_strrev(char* szT)
{
	if ( !szT )                 // handle null passed strings.
		return szT;
	int i = strlen(szT);
	int t = !(i%2)? 1 : 0;      // check the length of the string .
	for(int j = i-1 , k = 0 ; j > (i/2 -t) ; j-- )
	{
		char ch  = szT[j];
		szT[j]   = szT[k];
		szT[k++] = ch;
	}
	return szT;
}

int LinuxIO::_mkdir(const char *dirname)
{
	const int  MY_MASK = 0777;
	return mkdir(dirname, MY_MASK);
}
#endif
