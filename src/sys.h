#include <sys/resource.h>
#include <sys/time.h>

// Stolen from https://github.com/lh3/miniprot/blob/69487c93a7680dbfa6bfe2e015645246f8728809/sys.c#L96
double mp_cputime(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

long mp_peakrss(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
#ifdef __linux__
	return r.ru_maxrss * 1024;
#else
	return r.ru_maxrss;
#endif
}