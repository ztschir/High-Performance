/*
  This routine has been furnished by Kazushige Goto (TACC, UT-Austin).
  We use it with his permission.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>


// This function is only used internally, which is why its prototype is not
// with the other "util" directory function prototypes.
void   detect_clocks( void );

// A global variable used when dclock() is defined in terms of gettimeofday().
double gtod_ref_time_sec = 0.0;


// --- Begin portable dclock() definitions -------------------------------------
#ifdef FLA_ENABLE_PORTABLE_DCLOCK

double dclock()
{
  double the_time, norm_sec;
  struct timeval tv;

  gettimeofday( &tv, NULL );

  // If this is the first invocation of through dclock(), then initialize the 
  // "reference time" global variable to the seconds field of the tv struct.
  if( gtod_ref_time_sec == 0.0 )
    gtod_ref_time_sec = ( double ) tv.tv_sec;

  // Normalize the seconds field of the tv struct so that it is relative to the
  // "reference time" that was recorded during the first invocation of dclock().
  norm_sec = ( double ) tv.tv_sec - gtod_ref_time_sec;

  // Compute the number of seconds since the reference time.
  the_time = norm_sec + tv.tv_usec * 1.0e-6;

  return the_time;
}

#else
// --- Begin possibly non-portable dclock() definitions ------------------------


// Global variables that are used only for non-portable dclock() definitions.
static double        clocks             = 0.0;
#ifdef __i386__
static unsigned int  initialclockoffset = 0;
#endif

// --- Begin ia64 section ------------------------------------------------------
#ifdef __ia64__


static __inline unsigned long rdtsc( void );
static __inline unsigned long rdtsc()
{
  unsigned long clocks;

#ifdef __INTEL_COMPILER
  clocks = __getReg(_IA64_REG_AR_ITC);  
#else
  __asm__ __volatile__ ("mov %0=ar.itc" : "=r"(clocks));
#endif

  return clocks;
}


double dclock()
{
  unsigned long totalclocks;

  if ( clocks == 0.0 ) detect_clocks();
  totalclocks = rdtsc();

  return (double)totalclocks / clocks;
}


// --- End ia64 section --------------------------------------------------------
#else
// --- Begin i386 section ------------------------------------------------------
#ifdef __i386__

// rdtsc was once declared static. neccessary?
static inline void rdtsc( unsigned int *high, unsigned int *low )
{
  asm("rdtsc" : "=a" (*low), "=d"(*high): : "cc");
}


double dclock()
{
  unsigned int high, low;
  unsigned long long totalclocks;

  if (!clocks) detect_clocks();
  rdtsc(&high, &low);
  high -= initialclockoffset;
  totalclocks = (((unsigned long long)high) << 32)
                | (unsigned long long)low;

  return (double)totalclocks / clocks;
}


// --- End i386 section --------------------------------------------------------
#else
// --- Begin non-ia64, non-i386 section ----------------------------------------


double dclock()
{
  double the_time, norm_sec;
  struct timeval tv;

  gettimeofday( &tv, NULL );

  // If this is the first invocation of through dclock(), then initialize the 
  // "reference time" global variable to the seconds field of the tv struct.
  if( gtod_ref_time_sec == 0.0 )
    gtod_ref_time_sec = ( double ) tv.tv_sec;

  // Normalize the seconds field of the tv struct so that it is relative to the
  // "reference time" that was recorded during the first invocation of dclock().
  norm_sec = ( double ) tv.tv_sec - gtod_ref_time_sec;

  // Compute the number of seconds since the reference time.
  the_time = norm_sec + tv.tv_usec * 1.0e-6;

  return the_time;
}

#endif
#endif
// --- End non-ia64, non-i386 section ------------------------------------------
// --- Return to possibly non-portable dclock() definitions --------------------


// This implementation of detect_clocks() works only on operating systems that
// support the /proc filesystem (in particular, GNU/Linux).
// We place this function here at the end so that rdtsc() has been defined
// (and prototyped) by the time we get here. Otherwise the compiler spits out
// warnings.
void detect_clocks()
{
  FILE *infile;
  char buffer[256], *p;
#ifdef __i386__
  unsigned int high, low;
#endif

  if ( clocks == 0.0 ){
    p = (char *)NULL;
    infile = fopen("/proc/cpuinfo", "r");
    while (fgets(buffer, sizeof(buffer), infile)){
      if (!strncmp("cpu MHz", buffer, 6)){
	p = strchr(buffer, ':') + 1;
	break;
      }
    }
    clocks = 1.e6 * atof(p);
#ifdef __i386__
    rdtsc(&high, &low);
    initialclockoffset = high;
#endif
  }
}

// --- End possibly non-portable dclock() definitions --------------------------
#endif
