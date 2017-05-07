#ifndef _PLATFORM_H
#define _PLATFORM_H

#include <ctime>

#if defined(__WIN32__) /* Borland */ || \
  defined(_WIN32) /* MSC et al */

#define WIN32

#if !defined(_MSC_VER) || defined(_MSC_EXTENSIONS)

# pragma warning(push, 0)
# include <windows.h>
# pragma warning(pop)
#undef min
#undef max

#else

  // N.B. Including <windows.h> disallows us to compile with /Za
extern "C" {
  unsigned __stdcall GetCurrentProcessId(void);
  unsigned __stdcall GetTickCount(void);
}

#endif // !(_MSC_VER && _MSC_EXTENSIONS)

namespace {

inline static clock_t
getseed()
{ return clock_t(GetTickCount() ^ (GetCurrentProcessId() << 8)); }

} // <anonymous>

#elif defined(__APPLE__) && defined(__MACH__) /* Mac OS X */ || \
  defined(__unix) || \
  defined(__unix__) || \
  defined(unix) /* Most Unices define one of the three. */

# include <unistd.h>
# include <sys/time.h>

namespace {

inline static clock_t
getseed()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return clock_t((tv.tv_usec / 1000UL) | (unsigned(getpid()) << 8));
}

} // <anonymous>

#else
# error: Unsupported platform.
#endif

#endif // _PLATFORM_H
