#include <stdio.h>
#include <stdarg.h>

// Overload printf to redirect outputs to Scilab
#include <MALLOC.h>
#include <sciprint.h>

#define BUFFER_SIZE 4096

#ifdef WIN32
extern "C" int printf(const char * format,...)
#else
// this function is used with linker option -Wl,-wrap,printf to overload printf
extern "C" int __wrap_printf (char *format, ...)
#endif
{
  char * BUFFER = (char *)MALLOC(BUFFER_SIZE*sizeof(char));
  va_list args;
  va_start(args, format);
  vsprintf(BUFFER, format, args);
  sciprint("SYMPHONY: %s", BUFFER);
  FREE(BUFFER);
  return 0;
}
