#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "fscl.h"

static int verbosity_level;

void configure_logmsg(int level) {

  verbosity_level = level;
  return;
}

/* Selective verbosity capabilities -- adapted from SLAN code */
void logmsg(int priority, volatile char *s, ...) {
  va_list ap;
  char long_string[4096];

  if (priority <= verbosity_level) {
    va_start(ap, s);
    vsnprintf(long_string ,4096, (char *) s, ap);
    va_end(ap);
    
    fprintf(stderr, "%s\n", long_string);      
  } 

  if (priority == MSG_FATAL) exit(1);
}

void cr_logmsg(int priority, volatile char *s, ...) {
  va_list ap;
  char long_string[4096];

  if (priority <= verbosity_level) {
    va_start(ap, s);
    vsnprintf(long_string ,4096, (char *) s, ap);
    va_end(ap);
    
    fprintf(stderr, "\r%-79.79s", long_string);      
  } 

  if (priority == MSG_FATAL) exit(1);
}
