#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include <pthread.h>

#include "fscl.h"

static int verbosity_level;
static pthread_mutex_t logmsg_lock;

void configure_logmsg(int level) {

  pthread_mutex_init(&logmsg_lock, NULL);
  verbosity_level = level;
  return;
}

/* Selective verbosity capabilities -- adapted from SLAN code */
void logmsg(int priority, volatile char *s, ...) {
  va_list ap;
  char long_string[MAX_LOGMSG_LENGTH];

  pthread_mutex_lock(&logmsg_lock);
  if (priority <= verbosity_level) {
    va_start(ap, s);
    vsnprintf(long_string ,MAX_LOGMSG_LENGTH, (char *) s, ap);
    va_end(ap);
    
    fprintf(stderr, "%s\n", long_string);      
  } 
  pthread_mutex_unlock(&logmsg_lock);

  if (priority == MSG_FATAL) exit(1);
}

void cr_logmsg(int priority, volatile char *s, ...) {
  va_list ap;
  char long_string[MAX_LOGMSG_LENGTH];

  pthread_mutex_lock(&logmsg_lock);
  if (priority <= verbosity_level) {
    va_start(ap, s);
    vsnprintf(long_string , MAX_LOGMSG_LENGTH, (char *) s, ap);
    va_end(ap);
    
    fprintf(stderr, "\33[2K\r%-79.79s", long_string);
  } 
  pthread_mutex_unlock(&logmsg_lock);

  if (priority == MSG_FATAL) exit(1);
}
