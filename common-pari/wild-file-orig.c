#include "c-target.h"
#include <stdio.h>

#define pp      plist[pnum]
#define mystart pdata[pnum][0][0]
#define myend   pdata[pnum][0][1]
#define mypoint pdata[pnum][0][2]
#define val     3

/*
  This is a simple version of the wild file.
  All of the action takes place on init/next 1 which sets up all of the congruence
  values.
*/


/* Let's START *************************************************************/
int w_file_init(long pnum, long which, long what) {
  long j;

  /*
    We have to give the modulus and the value
     Input of which = 0 says to return for a_1
  */
  /*
    This file is unusual in that when "which" is 0, "what" tells us which global case
    we are accessing.
  */
  switch(which) {
  case 0:
    mystart = cases[pnum][what].start;
    myend = cases[pnum][what].end;

    for(j=1;j<=DEG;j++) {
      pdata[pnum][j][0] = (long) intpow((long long) pp, caselist[mystart][0]);
      /* all moduli use this value */
    }
    pflags[pnum][1] = 0;
    for(j=2;j<=DEG;j++) {
      pflags[pnum][j] = 1; /* Don't call these at all */
    }
    break;
  case 1:
    mypoint = mystart;
    for(j=1; j<=DEG; j++)
      pdata[pnum][j][1] = caselist[mypoint][j];
    mypoint++;
    break;
  case 2: /* No one else should ever be called. */
  default:
    printf("*******    Error in wild file case %ld, p=%ld\n", which,pp);
    break;
  }
  return 1;
}

/* increment counter number pnum modulo p
     returns 1 if ok */
long w_file_next(long pnum, long which) {
  long j;
  
  if(mypoint>=myend || which != 1)
    return 0;
  
  for(j=1; j<=DEG; j++)
    pdata[pnum][j][1] = caselist[mypoint][j];
  mypoint++;
  
  return 1;
}


/*

We could have just looked at one thing
Inner loops could have moved the pointer forward so we are already on the next good one
They could have also moved us beyond what we can use right now

*/
