#include "c-target.h"
#include <stdio.h>

#define pp      plist[pnum]
#define mystart pdata[pnum][0][0]
#define myend   pdata[pnum][0][1]
#define mypoint pdata[pnum][0][2]
#define holda   pdata[pnum][0][3]
#define holdb   pdata[pnum][0][4]
#define holdc   pdata[pnum][0][5]
#define MYMODU  pdata[pnum][0][6]
#define val     3
#define jump    4

/*
  This is a simple version of the wild file.
  All of the action takes place on init/next 1 which sets up all of the congruence
  values.
*/


/* Let's START *************************************************************/
int w_file_init(long pnum, long which, long what) {
  long j, thispoint;

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

    MYMODU = (long) intpow((long long) pp, caselist[mystart][0]);
    for(j=3;j<=DEG;j++) {
      pdata[pnum][j][0] = MYMODU;
      /* all other moduli use this value */
    }
    
    pdata[pnum][1][0] = 1; pdata[pnum][1][1] = 0;
    pdata[pnum][2][0] = 1; pdata[pnum][2][1] = 0;
    pdata[pnum][DEG][0]=1; pdata[pnum][DEG][1] = 0; 

    for(j=1;j<=DEG;j++)  /* Initialize everyone */
      pflags[pnum][j] = CALL_NONE;
    pflags[pnum][1] = CALL_INIT;
    pflags[pnum][2] = CALL_INIT;
    pflags[pnum][3] = CALL_BOTH;
    pflags[pnum][DEG] = CALL_INIT; /* Don't do this for degree 3? */
    break;
  case 1: /* Nothing to do */
    break;
  case 2:
    thispoint = mystart;
    what = truemod(A1, MYMODU);
    while((thispoint<myend) && (caselist[thispoint][1]<what))
      thispoint++;
    if((thispoint==myend) || (caselist[thispoint][1] !=what))
      return 0;
    holda = thispoint;
    break;
  case DEG:
    thispoint = holda;
    what = truemod(A2, MYMODU);
    while((thispoint<myend) &&
          (caselist[thispoint][1]==caselist[holda][1])
          && (caselist[thispoint][2]<what))
      thispoint++;
    if((thispoint==myend) || (caselist[thispoint][1] !=caselist[holda][1])
       || (caselist[thispoint][2] !=what))
      return 0;
    holdb = thispoint;
    break;
  case 3: /* This should not be here if deg==3? */
    thispoint = holdb;
    what = truemod(An, MYMODU);
    while((thispoint<myend) &&
          (caselist[thispoint][1]==caselist[holda][1])
          && (caselist[thispoint][2]==caselist[holdb][2])
          && (caselist[thispoint][DEG]<what))
      thispoint++;
    if((thispoint==myend) || (caselist[thispoint][1] !=caselist[holda][1])
       || (caselist[thispoint][2] !=caselist[holdb][2])
       || (caselist[thispoint][DEG] !=what))
      return 0;
    holdc = thispoint;
    mypoint = thispoint;

    for(j=3; j<DEG; j++)
      pdata[pnum][j][1] = caselist[mypoint][j];
    break;
  default:
    printf("*******    Error in init wild file case %ld, p=%ld\n", which,pp);
    break;
  }
  return 1;
}

/* increment counter number pnum modulo p
     returns 1 if ok */
long w_file_next(long pnum, long which) {
  long j;

  if(which != 3) {
    printf("*******    Error in next wild file case %ld, p=%ld\n", which,pp);
    return 1;
  }
  (mypoint)++;
  
  if(mypoint>=myend)
    return 0;

  if((caselist[mypoint][1] !=caselist[holda][1]) ||
     (caselist[mypoint][2] !=caselist[holdb][2]) ||
     (caselist[mypoint][DEG] != caselist[holdc][DEG]))
    return 0;
  
  for(j=3; j<DEG; j++)
    pdata[pnum][j][1] = caselist[mypoint][j];
  /* 
  printf("\nHolding [%ld, %ld,...%ld] ", caselist[holda][1],
         caselist[holdb][2], caselist[holdc][DEG]);
  printf("Using Congruence (%ld, %ld, %ld, %ld, %ld, %ld) ",
         caselist[mypoint][1],caselist[mypoint][2],caselist[mypoint][3],
         caselist[mypoint][4],caselist[mypoint][5],caselist[mypoint][6]);
  */
  
  return 1;
}


/*

We could have just looked at one thing
Inner loops could have moved the pointer forward so we are already on the next good one
They could have also moved us beyond what we can use right now

*/
