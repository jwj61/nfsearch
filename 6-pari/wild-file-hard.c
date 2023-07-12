#include "c-target.h"
#include <stdio.h>

#define pp      plist[pnum]
#define mystart pdata[pnum][0][0]
#define myend   pdata[pnum][0][1]
#define mypoint pdata[pnum][0][2]
#define val     3

static long locorder[7];
static long lookups[7];

/* Let's START *************************************************************/
int w_file_init(long pnum, long which, long what) {
  long j;
  long lessthan, ans;

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
    mypoint = mystart;
    /* Delete all of these */
    pdata[pnum][0][0] = 1;
    pdata[pnum][0][1] = 1;
    pdata[pnum][1][0] = 1; /* modulus */
    pdata[pnum][1][1] = 0; /* starting value */
    pdata[pnum][2][0] = 1;
    pdata[pnum][2][1] = 0;
    /* End future delete */

    for(j=0;j<=DEG;j++) {
      pdata[pnum][j][0] = (long) intpow((long long) pp, caselist[mystart][0]); /* all moduli are this value */
      pflags[pnum][j] = 0;
    }
    locorder[0] = 1; locorder[1] = 2; locorder[2] = 6;
    locorder[3] = 3; locorder[4] = 4; locorder[5] = 5;
    lookups[1] = 0; lookups[2] = 1; lookups[3] = 3;
    lookups[4] = 4; lookups[5] = 5; lookups[6] = 2;    
    break;
  case 1:
  case 2:  /* c */
  case 6:
  case 3:
  case 5:
  case 4: /* Init loop is now on 4 */
    /*
      record incoming entry
      loop through list including current until a match is found
      return 0 if no matches, if there is a match, hold value so next can tell when we are
      on a new value?
     */
    pdata[pnum][which][val] = what;
    lessthan = -1;
    while(lessthan<0 && mypoint < myend) {
      ans = 0;
      for(j=0; (!ans) && j<lookups[which]; j++) {
        if(caselist[mypoint][locorder[j]] > pdata[pnum][locorder[j]][val])
          ans = 1; /* gone too far */
        else
          if(caselist[mypoint][locorder[j]] < pdata[pnum][locorder[j]][val])
            ans = 2; /* skip over this one */
        }
      switch(ans) {
      case 1:
        lessthan = 1; /* stop because we are too far */
        break;
      case 2:
        mypoint++;
        break;
      case 0: /* winner! */
        lessthan = 0;
        break;
      }
    }
    if(mypoint>=myend)
      return -1;
    
    break;
  default:
    break;
  }
  return 1;
}

/* increment counter number pnum modulo p
     returns 1 if ok */
long w_file_next(long pnum, long which) {
  long j, lessthan, ans, oldvalue;

  if(mypoint>=myend)
    return 0;
  /* Are we sitting on an old value? */
  lessthan = -1;
  ans = 0;
  for(j=0; (!ans) && j<lookups[which]; j++) {
    if(caselist[mypoint][locorder[j]] > pdata[pnum][locorder[j]][val])
      return -1;
      ans = 1; /* gone too far */
    else
      if(caselist[mypoint][locorder[j]] < pdata[pnum][locorder[j]][val]) {
        printf("Error in wild module next call");
        return -1;
      }
  } /* If we exit here, we are currently on what appears to be usable */
  oldvalue = case;
  
  
}


/*

We could have just looked at one thing
Inner loops could have moved the pointer forward so we are already on the next good one
They could have also moved us beyond what we can use right now

*/
