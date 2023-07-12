#include "c-target.h"
#include <stdio.h>

#define pp plist[pnum]
#define bb pdata[pnum][8][0]
#define aa pdata[pnum][8][1]
#define cc pdata[pnum][8][2]
#define dd pdata[pnum][8][3]
#define ee pdata[pnum][8][4]
#define ff pdata[pnum][8][5]

/* Assume p>3
   i2 = (p+1)/2;
   i3 = ((p%3) == 1) ?  (2*p+1)/3 : (p+1)/3;
*/

/* Let's START *************************************************************/
int t211111init(long pnum, long which, long what) {
  long long j;

  /*
    We have to give the modulus and the value
     Input of which = 0 says to return for a_1
  */
  switch(which) {
  case 0:
    pdata[pnum][0][0] = 1;
    pdata[pnum][0][1] = 6;
    pdata[pnum][1][0] = 1; /* modulus */
    pdata[pnum][1][1] = 0; /* starting value */
    pdata[pnum][2][0] = 1;
    pdata[pnum][2][1] = 0;
    pdata[pnum][DEG][0] = 1;
    pdata[pnum][DEG][1] = 0;
    pdata[pnum][3][0] = 1;
    pdata[pnum][3][1] = 0;
    pdata[pnum][4][0] = 1;
    pdata[pnum][4][1] = 0;
    pdata[pnum][5][0] = pp;
    pdata[pnum][6][0] = pp;

    pflags[pnum][1] = CALL_NONE;
    for(j=2; j<=7; j++)
      pflags[pnum][j] = CALL_INIT;
    pflags[pnum][5] = CALL_BOTH;
    pflags[pnum][6] = CALL_NONE;
    break;
    /*

if a == 0, loop on d (rest are 0)
if a != 0, any constant term is ok, but this determines d
    */
  case 2:  /* c */
    bb = what; /* for future reference */
    break;
  case DEG:
    ff = what;
    break;
  case 3:
    ee = what; /* for future reference */
    break;
  case 4:
    cc = what;
    break;
  case 5:
    dd = what;
    if(ee==0) {
      aa=0;
      pdata[pnum][5][0] = 1; /* Sets up the a=0 part */
      pdata[pnum][5][1] = 0;
      pdata[pnum][6][1] = 0;
    } else {
      aa=1;
      pdata[pnum][5][0] = pp; /* Sets up the a=1 part */
      pdata[pnum][5][1] = truemod(-5*bb + -3*cc + 2*dd + ee + 4*ff + 6, pp);
      pdata[pnum][6][1] = truemod(-4*bb + -2*cc + dd + 2*ee + 3*ff + 5, pp);
    }
    break;
  default:
    /*    pdata[pnum][which][0] = pp; */
    fprintf(stderr,"Error in t211111 init.\n");
    break;
  }
  return 1;
}

/* increment counter number pnum modulo p
     returns 1 if ok */
long t211111next(long pnum, long which) {
  long long a2, a3, ainv, ainv2, a4,a5;

  if(which != 5) {
    fprintf(stderr,"Error in t211111 next.\n");
    return 0;
  }

  if(aa==0) {
    pdata[pnum][5][0] = pp;
  }

  aa++;
  if(aa==pp)
    return 0;

  
  a2 = (aa*aa) % pp;
  a3 = (a2*aa) % pp;
  a4 = (a3*aa) % pp;
  a5 = (a4*aa) % pp;
  ainv = minv(aa, pp);
  ainv2 = (ainv*ainv) % pp;
    
  pdata[pnum][5][1] = truemod(6*a5 - 5*bb*a4 + 4*ff*a3 - 3*cc*a2 + 2*dd*aa
                              + ee*ainv2, pp);
  pdata[pnum][6][1] = truemod(5*a5*aa - 4*bb*a5 + 3*ff*a4 - 2*cc*a3 + dd*a2
                              + 2*ee*ainv, pp);
  
  return 1;
}


