#include "c-target.h"
#include <stdio.h>

#define pp plist[pnum]
#define bb pdata[pnum][8][0]
#define aa pdata[pnum][8][1]
#define cc pdata[pnum][8][2]
#define dd pdata[pnum][8][3]
#define ee pdata[pnum][8][4]

/* Assume p>3
   i2 = (p+1)/2;
   i3 = ((p%3) == 1) ?  (2*p+1)/3 : (p+1)/3;
*/

/* Let's START *************************************************************/
int t31111init(long pnum, long which, long what) {
  long long j;

  /*
    We have to give the modulus and the value
     Input of which = 0 says to return for a_1
  */
  switch(which) {
  case 0:
    pdata[pnum][0][0] = 2;
    pdata[pnum][0][1] = 5;
    pdata[pnum][1][0] = 1; /* modulus */
    pdata[pnum][1][1] = 0; /* starting value */
    pdata[pnum][2][0] = 1;
    pdata[pnum][2][1] = 0;
    pdata[pnum][DEG][0] = 1;
    pdata[pnum][DEG][1] = 0;
    pdata[pnum][3][0] = 1;
    pdata[pnum][3][1] = 0;
    pdata[pnum][4][0] = pp;
    pdata[pnum][5][0] = pp;
    pdata[pnum][6][0] = pp;

    for(j=1; j<=DEG; j++)
      pflags[pnum][j] = CALL_NONE;
    pflags[pnum][2] = CALL_INIT;
    pflags[pnum][3] = CALL_INIT;
    pflags[pnum][DEG] = CALL_INIT;
    pflags[pnum][4] = CALL_BOTH;
    break;
    /*

if a == 0, loop on d (rest are 0)
if a != 0, any constant term is ok, but this determines d
    */
  case 2:  /* c */
    bb = what; /* for future reference */
    break;
  case DEG:
    cc = what;
    break;
  case 3:
    ee = what; /* for future reference */
    break;
  case 4:
    dd = what;
    if(ee==0) {
      aa=0;
      pdata[pnum][4][0] = 1; /* Sets up the a=0 part */
      pdata[pnum][4][1] = 0;
      pdata[pnum][5][1] = 0;
      pdata[pnum][6][1] = 0;
    } else {
      aa=1;
      pdata[pnum][4][0] = pp; /* Sets up the a=1 part */
      pdata[pnum][4][1] = truemod(10*bb -6*cc + 3*dd + ee - 15, pp);
      pdata[pnum][5][1] = truemod(15*bb -8*cc + 3*dd + 3*ee - 24, pp);
      pdata[pnum][6][1] = truemod(6*bb -3*cc +dd + 3*ee - 10, pp);
    }
    break;
  default:
    /*    pdata[pnum][which][0] = pp; */
    fprintf(stderr,"Error in t31111 init.\n");
    break;
  }
  return 1;
}

/* increment counter number pnum modulo p
     returns 1 if ok */
long t31111next(long pnum, long which) {
  long long a2, a3, ainv, ainv2, ainv3, a4,a5;

  if(which != 4) {
    fprintf(stderr,"Error in t31111 next.\n");
    return 0;
  }

  if(aa==0) {
    pdata[pnum][4][0] = pp;
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
  ainv3 = (ainv2*ainv) % pp;
    
  pdata[pnum][4][1] = truemod(-15*a4 + 10*bb*a3 - 6*cc*a2 + 3*dd*aa
                              + ee*ainv3, pp);
  pdata[pnum][5][1] = truemod(-24*a5 + 15*bb*a4 - 8*cc*a3 + 3*dd*a2
                              + 3*ee*ainv2, pp);
  pdata[pnum][6][1] = truemod(-10*a5*aa + 6*bb*a5 - 3*cc*a4 + dd*a3
                              +3* ee*ainv, pp);
  
  /*
      (-15*a^7 + 10*b*a^6 - 6*c*a^5 + 3*d*a^4 + e)/a^3,
  (-24*a^7 + 15*b*a^6 - 8*c*a^5 + 3*d*a^4 + 3*e)/a^2,
  (-10*a^7 + 6*b*a^6 - 3*c*a^5 + d*a^4 + 3*e)/a,
  */
  
  return 1;
}


