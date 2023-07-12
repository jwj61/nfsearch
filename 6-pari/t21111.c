#include "c-target.h"
#include <stdio.h>

#define pp plist[pnum]
#define bb pdata[pnum][7][0]
#define aa pdata[pnum][7][1]
#define cc pdata[pnum][7][2]
#define dd pdata[pnum][7][3]
#define ee pdata[pnum][7][4]

/* Let's START *************************************************************/
int t21111init(long pnum, long which, long what) {
  long long j;

  /*
    We have to give the modulus and the value
     Input of which = 0 says to return for a_1
  */
  switch(which) {
  case 0:
    pdata[pnum][0][0] = 1;
    pdata[pnum][0][1] = 1;
    pdata[pnum][1][0] = 1; /* modulus */
    pdata[pnum][1][1] = 0; /* starting value */
    pdata[pnum][2][0] = 1;
    pdata[pnum][2][1] = 0;
    pdata[pnum][6][0] = 1;
    pdata[pnum][6][1] = 0;
    pdata[pnum][3][0] = 1;
    pdata[pnum][3][1] = 0;
    pdata[pnum][4][0] = pp;
    pdata[pnum][5][0] = pp;
    
    
    for(j=1; j<=5; j++)
      pflags[pnum][j] = CALL_NONE;
    pflags[pnum][2] = CALL_INIT;
    pflags[pnum][6] = CALL_INIT; 
    pflags[pnum][3] = CALL_INIT;
    pflags[pnum][4] = CALL_BOTH;
    break;
  case 2:  /* c */
    bb = what; /* for future reference */
    break;
  case 6:
    ee = what;
    break;
  case 3:
    dd = what; /* for future reference */
    break;
  case 4: /* Init loop is now on 4 */
    cc = what;
    if(dd==0) {
      aa=0;
      pdata[pnum][4][0] = 1; /* Sets up the a=0 part */
      pdata[pnum][4][1] = 0;
      pdata[pnum][5][1] = 0;
    } else {
      aa=1;
      pdata[pnum][4][0] = pp;
      pdata[pnum][4][1] = truemod(4*bb -3*ee + 2*cc + dd - 5, pp); /* 4*b + -3*e + 2*c + d - 5 */
      pdata[pnum][5][1] = truemod(3*bb -2*ee + cc + 2*dd - 4, pp); /* 3*b + -2*e + c + 2*d - 4  */
      
    }
    break;
  default:
    /*    pdata[pnum][which][0] = pp; */
    fprintf(stderr,"Error in t21111 init.\n");
    break;
  }
  return 1;
}

/* increment counter number pnum modulo p
     returns 1 if ok */
long t21111next(long pnum, long which) {
  long long a2, a3, ainv, ainv2, a4;

  if(which != 4) {
    fprintf(stderr,"Error in t21111 next.\n");
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
  a4= (a2*a2) % pp;
  ainv = minv(aa, pp);
  ainv2 = (ainv*ainv) % pp;
  
  pdata[pnum][4][1] = truemod(4*bb*a3-5*a2*a2-3*ee*a2+2*cc*aa+dd*ainv2, pp);
  pdata[pnum][5][1] = truemod(3*bb*a4 -4*a3*a2-2*ee*a3+cc*a2
                              +2*dd*ainv, pp);
  return 1;
}


/* 

21111

From a,b,c,d,e:
 [b, e, c, a, 0, 0]

 or
 
 [b, e, c,
   (-5*a^6 + 4*b*a^5 - 3*e*a^4 + 2*c*a^3 + d)/a^2,
   (-4*a^6 + 3*b*a^5 - 2*e*a^4 + c*a^3 + 2*d)/a,
   d]
*/
