#include "c-target.h"
#include <stdio.h>

#define pp plist[pnum]
#define bb pdata[pnum][8][0]
#define aa pdata[pnum][8][1]
#define cc pdata[pnum][8][2]
#define dd pdata[pnum][8][3]

/* Assume p>3
   i2 = (p+1)/2;
   i3 = ((p%3) == 1) ?  (2*p+1)/3 : (p+1)/3;
*/

/* Let's START *************************************************************/
int t4111init(long pnum, long which, long what) {
  long long j;

  /*
    We have to give the modulus and the value
     Input of which = 0 says to return for a_1
  */
  switch(which) {
  case 0:
    pdata[pnum][0][0] = 3;
    pdata[pnum][0][1] = 4;
    pdata[pnum][1][0] = 1; /* modulus */
    pdata[pnum][1][1] = 0; /* starting value */
    pdata[pnum][2][0] = 1;
    pdata[pnum][2][1] = 0;
    pdata[pnum][DEG][0] = 1;
    pdata[pnum][DEG][1] = 0;
    pdata[pnum][3][0] = pp;
    pdata[pnum][4][0] = pp;
    pdata[pnum][5][0] = pp;
    pdata[pnum][6][0] = pp;

    for(j=1; j<=DEG; j++)
      pflags[pnum][j] = CALL_NONE;
    pflags[pnum][2] = CALL_INIT;
    pflags[pnum][DEG] = CALL_INIT; 
    pflags[pnum][3] = CALL_BOTH;
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
    dd = what; /* for future reference */
    if(dd==0) {
      aa=0;
      pdata[pnum][3][0] = 1; /* Sets up the a=0 part */
      pdata[pnum][3][1] = 0;
      pdata[pnum][4][1] = 0;
      pdata[pnum][5][1] = 0;
      pdata[pnum][6][1] = 0;
    } else {
      aa=1;
      pdata[pnum][3][0] = pp; /* Sets up the a=1 part */
      pdata[pnum][3][1] = truemod(-10*bb + 4*cc + dd + 20, pp);
      pdata[pnum][4][1] = truemod(-20*bb + 6*cc + 4*dd + 45, pp);
      pdata[pnum][5][1] = truemod(-15*bb+ 4*cc + 6*dd + 36, pp);
      pdata[pnum][6][1] = truemod(-4*bb + cc + 4*dd + 10, pp);
    }
    break;
  default:
    /*    pdata[pnum][which][0] = pp; */
    fprintf(stderr,"Error in t4111 init.\n");
    break;
  }
  return 1;
}

/* increment counter number pnum modulo p
     returns 1 if ok */
long t4111next(long pnum, long which) {
  long long a2, a3, ainv, ainv2, ainv3, a4,a5, ainv4;

  if(which != 3) {
    fprintf(stderr,"Error in t4111 next.\n");
    return 0;
  }

  if(aa==0) {
    pdata[pnum][3][0] = pp;
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
  ainv4 = (ainv3*ainv) % pp;
    
  pdata[pnum][3][1] = truemod(-10*a2*bb + 4*aa*cc + 20*a3
                              +dd*ainv4, pp);
  pdata[pnum][4][1] = truemod(-20*a3*bb + 6*a2*cc + 45*a4
                              +4*dd*ainv3, pp);
  pdata[pnum][5][1] = truemod(-15*a4*bb + 4*a3*cc + 36*a5 + 6*dd*ainv2, pp);
  pdata[pnum][6][1] = truemod(-4*a5*bb + a4*cc + 10*a4*a2  +4*dd*ainv, pp);
  
  return 1;
}


