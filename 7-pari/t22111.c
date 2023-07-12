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
int t22111init(long pnum, long which, long what) {
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
    cc = what; /* for future reference */
    break;
  case DEG:
    dd = what;
    break;
  case 3:
    ee = what; /* for future reference */
    if(ee==0) {
      ff=0;
      aa=0;
      bb=1;
      pdata[pnum][3][1] = 0;
      pdata[pnum][4][1] = 0;
      pdata[pnum][5][1] = 0;
      pdata[pnum][6][1] = 0;
    } else {
      aa=0;
      bb=1;
      ff=pp;
      pdata[pnum][3][1] = truemod(2*cc + ee, pp);
      pdata[pnum][4][1] = truemod(-3 + 2 *dd, pp);
      pdata[pnum][5][1] = truemod(cc + 2 *ee, pp);
      pdata[pnum][6][1] = truemod(-2 + dd, pp);
    }
    break;
  default:
    /*    pdata[pnum][which][0] = pp; */
    fprintf(stderr,"Error in t22111 init.\n");
    break;
  }
  return 1;
}

/* increment counter number pnum modulo p
     returns 1 if ok */
long t22111next(long pnum, long which) {
  long long a2, a3, binv, binv2,a4,a5, b2;

  if(which != 3) {
    fprintf(stderr,"Error in t22111 next.\n");
    return 0;
  }

  aa++;
  if((ee==0) && (ff<pp)) {
    if(aa==pp) {
      aa=0;
      ff++;
    }
    if(ff<pp) {
      
      a2 = (aa*aa) % pp;
      a3 = (a2*aa) % pp;
      a4 = (a3*aa) % pp;
      a5 = (a4*aa) % pp;
      
      pdata[pnum][3][1] = ff;
      pdata[pnum][4][1] = truemod(-5*a4 + 4*cc*a3 - 3*dd*a2 + 2*ff*aa, pp);
      pdata[pnum][5][1] = truemod(-4*a5 + 3*cc*a4 - 2*dd*a3 + ff*a2, pp);
      return 1;
    }
  }

  if(aa==pp) {
    (bb)++;
    if(bb == pp) {
      return 0; /* We are done with the double loop */
    } /* Really an else now */
    aa = 0;
  }
  
  a2 = (aa*aa) % pp;
  a3 = (a2*aa) % pp;
  a4 = (a3*aa) % pp;
  a5 = (a4*aa) % pp;
  b2 = (bb*bb) % pp;
  
  binv = minv(bb, pp);
  binv2 = (binv*binv) % pp;

  pdata[pnum][3][1] = truemod(4*a3 - 3*cc*a2 + (-6*bb + 2*dd)*aa +
                              2*cc*bb + ee*binv2, pp);
  pdata[pnum][4][1] = truemod(3*a4 - 2*cc*a3 + dd*a2 +
                              ((-2*cc*bb + 2*ee*binv2)%pp)*aa +
                              -3*b2 + 2*dd*bb, pp);
  pdata[pnum][5][1] = truemod(6*bb*a3 + ((-4*cc*bb + ee*binv2)%pp)*a2 +
                              ((-6*b2 + 2*dd*bb)%pp)*aa +
                              cc*b2 +  2*ee*binv, pp);
  pdata[pnum][6][1] = truemod(3*b2*a2 + ((-2*cc*b2 + 2*ee*binv)%pp)*aa +
    -2*b2*bb + dd*b2, pp);
  
  return 1;
}


