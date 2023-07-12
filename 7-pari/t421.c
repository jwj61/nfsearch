#include "c-target.h"

#define pp plist[pnum]
#define bb pdata[pnum][8][0]
#define aa pdata[pnum][8][1]
#define cc pdata[pnum][8][2]

/* Assume p>3
   i2 = (p+1)/2;
   i3 = ((p%3) == 1) ?  (2*p+1)/3 : (p+1)/3;
*/

/* Let's START *************************************************************/
int t421init(long pnum, long which, long what) {
  long jj;

  /*
    We have to give the modulus and the value
     Input of which = 0 says to return for a_1
  */
  
  switch(which) {
  case 0:
    pdata[pnum][0][0] = 4;
    pdata[pnum][0][1] = 3;
    pdata[pnum][1][0] = 1; /* modulus */
    pdata[pnum][1][1] = 0; /* starting value */

    pflags[pnum][1] = CALL_NONE;
    pflags[pnum][2] = CALL_BOTH;
    for(jj=3;jj<=DEG; jj++)
      pflags[pnum][jj] = CALL_NONE;
    for(jj=2;jj<=DEG; jj++)
      pdata[pnum][jj][0] = pp;
    
    break;
    
  case 2:  /* c */
    cc = what; /* for future reference */
    aa = 0; /* ready for double loop */
    bb = 0;
    for(jj=2; jj<=DEG; jj++) {
      pdata[pnum][jj][1] = 0;
    }
    
    break;
  default:
    /*    pdata[pnum][which][0] = pp; */
    fprintf(stderr,"Error in t421 init.\n");
    break;
  }
  return 1;
}

/* increment counter number pnum modulo p
     returns 1 if ok */
long t421next(long pnum, long which) {
  long long a2, a3, b2, b3, a4, a5, bc, b2c;
  
  /* Should only be called with which ==2 */
  if(which != 2) {
    fprintf(stderr,"Error in t421 next.\n");
    return 0;
  }
  
  (aa)++;
  if(aa == pp) {
    (bb)++;
    if(bb == pp) {
      return 0; /* We are done with the double loop */
    } /* Really an else now */
    aa = 0;
  }
  /* We now have a new aa and a new bb, time to set everyone: */
  a2 = (aa*aa) % pp;  a3 = (a2 * aa) % pp; a4 = (a3*aa) % pp;
  a5 = (a4*aa) % pp;
  b2 = (bb*bb) % pp;  b3 = (b2 * bb) % pp;
  bc = (bb*cc) % pp;
  b2c = (b2*cc) % pp;
  
  pdata[pnum][2][1] = truemod(-10*a2 + (-8*bb + 4*cc)*aa +
                              (-3*b2 + 2*cc*bb),pp);
  pdata[pnum][3][1] = truemod(-20*a3 + (-32*bb +6*cc)*a2 + (-16*b2 + 8*bc)*aa
                              + (-2*b3 + cc*b2),pp);
  pdata[pnum][4][1] = truemod(-15*a4 + (-48*bb +4*cc)*a3 +
                              (-34*b2 + 12*bc)*a2 +
                              (-8*b3 + 4*b2c)*aa,pp);
  pdata[pnum][5][1] = truemod(-4*a5 +(-32*bb + cc)*a4 +
                              (-36*b2 + 8*bc)*a3 +
                              (-12*b3 + 6*b2c)*a2,pp);
  pdata[pnum][6][1] = truemod(-8*bb*a5 + (-19*b2 + 2*bc)*a4
                              + (-8*b3 + 4*b2c)*a3,pp);
  pdata[pnum][7][1] = truemod(-4*b2*a5 + (-2*b3 + b2c)*a4,pp);
  
  return 1;
}

