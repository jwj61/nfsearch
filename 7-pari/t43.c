#include "c-target.h"

#define pp plist[pnum]
#define aa pdata[pnum][8][0]
#define bb pdata[pnum][8][1]
#define i3 pdata[pnum][8][2]
#define i9 pdata[pnum][8][3]
#define i27 pdata[pnum][8][4]

/* basic functions:

   init - initialize values given number of known coefs, starting vector of these coefs, p
   next - give the next value after an init has been called, returns 0 if there are no more
*/


/* Initialize a variable counter
    produce next congruence for that counter
*/

int t43init(long pnum, long which, long what) {
  long j;

  /*
    We have to give the modulus and the value
     Input of which = 0 says to return for a_1
  */
  switch(which) {
  case 0:
    pdata[pnum][which][0] = 5;
    pdata[pnum][which][1] = 2;
    for(j=1; j<=DEG; j++)
      pflags[pnum][j] = CALL_NONE;
    pflags[pnum][2] = CALL_BOTH;
    pdata[pnum][1][0] = 1; /* modulus */
    pdata[pnum][1][1] = 0; /* starting value */
    for(j=2; j<=DEG; j++)
      pdata[pnum][j][0] = pp;
    i3 = minv(3,pp);
    i9 = (i3*i3) % pp;
    i27 = (i3*i9) % pp;
    break;
  case 2: /* Capture b, loop through values of a */
    bb = what;  /* Store this for future reference */
    aa = 0;
    j = (bb*bb) % pp;
    pdata[pnum][2][1] = truemod(i3*j,pp); /* Consequences of a=0 */
    j = (j*bb) % pp;
    pdata[pnum][3][1] = truemod(i27*j,pp);
    pdata[pnum][4][1] = 0;
    pdata[pnum][5][1] = 0;
    pdata[pnum][6][1] = 0;
    pdata[pnum][7][1] = 0;
    break;
  default:
    fprintf(stderr,"Error in t43 init.\n");
    break;
  }
  return 1;
}

/* increment counter number pnum modulo p
     returns 1 if ok */
long t43next(long pnum, long which) { /* This should never be called */
  long long a2, ab, a4, a5, a3, b2, a6, a7, ib3, ai3, b3;
  long long ib3square, ai3square, tmp1, tmp2;
  
  if(which != 2) {
    fprintf(stderr,"Error in t43 next.\n");
    return 0;
  }
  aa++;
  if(aa==pp)
    return 0;
  a2 = (aa*aa) % pp; /* a^2 */
  a3 = (aa*a2) % pp;
  b2 = (bb*bb) % pp;
  b3 = (bb*b2) % pp;
  ab = (bb*aa) % pp; /* ab */
  a4 = (a2*a2) % pp; /* a^4 */
  a5 = (a4*aa) % pp; /* a^5 */
  a6 = (a5*aa) % pp; /* a^5 */
  a7 = (a6*aa) % pp; /* a^5 */
  ib3 = (bb*i3) % pp;
  ai3 = (aa*i3) % pp;
  ai3square = (ai3*ai3) % pp;
  ib3square = (ib3*ib3) % pp;
  
  pdata[pnum][2][1] = truemod(-14*i3*a2 + 4*i3*ab + 3*ib3*ib3, pp);
  pdata[pnum][3][1] = truemod(-28*i27*a3 - 26*bb*ai3square +
                              8*ib3square*aa+ ib3square*ib3, pp);
  tmp1 = (b2*aa) % pp;
  tmp1 = (tmp1*i27) % pp;
  pdata[pnum][4][1] = truemod(203*i27*a4 - 44*ab*ai3square
                              + 2*b2*ai3square + 4*bb* tmp1, pp);
  tmp1 = (4*i3) % pp;
  tmp1 = (tmp1*b2) % pp;
  pdata[pnum][5][1] = truemod(28*i9*a5 + bb*a4 - tmp1*a3
                              + 2*b3*ai3square, pp);
  tmp1 = (40*i9) % pp;
  tmp1 = (tmp1 * bb) % pp;
  tmp2 = (ib3square * ib3) % pp;
  pdata[pnum][6][1] = truemod(-112*i27*a6 + tmp1*a5
                              - 13*ib3square*a4 + 4*tmp2*a3, pp);

  tmp1 = (bb*i9) % pp;
  tmp2 = (ib3square*ib3) % pp; /* b^3 / 9 */
  pdata[pnum][7][1] = truemod(-64*i27*a7 + 16*tmp1*a6
                              - 4*ib3square*a5 +tmp2*a4, pp);
  return 1;
}
