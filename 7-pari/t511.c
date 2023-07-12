#include "c-target.h"

#define pp plist[pnum]
#define aa pdata[pnum][8][0]
#define bb pdata[pnum][8][1]
#define cc pdata[pnum][8][2]

/* basic functions:

   init - initialize values given number of known coefs, starting vector of these coefs, p
   next - give the next value after an init has been called, returns 0 if there are no more
*/


/* Initialize a variable counter
    produce next congruence for that counter
*/

int t511init(long pnum, long which, long what) {
  long j;

  /*
    We have to give the modulus and the value
     Input of which = 0 says to return for a_1
  */
  switch(which) {
  case 0:
    pdata[pnum][which][0] = 4;
    pdata[pnum][which][1] = 3;
    for(j=1; j<=DEG; j++)
      pflags[pnum][j] = CALL_NONE;
    pflags[pnum][2] = CALL_INIT;
    pflags[pnum][DEG] = CALL_BOTH;
    pdata[pnum][1][0] = 1; /* modulus */
    pdata[pnum][1][1] = 0; /* starting value */
    pdata[pnum][2][0] = 1; /* modulus */
    pdata[pnum][2][1] = 0; /* starting value */
    for(j=3; j<=DEG; j++)
      pdata[pnum][j][0] = pp;
    break;
  case 2:
    bb = what;
    break;
  case DEG: /* Capture b, loop through values of a */
    cc = what;  /* Store this for future reference */
    aa = 0;
    pdata[pnum][3][1] = 0;
    pdata[pnum][4][1] = 0;
    pdata[pnum][5][1] = 0;
    pdata[pnum][6][1] = 0;
    pdata[pnum][7][1] = 0;
    break;
  default:
    fprintf(stderr,"Error in t511 init.\n");
    break;
  }
  return 1;
}

/* increment counter number pnum modulo p
     returns 1 if ok */
long t511next(long pnum, long which) { /* This should never be called */
  long long a2, ab, a4, a5, a3, b2, a6, a7;
  
  if(which != DEG) {
    fprintf(stderr,"Error in t511 next.\n");
    return 0;
  }
  aa++;
  if(aa==pp)
    return 0;
  a2 = (aa*aa) % pp; /* a^2 */
  a3 = (aa*a2) % pp;
  b2 = (bb*bb) % pp;
  ab = (bb*aa) % pp; /* ab */
  a4 = (a2*a2) % pp; /* a^4 */
  a5 = (a4*aa) % pp; /* a^5 */
  a6 = (a5*aa) % pp; /* a^5 */
  a7 = (a6*aa) % pp; /* a^5 */
  
  pdata[pnum][3][1] = truemod(35*a3 - 15*bb*a2 + 5*cc*aa, pp);
  pdata[pnum][4][1] = truemod(105*a4 - 40*bb*a3 + 10*cc*a2, pp);
  pdata[pnum][5][1] = truemod(126*a5 - 45*bb*a4 + 10*cc*a3, pp);
  pdata[pnum][6][1] = truemod(70*a6 - 24*bb*a5 + 5*cc*a4, pp);
  pdata[pnum][7][1] = truemod(15*a7 -5*bb*a6 + cc*a5, pp);
  return 1;
}
