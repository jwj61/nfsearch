#include "c-target.h"

#define pp plist[pnum]
#define aa pdata[pnum][8][0]
#define bb pdata[pnum][8][1]
#define i2 pdata[pnum][8][2]
#define i4 pdata[pnum][8][3]

/* basic functions:

   init - initialize values given number of known coefs, starting vector of these coefs, p
   next - give the next value after an init has been called, returns 0 if there are no more
*/


/* Initialize a variable counter
    produce next congruence for that counter
*/

int t52init(long pnum, long which, long what) {
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
    i2 = (pp+1)/2;
    i4 = (i2*i2) % pp;
    break;
  case 2: /* Capture b, loop through values of a */
    bb = what % pp;  /* Store this for future reference */
    aa = 0;
    j = (bb*i2) % pp;
    pdata[pnum][2][1] = truemod(j*j,pp); /* Consequences of a=0 */
    pdata[pnum][3][1] = 0;
    pdata[pnum][4][1] = 0;
    pdata[pnum][5][1] = 0;
    pdata[pnum][6][1] = 0;
    pdata[pnum][7][1] = 0;
    break;
  default:
    fprintf(stderr,"Error in t52 init.\n");
    break;
  }
  return 1;
}

/* increment counter number pnum modulo p
     returns 1 if ok */
long t52next(long pnum, long which) { /* This should never be called */
  long long a2, ab, a4, a5, a3, b2, a6, a7, ib2;
  long long ib2square;
  
  if(which != 2) {
    fprintf(stderr,"Error in t52 next.\n");
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
  ib2 = (bb*i2) % pp;
  ib2square = (ib2*ib2) % pp;
  
  pdata[pnum][2][1] = truemod(-35*i4*a2 + 5*i2*ab + i4*b2, pp);
  pdata[pnum][3][1] = truemod(-35*i4*a3 - 5*ib2*a2 + 5*ib2square*aa, pp);
  pdata[pnum][4][1] = truemod(35*i2*a4 - 15*bb*a3 + 10*ib2square*a2, pp);
  pdata[pnum][5][1] = truemod(77*i2*a5 - 20*bb*a4 + 10*ib2square*a3, pp);
  pdata[pnum][6][1] = truemod(105*i4*a6 - 23*ib2*a5 + 5*ib2square*a4, pp);
  pdata[pnum][7][1] = truemod(25*i4*a7 - 5*ib2*a6 +ib2square*a5, pp);
  return 1;
}
