#include "c-target.h"

#define pp plist[pnum]
#define aa pdata[pnum][8][0]
#define i7 pdata[pnum][8][1]

/* Assume p >5
   i2 = (p+1)/2;
   i3 = ((p%3) == 1) ?  (2*p+1)/3 : (p+1)/3;
*/

/* basic functions:

   init - initialize values given number of known coefs, starting vector of these coefs, p
   next - give the next value after an init has been called, returns 0 if there are no more
*/


/* Initialize a variable counter
    produce next congruence for that counter
*/

int t7init(long pnum, long which, long what) {
  long j;
  long long aaa;

  /*
    We have to give the modulus and the value
     Input of which = 0 says to return for a_1
  */
  switch(which) {
  case 0:
    pdata[pnum][which][0] = 6;
    pdata[pnum][which][1] = 1;
    i7 = (long) minv((long long) 7, (long long) pp);
    
    for(j=1; j<=7; j++)
      pflags[pnum][j] = CALL_NONE;
    pflags[pnum][2] = CALL_INIT;
    pdata[pnum][1][0] = 1; /* modulus */
    pdata[pnum][1][1] = 0; /* starting value */
    for(j=2; j<=DEG; j++)
      pdata[pnum][j][0] = pp;
    break;
  case 2: /* Take the value of a and use it to initialize the rest */
    aa = (what*i7) % pp;  /* Store this for future reference */
    
    aaa = (aa*aa) % pp; /* a^2 */
    pdata[pnum][2][1] = (21*aaa) % pp;
    aaa = (aaa*aa) % pp; /* a^3 */
    pdata[pnum][3][1] = (35*aaa) % pp;
    aaa = (aaa*aa) % pp; /* Now a^4 */
    pdata[pnum][4][1] = (35*aaa) % pp;
    aaa = (aaa*aa) % pp; /* Now a^5 */
    pdata[pnum][5][1] = (21*aaa) % pp;
    aaa = (aaa*aa) % pp; /* Now a^6 */
    pdata[pnum][6][1] = (7*aaa) % pp;
    pdata[pnum][7][1] = (aa*aaa) % pp;
    
    break;
  default:
    break;
  }
  return 1;
}

/* increment counter number pnum modulo p
     returns 1 if ok */
long t7next(long pnum, long which) { /* This should never be called */
  return 0;
}
