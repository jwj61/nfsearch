#include "c-target.h"

#define pp plist[pnum]
#define aa pdata[pnum][7][0]
#define bb pdata[pnum][7][1]
#define i2 pdata[pnum][7][2]

/* basic functions:

   init - initialize values given number of known coefs, starting vector of these coefs, p
   next - give the next value after an init has been called, returns 0 if there are no more
*/


/* Initialize a variable counter
    produce next congruence for that counter
*/

int t22init(long pnum, long which, long what) {
  long j;
  long long a2, a4;

  /*
    We have to give the modulus and the value
     Input of which = 0 says to return for a_1
  */
  switch(which) {
  case 0:
    pdata[pnum][which][0] = 2;
    pdata[pnum][which][1] = 2;
    for(j=1; j<=DEG; j++)
      pflags[pnum][j] = CALL_NONE;
    pflags[pnum][2] = CALL_INIT;
    pflags[pnum][4] = CALL_INIT;
    pdata[pnum][1][0] = 1; /* modulus */
    pdata[pnum][1][1] = 0; /* starting value */
    pdata[pnum][2][0] = 1; /* modulus */
    pdata[pnum][2][1] = 0; /* starting value */
    for(j=3; j<=DEG; j++)
      pdata[pnum][j][0] = pp;
    i2 = (pp+1)/2;
    break;
  case 2: /* Capture b, loop through values of a */
    aa = (i2*what) % pp;  /* Store this for future reference */
    break;
  case 4:    
    bb = (i2*what) % pp;
    a2 = (aa*aa) % pp;
    a4 = (a2*i2) % pp; /* a^2 /2 */
    a4 = (a4*a4) % pp; /* Now square it */
    pdata[pnum][3][1] = truemod(-aa*a2 + 2*bb*aa, pp);
    pdata[pnum][4][1] = truemod(a4 - bb*a2 + bb*bb, pp);
    break;
  default:
    fprintf(stderr,"Error in t22 init.\n");
    break;
  }
  return 1;
}

/* increment counter number pnum modulo p
     returns 1 if ok */
long t22next(long pnum, long which) { 
  
  fprintf(stderr,"Error in t22 next, which = %ld\n", which);
  return 0;
}

/*

 1^2 1^2
 2*a, 2*b, -a^3 + 2*b*a, 1/4*a^4 - b*a^2 + b^2

 1^2 1 1
 b, c,    4*a^3 - 3*b*a^2 + 2*c*a,    3*a^4 - 2*b*a^3 + c*a^2

 
*/
