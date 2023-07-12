#include "c-target.h"

#define pp plist[pnum]
#define aa pdata[pnum][7][0]
#define i2 pdata[pnum][7][1]
#define i4 pdata[pnum][7][2]

/* basic functions:

   init - initialize values given number of known coefs, starting vector of these coefs, p
   next - give the next value after an init has been called, returns 0 if there are no more
*/


/* Initialize a variable counter
    produce next congruence for that counter
*/

int t4init(long pnum, long which, long what) {
  long j;
  long long aaa;

  /*
    We have to give the modulus and the value
     Input of which = 0 says to return for a_1
  */
  switch(which) {
  case 0:
    pdata[pnum][which][0] = 3;
    pdata[pnum][which][1] = 1;
    for(j=1; j<=4; j++)
      pflags[pnum][j] = CALL_NONE;
    pflags[pnum][2] = CALL_INIT;
    pdata[pnum][1][0] = 1; /* modulus */
    pdata[pnum][1][1] = 0; /* starting value */
    for(j=2; j<=DEG; j++)
      pdata[pnum][j][0] = pp;
    i2 = (pp+1)/2;
    i4 = (i2*i2) % pp;
    break;
  case 2: /* Take the value of a and use it to initialize the rest */
    aa = (i4*what) % pp;  /* Store this for future reference */
    aaa = (aa*aa) % pp; /* a^2 */
    pdata[pnum][2][1] = (6*aaa) % pp;
    aaa = (aa*aaa) % pp; /* Now as a^3 */
    pdata[pnum][3][1] = (4*aaa) % pp;
    pdata[pnum][4][1] = (aa*aaa) % pp;
    break;
  default:
    fprintf(stderr,"Error in t4 init.\n");
    break;
  }
  return 1;
}

/* increment counter number pnum modulo p
     returns 1 if ok */
long t4next(long pnum, long which) { /* This should never be called */
  fprintf(stderr,"Error in t4 next, which = %ld\n", which);
  return 0;
}

/*
 1^4
 4*a, 6*a^2, 4*a^3, a^4

 1^3 1
 b,   -6*a^2 + 3*b*a,   -8*a^3 + 3*b*a^2,   -3*a^4 + b*a^3

 1^2 1^2
 2*a, 2*b, -a^3 + 2*b*a, 1/4*a^4 - b*a^2 + b^2

 1^2 1 1
 b, c,    4*a^3 - 3*b*a^2 + 2*c*a,    3*a^4 - 2*b*a^3 + c*a^2

 
*/
