#include "c-target.h"

#define pp plist[pnum]
#define aa pdata[pnum][7][0]
#define i3 pdata[pnum][7][1]

/* basic functions:

   init - initialize values given number of known coefs, starting vector of these coefs, p
   next - give the next value after an init has been called, returns 0 if there are no more
*/


/* Initialize a variable counter
    produce next congruence for that counter
*/

int t3init(long pnum, long which, long what) {
  long j;
  long long aaa;

  /*
    We have to give the modulus and the value
     Input of which = 0 says to return for a_1
  */
  switch(which) {
  case 0:
    pdata[pnum][which][0] = 2;
    pdata[pnum][which][1] = 1;
    for(j=1; j<=3; j++)
      pflags[pnum][j] = CALL_NONE;
    pflags[pnum][2] = CALL_INIT;
    pdata[pnum][1][0] = 1; /* modulus */
    pdata[pnum][1][1] = 0; /* starting value */
    for(j=2; j<=DEG; j++)
      pdata[pnum][j][0] = pp;
    i3 = (long) minv(3, pp);
    break;
  case 2: /* Take the value of a and use it to initialize the rest */
    aa = (i3*what) % pp;  /* Store this for future reference */
    aaa = (aa*aa) % pp; /* a^2 */
    pdata[pnum][2][1] = (3*aaa) % pp;
    pdata[pnum][3][1] = (aa*aaa) % pp;
    
    break;
  default:
    break;
  }
  return 1;
}

/* increment counter number pnum modulo p
     returns 1 if ok */
long t3next(long pnum, long which) { /* This should never be called */
  return 0;
}

/*
 1^3
 3*a,    3*a^2,    a^3

 1^2 1
 b,    -3*a^2 + 2*b*a,    -2*a^3 + b*a^2
*/
