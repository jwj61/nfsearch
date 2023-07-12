#include "c-target.h"

#define pp plist[pnum]
#define aa pdata[pnum][7][0]
#define bb pdata[pnum][7][1]

/* basic functions:

   init - initialize values given number of known coefs, starting vector of these coefs, p
   next - give the next value after an init has been called, returns 0 if there are no more
*/


/* Initialize a variable counter
    produce next congruence for that counter
*/

int t31init(long pnum, long which, long what) {
  long j;

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
    pflags[pnum][2] = CALL_BOTH;
    pdata[pnum][1][0] = 1; /* modulus */
    pdata[pnum][1][1] = 0; /* starting value */
    for(j=2; j<=DEG; j++)
      pdata[pnum][j][0] = pp;
    break;
  case 2: /* Capture b, loop through values of a */
    bb = what;  /* Store this for future reference */
    aa = 0;
    pdata[pnum][2][1] = 0; /* Consequences of a=0 */
    pdata[pnum][3][1] = 0;
    pdata[pnum][4][1] = 0;
    break;
  default:
    fprintf(stderr,"Error in t31 init.\n");
    break;
  }
  return 1;
}

/* increment counter number pnum modulo p
     returns 1 if ok */
long t31next(long pnum, long which) { 
  long long a2, ab;
  
  if(which != 2) {
    fprintf(stderr,"Error in t31 next, which = %ld\n", which);
    return 0;
  }
  aa++;
  if(aa==pp)
    return 0;
  a2 = (aa*aa) % pp;
  ab = (bb*aa) % pp;
  
  pdata[pnum][2][1] = truemod(-6*a2 + 3*ab, pp);
  pdata[pnum][3][1] = truemod(-8*a2*aa + 3*bb*a2, pp);
  pdata[pnum][4][1] = truemod(-3*a2*a2 + ab*a2, pp);
  return 1;
}

/*
 1^3 1
 b,   -6*a^2 + 3*b*a,   -8*a^3 + 3*b*a^2,   -3*a^4 + b*a^3

 1^2 1^2
 2*a, 2*b, -a^3 + 2*b*a, 1/4*a^4 - b*a^2 + b^2

 1^2 1 1
 b, c,    4*a^3 - 3*b*a^2 + 2*c*a,    3*a^4 - 2*b*a^3 + c*a^2

 
*/
