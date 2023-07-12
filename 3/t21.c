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

int t21init(long pnum, long which, long what) {
  long j;

  /*
    We have to give the modulus and the value
     Input of which = 0 says to return for a_1
  */
  switch(which) {
  case 0:
    pdata[pnum][which][0] = 1;
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
    break;
  default:
    fprintf(stderr,"Error in t21 init.\n");
    break;
  }
  return 1;
}

/* increment counter number pnum modulo p
     returns 1 if ok */
long t21next(long pnum, long which) { 
  long long a2, ab;
  
  if(which != 2) {
    fprintf(stderr,"Error in t21 next, which = %ld\n", which);
    return 0;
  }
  aa++;
  if(aa==pp)
    return 0;
  a2 = (aa*aa) % pp;
  ab = (bb*aa) % pp;
  
  pdata[pnum][2][1] = truemod(-3*a2 + 2*ab, pp);
  pdata[pnum][3][1] = truemod(-2*a2*aa + bb*a2, pp);
  return 1;
}

/*
1^2 1
 b,    -3*a^2 + 2*b*a,    -2*a^3 + b*a^2
*/
