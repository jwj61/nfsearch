#include "c-target.h"

#define pp plist[pnum]
#define aa pdata[pnum][8][0]
#define bb pdata[pnum][8][1]

/* basic functions:

   init - initialize values given number of known coefs, starting vector of these coefs, p
   next - give the next value after an init has been called, returns 0 if there are no more
*/


/* Initialize a variable counter
    produce next congruence for that counter
*/

int t61init(long pnum, long which, long what) {
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
    break;
  case 2: /* Capture b, loop through values of a */
    bb = what % pp;  /* Store this for future reference */
    aa = 0;
    pdata[pnum][2][1] = 0; /* Consequences of a=0 */
    pdata[pnum][3][1] = 0;
    pdata[pnum][4][1] = 0;
    pdata[pnum][5][1] = 0;
    pdata[pnum][6][1] = 0;
    pdata[pnum][7][1] = 0;
    break;
  default:
    fprintf(stderr,"Error in t61 init.\n");
    break;
  }
  return 1;
}

/* increment counter number pnum modulo p
     returns 1 if ok */
long t61next(long pnum, long which) { /* This should never be called */
  long long a2, ab, a4, a5;
  
  if(which != 2) {
    fprintf(stderr,"Error in t61 next.\n");
    return 0;
  }
  aa++;
  if(aa==pp)
    return 0;
  a2 = (aa*aa) % pp; /* a^2 */
  ab = (bb*aa) % pp; /* ab */
  a4 = (a2*a2) % pp; /* a^4 */
  a5 = (a4*aa) % pp; /* a^5 */
  
  pdata[pnum][2][1] = truemod(-21*a2 + 6*ab, pp);
  pdata[pnum][3][1] = truemod(-70*a2*aa + 15*ab*aa, pp);
  pdata[pnum][4][1] = truemod(-105*a4 + 20*ab*a2, pp);
  pdata[pnum][5][1] = truemod(-84*a5 +15*bb*a4, pp);
  pdata[pnum][6][1] = truemod(-35*a4*a2 + 6*ab*a4, pp);
  pdata[pnum][7][1] = truemod(-6*a5*a2 + ab*a5, pp);
  return 1;
}

/*
  1^5
  [5*a, 10*a^2, 10*a^3, 5*a^4, a^5]

  1^4 1
  [b, -10*a^2 + 4*b*a, -20*a^3 + 6*b*a^2, -15*a^4 + 4*b*a^3, -4*a^5 + b*a^4]

  1^3 1 1
  [b, c, 10*a^3 - 6*b*a^2 + 3*c*a, 15*a^4 - 8*b*a^3 + 3*c*a^2, 6*a^5 - 3*b*a^4 + c*a^3]

   1^3 1^2
   [b, -15/4*a^2 + 3/2*b*a + 1/4*b^2, -5/4*a^3 - 3/2*b*a^2 + 3/4*b^2*a, 15/4*a^4 - 7/2*b*a^3 + 3/4*b^2*a^2, 9/4*a^5 - 3/2*b*a^4 + 1/4*b^2*a^3]

   1^2 1^2 1
   [c, b, -5*a^3 + 6*c*a^2 + (-b - 2*c^2)*a + c*b, -15/4*a^4 + 4*c*a^3 + (-1/2*b - c^2)*a^2 + 1/4*b^2, -9/2*a^5 + 33/4*c*a^4 + (-3*b - 5*c^2)*a^3 + (7/2*c*b + c^3)*a^2 + (-1/2*b^2 - c^2*b)*a + 1/4*c*b^2]

   1^2 1 1 1
   if a!= 0,
   [b, c, (4*a^5 - 3*b*a^4 + 2*c*a^3 + d)/a^2, (3*a^5 - 2*b*a^4 + c*a^3 + 2*d)/a, d]
   else
   x^5 + b*x^4 + c*x^3 + d*x^2
   
*/
