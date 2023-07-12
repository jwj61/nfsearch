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

int t32init(long pnum, long which, long what) {
  long j;
  long long tmp;

  /*
    We have to give the modulus and the value
     Input of which = 0 says to return for a_1
  */
  switch(which) {
  case 0:
    pdata[pnum][which][0] = 3;
    pdata[pnum][which][1] = 2;
    i2 = (pp+1)/2;
    for(j=1; j<=5; j++)
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
    tmp = (bb*i2) % pp;
    pdata[pnum][2][1] = (tmp*tmp) % pp; /* Consequences of a=0 */
    pdata[pnum][3][1] = 0;
    pdata[pnum][4][1] = 0;
    pdata[pnum][5][1] = 0;
    break;
  default:
    fprintf(stderr,"Error in t32 init.\n");
    break;
  }
  return 1;
}

/* increment counter number pnum modulo p
     returns 1 if ok */
long t32next(long pnum, long which) { /* This should never be called */
  long long a2, ab, a4, a3, b2, i4;
  
  if(which != 2) {
    fprintf(stderr,"Error in t32 next.\n");
    return 0;
  }
  aa++;
  if(aa==pp)
    return 0;
  a2 = (aa*aa) % pp;
  a3 = (aa*a2) % pp;
  ab = (bb*aa) % pp;
  a4 = (a2*a2) % pp;
  b2 = (bb*bb) % pp;
  i4 = (i2*i2) % pp;
  
  pdata[pnum][2][1] = truemod(-15*i4*a2 + 3*i2*ab + i4*b2, pp);
  pdata[pnum][3][1] = truemod(-5*i4*a3 - 3*i2*bb*a2 + 3*i4*b2*aa, pp);
  pdata[pnum][4][1] = truemod(15*i4*a4 - 7*i2*ab*a2 + 3*i4*b2*a2, pp);
  pdata[pnum][5][1] = truemod(9*i4*a4*aa - 3*i2*bb*a4 + i4*b2*a3, pp);
  return 1;
}

/*
  1^3 1 1
  [b, c, 10*a^3 - 6*b*a^2 + 3*c*a, 15*a^4 - 8*b*a^3 + 3*c*a^2, 6*a^5 - 3*b*a^4 + c*a^3]

   1^3 1^2
   [b,
   -15/4*a^2 + 3/2*b*a + 1/4*b^2,
   -5/4*a^3 - 3/2*b*a^2 + 3/4*b^2*a,
   15/4*a^4 - 7/2*b*a^3 + 3/4*b^2*a^2,
   9/4*a^5 - 3/2*b*a^4 + 1/4*b^2*a^3]

   1^2 1^2 1
   [c, b, -5*a^3 + 6*c*a^2 + (-b - 2*c^2)*a + c*b, -15/4*a^4 + 4*c*a^3 + (-1/2*b - c^2)*a^2 + 1/4*b^2, -9/2*a^5 + 33/4*c*a^4 + (-3*b - 5*c^2)*a^3 + (7/2*c*b + c^3)*a^2 + (-1/2*b^2 - c^2*b)*a + 1/4*c*b^2]

   1^2 1 1 1
   if a!= 0,
   [b, c, (4*a^5 - 3*b*a^4 + 2*c*a^3 + d)/a^2, (3*a^5 - 2*b*a^4 + c*a^3 + 2*d)/a, d]
   else
   x^5 + b*x^4 + c*x^3 + d*x^2
   
*/
