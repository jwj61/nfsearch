#include "c-target.h"

#define pp plist[pnum]
#define aa pdata[pnum][7][0]
#define bb pdata[pnum][7][1]
#define cc pdata[pnum][7][2]
#define i2 pdata[pnum][7][3]
#define i4 pdata[pnum][7][4]

/* basic functions:

   init - initialize values given number of known coefs, starting vector of these coefs, p
   next - give the next value after an init has been called, returns 0 if there are no more
*/


/* Initialize a variable counter
    produce next congruence for that counter
*/

int t221init(long pnum, long which, long what) {
  long j;
  long long b2, tmp;

  /*
    We have to give the modulus and the value
     Input of which = 0 says to return for a_1
  */
  switch(which) {
  case 0:
    pdata[pnum][which][0] = 2; /* disc exponent */
    pdata[pnum][which][1] = 3; /* Number of prime factors */
    i2 = (pp+1)/2;
    i4 = (i2*i2) % pp;
    for(j=1; j<=5; j++)
      pflags[pnum][j] = CALL_NONE;
    pflags[pnum][2] = CALL_INIT;
    pflags[pnum][5] = CALL_BOTH;
    pdata[pnum][1][0] = 1; /* modulus */
    pdata[pnum][1][1] = 0; /* starting value */
    pdata[pnum][2][0] = 1;
    pdata[pnum][2][1] = 0;
    for(j=3; j<=DEG; j++)
      pdata[pnum][j][0] = pp;
    break;
  case 2: /* Capture b, loop through values of a */
    cc = what % pp;  /* Store this for future reference */
    break;
  case 5:
    bb = what % pp;
    aa = 0;
    pdata[pnum][3][1] = (cc*bb) % pp; /* Consequences of a=0 */
    b2 = (bb*bb) % pp;
    pdata[pnum][4][1] = (i4*b2) % pp;
    tmp = (i4*b2) % pp;
    pdata[pnum][5][1] = (tmp*cc) % pp;
    break;
  default:
    fprintf(stderr,"Error in t221 init.\n");
    break;
  }
  return 1;
}

/* increment counter number pnum modulo p
     returns 1 if ok */
long t221next(long pnum, long which) { /* This should never be called */
  long long a2, a4, a3, c2, b2, i2bb, a5, i4a4, coef1, coef2, coef3;
  
  if(which != 5) {
    fprintf(stderr,"Error in t221 next.\n");
    return 0;
  }
  aa++;
  if(aa==pp)
    return 0;
  a2 = (aa*aa) % pp;
  a3 = (aa*a2) % pp;
  a4 = (a2*a2) % pp;
  a5 = (a4*aa) % pp;
  c2 = (cc*cc) % pp;
  b2 = (bb*bb) % pp;
  i2bb = (bb*i2) % pp;
  i4a4 = (i4*a4) % pp;
  coef1 = (7*cc*i2bb + c2*cc) % pp;
  coef2 = (-i2*b2 - c2*bb) % pp;
  coef3 = (i4*cc) % pp;
  
  pdata[pnum][3][1] = truemod(-5*a3 + 6*cc*a2 + (-bb - 2*c2)*aa + cc*bb, pp);
  pdata[pnum][4][1] = truemod(-15*i4*a4 + 4*cc*a3 + (-i2bb - c2)*a2 + i4*b2, pp);
  pdata[pnum][5][1] = truemod(-9*i2*a5 + 33*i4a4*cc + (-3*bb - 5*c2)*a3 + coef1*a2 + coef2*aa + coef3*b2, pp);
  return 1;
}

/*
   1^2 1^2 1
   [c, b,
   -5*a^3 + 6*c*a^2 + (-b - 2*c^2)*a + c*b,
   -15/4*a^4 + 4*c*a^3 + (-1/2*b - c^2)*a^2 + 1/4*b^2,
   -9/2*a^5 + 33/4*c*a^4 + (-3*b - 5*c^2)*a^3 + (7/2*c*b + c^3)*a^2 + (-1/2*b^2 - c^2*b)*a + 1/4*c*b^2]

   1^2 1 1 1
   if a!= 0,
   [b, c, (4*a^5 - 3*b*a^4 + 2*c*a^3 + d)/a^2, (3*a^5 - 2*b*a^4 + c*a^3 + 2*d)/a, d]
   else
   x^5 + b*x^4 + c*x^3 + d*x^2
   
*/
