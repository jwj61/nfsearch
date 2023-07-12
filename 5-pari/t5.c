#include "c-target.h"

#define pp plist[pnum]
#define aa pdata[pnum][7][0]
#define i5 pdata[pnum][7][1]

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

int t5init(long pnum, long which, long what) {
  long j;
  long long aaa;

  /*
    We have to give the modulus and the value
     Input of which = 0 says to return for a_1
  */
  switch(which) {
  case 0:
    pdata[pnum][which][0] = 4;
    pdata[pnum][which][1] = 1;
    i5 = (long) minv(5, pp);    
    for(j=1; j<=5; j++)
      pflags[pnum][j] = CALL_NONE;
    pflags[pnum][2] = CALL_INIT;
    pdata[pnum][1][0] = 1; /* modulus */
    pdata[pnum][1][1] = 0; /* starting value */
    for(j=2; j<=DEG; j++)
      pdata[pnum][j][0] = pp;
    break;
  case 2: /* Take the value of a and use it to initialize the rest */
    aa = (what*i5) % pp;  /* Store this for future reference */
    aaa = (aa*aa) % pp; /* a^2 */
    pdata[pnum][2][1] = (10*aaa) % pp;
    aaa = (aaa*aa) % pp; /* a^3 */
    pdata[pnum][3][1] = (10*aaa) % pp;
    aaa = (aaa*aa) % pp; /* Now a^4 */
    pdata[pnum][4][1] = (5*aaa) % pp;
    pdata[pnum][5][1] = (aa*aaa) % pp;
    
    break;
  default:
    break;
  }
  return 1;
}

/* increment counter number pnum modulo p
     returns 1 if ok */
long t5next(long pnum, long which) { /* This should never be called */
  return 0;
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
