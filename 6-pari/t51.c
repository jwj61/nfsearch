#include "c-target.h"

#define pp plist[pnum]
#define bb pdata[pnum][7][0]
#define aa pdata[pnum][7][1]

/* Assume p>3
   i2 = (p+1)/2;
   i3 = ((p%3) == 1) ?  (2*p+1)/3 : (p+1)/3;
*/

/* basic functions:

   init - initialize values given number of known coefs, starting vector of these coefs, p
   next - give the next value after an init has been called, returns 0 if there are no more
   
   do we need to store stuff globally?  There may be multiple instances of this.  Maybe send
   pnum to all functions to tell us which p this is for.
*/


/* Initialize a variable counter
    produce next congruence for that counter
*/

/* Let's START *************************************************************/
int t51init(long pnum, long which, long what) {
  long j, k;

  /*
    We have to give the modulus and the value
     Input of which = 0 says to return for a_1
  */
  switch(which) {
  case 0:
    pdata[pnum][which][0] = 4;
    pdata[pnum][which][1] = 2;
    pflags[pnum][1] = -1;
    pflags[pnum][2] = 0;
    pflags[pnum][6] = -1;
    for(j=3; j<=5; j++)
      pflags[pnum][j] = 1;
    break;
  case 1: /* b */
    pdata[pnum][which][0] = 1; /* modulus */
    pdata[pnum][which][1] = 0; /* starting value */
    break;
  case 2: /* -15*a^2 + 5*b*a */
    pdata[pnum][2][0] = pp;
    pdata[pnum][2][1] = 0;
    bb = what % pp; /* Store this for future reference, it is b */
    aa = 0; /* This is where we will keep a */
    break;
  case 6:
    pdata[pnum][3][0] = pp;
    j = (aa*aa) % pp;
    k = truemod(10*bb-40*aa, pp);    
    pdata[pnum][3][1] = (j*k) % pp;

    pdata[pnum][4][0] = pp;
    j = (aa*aa) % pp;
    j = (j*aa) % pp;
    k = truemod(10*bb-45*aa, pp);
    pdata[pnum][4][1] = (j*k) % pp;
    
    pdata[pnum][5][0] = pp;
    j = (aa*aa) % pp;
    j = (j*j) % pp;
    k = truemod(5*bb-24*aa, pp);
    pdata[pnum][5][1] = (j*k) % pp;
    
    pdata[pnum][6][0] = pp;
    j = (aa*aa) % pp;
    j = (j*j) % pp;
    j = (j*aa) % pp;
    k = truemod(bb-5*aa, pp);
    pdata[pnum][6][1] = (j*k) % pp;
    break;
  default:
    pdata[pnum][which][0] = pp;
    break;
  }
  return 1;
}

/* increment counter number pnum modulo p
     returns 1 if ok */
long t51next(long pnum, long which) {
  long long a, b;
  long p;
  
  if(which==2) {
    p = pp;
    b = pdata[pnum][7][0];
    a = ++(pdata[pnum][7][1]);
    if(a==p)
      return 0; /* Done with this case! */
    /* -15*a^2 + 5*b*a */
    pdata[pnum][which][1] = truemod(5*b*a-15*a*a, p);
    return 1;
  }
  /* in all other cases, there is only one case to consider */
  return 0;
}


/* 1^2 1111
   if a != 0,
 x^6 + b*x^5 + c*x^4 + d*x^3 + ((-5*a^6 + 4*b*a^5 - 3*c*a^4 + 2*d*a^3 + e)/a^2)*x^2 + ((-4*a^6 + 3*b*a^5 - 2*c*a^4 + d*a^3 + 2*e)/a)*x + e

 else x^6 + b*x^5 + c*x^4 + d*x^3 + e*x^2

 so const. term = 0, anything up to linear term.  Then, rest 0 or

  x^6 + b*x^5 + c*x^4 + d*x^3 + (-5*a^4 + 4*b*a^3 - 3*c*a^2 + 2*d*a)*x^2 +
    (-4*a^5 + 3*b*a^4 - 2*c*a^3 + d*a^2)*x

--------------------------------------
51
[b, -15*a^2 + 5*b*a,
  -40*a^3 + 10*b*a^2,
  -45*a^4 + 10*b*a^3,
  -24*a^5 + 5*b*a^4,
  -5*a^6 + b*a^5]

--------------------------------------
42

[b, -6*a^2 + 2*b*a + 1/4*b^2, -4*a^3 - 2*b*a^2 + b^2*a, 9*a^4 - 8*b*a^3 + 3/2*b^2*a^2,
  12*a^5 - 7*b*a^4 + b^2*a^3, 4*a^6 - 2*b*a^5 + 1/4*b^2*a^4]


--------------------------------------
411
    
[b, c, 20*a^3 - 10*b*a^2 + 4*c*a, 45*a^4 - 20*b*a^3 + 6*c*a^2,
  36*a^5 - 15*b*a^4 + 4*c*a^3, 10*a^6 - 4*b*a^5 + c*a^4]


--------------------------------------

33

 [a, b, -5/27*a^3 + 2/3*b*a, -1/9*b*a^2 + 1/3*b^2,
   1/81*a^5 - 2/27*b*a^3 + 1/9*b^2*a,
   -1/729*a^6 + 1/81*b*a^4 - 1/27*b^2*a^2 + 1/27*b^3]
 
 --------------------------------------
321
 
[c, -6*a^2 + (-6*b + 3*c)*a + (-3*b^2 + 2*c*b),
  -8*a^3 + (-18*b + 3*c)*a^2 + (-12*b^2 + 6*c*b)*a + (-2*b^3 + c*b^2),
  -3*a^4 + (-18*b + c)*a^3 + (-18*b^2 + 6*c*b)*a^2 + (-6*b^3 + 3*c*b^2)*a,
  -6*b*a^4 + (-12*b^2 + 2*c*b)*a^3 + (-6*b^3 + 3*c*b^2)*a^2,
  -3*b^2*a^4 + (-2*b^3 + c*b^2)*a^3]

  --------------------------------------
3111

 [b, c, (10*a^6 - 6*b*a^5 + 3*c*a^4 + d)/a^3,
   (15*a^6 - 8*b*a^5 + 3*c*a^4 + 3*d)/a^2,
   (6*a^6 - 3*b*a^5 + c*a^4 + 3*d)/a, d]

 and

 [b, c, d, 0, 0, 0]
 
  --------------------------------------
222

 [a, b, -1/8*a^3 + 1/2*b*a + 2*c, 1/64*a^4 - 1/8*b*a^2 + c*a + 1/4*b^2,
   -1/4*c*a^2 + c*b, c^2]

  --------------------------------------
2211

  [c, d, 4*a^3 - 3*c*a^2 + (-6*b + 2*d)*a + 2*c*b,
    3*a^4 - 2*c*a^3 + d*a^2 - 2*c*b*a + (-3*b^2 + 2*d*b),
    6*b*a^3 - 4*c*b*a^2 + (-6*b^2 + 2*d*b)*a + c*b^2,
    3*b^2*a^2 - 2*c*b^2*a + (-2*b^3 + d*b^2)]

  --------------------------------------
21111

 From a,b,c,d,e

 [b, e, c, d, 0, 0]

 or
 
 [b, e, c,
   (-5*a^6 + 4*b*a^5 - 3*e*a^4 + 2*c*a^3 + d)/a^2,
   (-4*a^6 + 3*b*a^5 - 2*e*a^4 + c*a^3 + 2*d)/a,
   d]
*/
