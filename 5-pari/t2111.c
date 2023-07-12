#include "c-target.h"

#define pp plist[pnum]
#define aa pdata[pnum][7][0]
#define bb pdata[pnum][7][1]
#define cc pdata[pnum][7][2]
#define dd pdata[pnum][7][3]

/* basic functions:

   init - initialize values given number of known coefs, starting vector of these coefs, p
   next - give the next value after an init has been called, returns 0 if there are no more
*/


/* Initialize a variable counter
    produce next congruence for that counter
*/

int t2111init(long pnum, long which, long what) {
  long j;

  /*
    We have to give the modulus and the value
     Input of which = 0 says to return for a_1
  */
  switch(which) {
  case 0:
    pdata[pnum][which][0] = 1; /* disc exponent */
    pdata[pnum][which][1] = 3; /* Number of prime factors */
    for(j=1; j<=5; j++)
      pflags[pnum][j] = CALL_NONE;
    pflags[pnum][2] = CALL_INIT;
    pflags[pnum][5] = CALL_INIT;
    pflags[pnum][3] = CALL_BOTH;
    pdata[pnum][1][0] = 1; /* modulus */
    pdata[pnum][1][1] = 0; /* starting value */
    pdata[pnum][2][0] = 1;
    pdata[pnum][2][1] = 0;
    pdata[pnum][5][0] = 1;
    pdata[pnum][5][1] = 0;
    
    pdata[pnum][3][0] = pp;
    pdata[pnum][4][0] = pp;
    break;
  case 2: /* Capture b */
    bb = what;  /* Store this for future reference */
    break;
  case 5:
    cc = what;
    break;
  case 3:
    dd = what;
    if(dd==0) {
      aa = 0;
      pdata[pnum][3][0] = 1; /* Sets up the a=0 part */
      pdata[pnum][3][1] = 0;
      pdata[pnum][4][1] = 0;
    } else {
      aa=1;
      pdata[pnum][3][0] = pp;
      pdata[pnum][3][1] = truemod(-3*bb + 2*cc + dd + 4, pp); /* Values for given d and a=1 */
      pdata[pnum][4][1] = truemod(-2*bb + cc + 2*dd + 3, pp);
    }
    break;
  default:
    fprintf(stderr,"Error in t2111 init.\n");
    break;
  }
  return 1;
}

/* increment counter number pnum modulo p
     returns 1 if ok */
long t2111next(long pnum, long which) { /* This should never be called */
  long long a2, a3, a4, ainv, a2inv;
  
  if(which != 3) {
    fprintf(stderr,"Error in t2111 next.\n");
    return 0;
  }
  if(aa==0) {
    pdata[pnum][3][0] = pp;
  }
  
  aa++;
  if(aa==pp)
    return 0;
  
  a2 = (aa*aa) % pp;
  a3 = (aa*a2) % pp;
  a4 = (a2*a2) % pp;
  ainv = minv(aa, pp);
  a2inv = (ainv*ainv) % pp;
  
  pdata[pnum][3][1] = truemod(4*a3 - 3*bb*a2 + 2*cc*aa+dd*a2inv, pp);
  pdata[pnum][4][1] = truemod(3*a4 - 2*bb*a3 + cc*a2 + 2*dd*ainv, pp);
  return 1;
}

/*
   1^2 1 1 1
   if a!= 0,
   [b, c, (4*a^5 - 3*b*a^4 + 2*c*a^3 + d)/a^2, (3*a^5 - 2*b*a^4 + c*a^3 + 2*d)/a, d]
   else
   x^5 + b*x^4 + c*x^3 + d*x^2
   
*/
